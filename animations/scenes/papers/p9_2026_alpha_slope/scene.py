"""Eldadi & Loeb (2026) -- a framework for applying the Loeb-Turner alpha-slope
test to archival photometry of trans-Neptunian objects.

Flux against heliocentric distance follows a power law, F ~ d^alpha: -4 for
reflected sunlight, -2 for a body that shines by itself. The crate fits alpha
by regression; here it recovers both slopes from the p9-core flux laws, and
shows how little calibration drift it takes to corrupt the fit when a body
barely changes distance. The census counts are the paper's.
Everything is from anim.json -> papers -> p9-2026-alpha-slope.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Polygon,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2026-alpha-slope"


def rows(items, font_size=15, buff=0.14):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


def framed_axes(x_range, y_range, x_label, y_label, x_ticks, y_ticks, x_length=7.4,
                y_length=4.0):
    """Axes with hand-placed tick labels (for logarithmic scales)."""
    ax = widgets.axes(x_range, y_range, x_length=x_length, y_length=y_length, shift_down=-0.2)
    marks = VGroup()
    for x, text in x_ticks:
        marks.add(layout.label(text, font_size=13, color=P.MUTED)
                  .next_to(ax.c2p(x, y_range[0]), DOWN, buff=0.12))
    for y, text in y_ticks:
        marks.add(layout.label(text, font_size=13, color=P.MUTED)
                  .next_to(ax.c2p(x_range[0], y), LEFT, buff=0.12))
    xl = layout.label(x_label, font_size=17, color=P.FG)
    xl.next_to(ax, DOWN, buff=0.42)
    yl = layout.label(y_label, font_size=15, color=P.FG).rotate(np.pi / 2)
    yl.next_to(ax, LEFT, buff=0.6)
    g = VGroup(ax, marks, xl, yl)
    g.shift(LEFT * 2.0)
    return ax, g


def band(ax, x0, x1, y0, y1, color, opacity=0.18):
    return Polygon(ax.c2p(x0, y0), ax.c2p(x1, y0), ax.c2p(x1, y1), ax.c2p(x0, y1),
                   stroke_width=0).set_fill(color, opacity=opacity)


class AlphaSlope2026(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        th = d["theory"]
        pl = d["pluto"]
        fn = d["funnel"]

        self.add(paper.scene_header(CRATE))

        # 1. two laws, two slopes
        laws = d["laws"]
        # Axes cross at their origin, so both panels plot offsets from the lower-left
        # corner and label the ticks with the true values.
        lx0, ly0 = np.log10(28), -2.6
        lx = np.log10(laws["distance_au"]) - lx0
        ticks = [30, 40, 60, 80, 120]
        ax, furniture = framed_axes(
            [0, np.log10(125) - lx0, 10], [0, 0.2 - ly0, 10],
            "heliocentric distance (AU)", "brightness, relative to 30 AU",
            [(np.log10(t) - lx0, str(t)) for t in ticks],
            [(0 - ly0, "1"), (-1 - ly0, "1/10"), (-2 - ly0, "1/100")])
        refl = widgets.curve(ax, lx, np.log10(laws["reflected"]) - ly0, color=P.GREEN)
        self_l = widgets.curve(ax, lx, np.log10(laws["thermal"]) - ly0, color=P.ORANGE)
        key = rows([
            ("reflected sunlight", P.GREEN),
            ("out to the body, and back", P.GREEN),
            (f"fitted slope {d['alpha_reflected']:.2f}", P.GREEN),
            ("shining by itself", P.ORANGE),
            ("one way only", P.ORANGE),
            (f"fitted slope {d['alpha_thermal']:.2f}", P.ORANGE),
        ])
        key[3:].shift(DOWN * 0.25)
        key.move_to([4.7, 1.3, 0])
        cap = layout.caption(
            "Follow a body as its distance changes: how fast it fades tells how it shines",
            font_size=21)
        self.play(FadeIn(furniture), FadeIn(cap), run_time=0.9)
        self.play(Create(refl), FadeIn(key[:3]), run_time=1.1)
        self.play(Create(self_l), FadeIn(key[3:]), run_time=1.1)
        timing.hold_to_read(self, cap, key, settle=0.7)
        self.play(FadeOut(VGroup(furniture, refl, self_l, key, cap)))

        # 2. the lever arm is short, so calibration decides the slope
        drift = np.array(pl["drift_mag"])
        alpha = np.array(pl["alpha"])
        keep = (alpha >= -7) & (alpha <= -1)
        ax2, labels2 = framed_axes(
            [0, 2, 10], [0, 6, 10],
            "calibration drift across the record (magnitudes)", "fitted slope",
            [(x + 1, f"{x:+.1f}" if x else "0") for x in (-1, -0.5, 0, 0.5, 1)],
            [(y + 7, f"{y}") for y in range(-7, 0)])

        def q(x, y):
            return ax2.c2p(x + 1, y + 7)

        tol = th["tolerance"]
        win_r = band(ax2, 0, 2, th["reflected"] - tol + 7, th["reflected"] + tol + 7, P.GREEN)
        win_s = band(ax2, 0, 2, th["self_luminous"] - tol + 7, th["self_luminous"] + tol + 7,
                     P.ORANGE)
        win_r_lab = layout.label("reads as reflected", font_size=13, color=P.GREEN)
        win_r_lab.move_to(q(0.62, th["reflected"] + 0.15))
        win_s_lab = layout.label("reads as self-luminous", font_size=13, color=P.ORANGE)
        win_s_lab.move_to(q(0.55, th["self_luminous"]))
        line = widgets.curve(ax2, drift[keep] + 1, alpha[keep] + 7, color=P.FG)
        true = Dot(q(0, pl["alpha_clean"]), radius=0.08, color=P.GREEN)
        flip = Dot(q(pl["drift_to_flip_mag"], th["self_luminous"] - tol), radius=0.08,
                   color=P.ORANGE)
        drop = DashedLine(q(pl["drift_to_flip_mag"], -7),
                          q(pl["drift_to_flip_mag"], th["self_luminous"] - tol),
                          color=P.ORANGE, stroke_width=2)
        zero = DashedLine(q(0, -7), q(0, -1), color=P.MUTED, stroke_width=1.2)
        lever = d["lever"]
        key2 = rows([
            (f"Pluto, {pl['range_au'][0]:.1f} to {pl['range_au'][1]:.0f} AU", P.FG),
            (f"0.1 mag of drift moves the slope "
             f"{abs(pl['alpha_per_tenth_mag']):.2f}", P.FG),
            (f"{abs(pl['drift_to_flip_mag']):.2f} mag turns reflected", P.ORANGE),
            ("into self-luminous", P.ORANGE),
            (f"a body whose distance changes {100 * (lever['span'][0] - 1):.0f}%:", P.FG),
            (f"0.1 mag moves the slope "
             f"{lever['alpha_error_per_tenth_mag'][0]:.1f}", P.RED),
        ], font_size=14)
        key2[2:].shift(DOWN * 0.2)
        key2[4:].shift(DOWN * 0.2)
        key2.move_to([4.75, 1.1, 0])
        cap2 = layout.caption(
            "Distant bodies barely change distance, so a small calibration drift "
            "rewrites the slope", font_size=21)
        order = np.argsort(drift[keep])
        dx, ay = drift[keep][order], alpha[keep][order]
        sweep = ValueTracker(0.0)
        probe = always_redraw(lambda: Dot(
            q(sweep.get_value(), float(np.interp(sweep.get_value(), dx, ay))),
            radius=0.07, color=P.FG))
        self.play(FadeIn(labels2), run_time=0.9)
        self.play(FadeIn(win_r), FadeIn(win_s), FadeIn(win_r_lab), FadeIn(win_s_lab),
                  Create(zero), Create(line), FadeIn(true), FadeIn(key2[:2]), FadeIn(cap2),
                  run_time=1.4)
        self.add(probe)
        self.play(sweep.animate.set_value(pl["drift_to_flip_mag"]), run_time=2.0)
        self.play(Create(drop), FadeIn(flip), FadeIn(key2[2:]), run_time=1.0)
        self.remove(probe)
        timing.hold_to_read(self, cap2, key2, settle=0.8)
        self.play(FadeOut(VGroup(labels2, win_r, win_s, win_r_lab, win_s_lab, zero, line, true,
                                 flip, drop, key2, cap2)))

        # 3. the census of the archive (the paper's counts)
        x0, full = -5.6, 8.2
        stages = [
            (fn["candidate_bins"], "bins: object × observatory × band", P.MUTED),
            (fn["pass_q1_q3"], "enough epochs and distance range", P.MUTED),
            (fn["pass_q4_q6"], "a clean, decisive fit", P.FG),
        ]
        bars = VGroup()
        for k, (n, text, col) in enumerate(stages):
            w = full * np.log10(n) / np.log10(stages[0][0])
            y = 2.2 - 0.95 * k
            bar = Rectangle(width=w, height=0.5, stroke_width=0).set_fill(col, opacity=0.35)
            bar.move_to([x0 + w / 2, y, 0])
            num = layout.label(f"{n:,}", font_size=18, color=P.FG, weight="BOLD")
            num.next_to(bar, RIGHT, buff=0.15)
            lab = layout.label(text, font_size=14, color=P.FG).next_to(num, RIGHT, buff=0.2)
            bars.add(VGroup(bar, num, lab))
        scale_note = layout.label("bar length: logarithmic", font_size=12, color=P.MUTED)
        scale_note.move_to([x0 + 1.0, 2.75, 0])
        parts = [
            (fn["reflected"], "reflected", P.GREEN),
            (fn["self_luminous"], "self-luminous", P.ORANGE),
            (fn["anomalous"], "neither", P.MUTED),
        ]
        split = VGroup()
        cursor = x0
        wide = 11.2
        for n, text, col in parts:
            w = wide * n / fn["pass_q4_q6"]
            seg = Rectangle(width=w, height=0.6, stroke_width=1, color=P.BG)
            seg.set_fill(col, opacity=0.75).move_to([cursor + w / 2, -0.82, 0])
            tag = layout.label(f"{n}  {text}", font_size=14, color=col, weight="BOLD")
            tag.next_to(seg, DOWN, buff=0.12)
            split.add(VGroup(seg, tag))
            cursor += w
        split_head = layout.label(f"the {fn['pass_q4_q6']} clean fits, by slope", font_size=14,
                                  color=P.FG)
        split_head.next_to(split[0][0], UP, buff=0.12, aligned_edge=LEFT)
        pluto = layout.label(
            f"Pluto: {pl['bins']} bins, {pl['recoveries']} recover the reflected slope",
            font_size=16, color=P.RED, weight="BOLD")
        pluto.move_to([0, -2.32, 0])
        origin = layout.label(
            f"all {fn['self_luminous']} self-luminous fits come from Pan-STARRS: "
            "calibration, not physics", font_size=16, color=P.ORANGE, weight="BOLD")
        origin.move_to([0, -1.9, 0])
        cap3 = layout.caption(
            f"The archive, {fn['tnos']} numbered trans-Neptunian objects, put through "
            "six quality cuts", font_size=21)
        self.play(FadeIn(bars, lag_ratio=0.3), FadeIn(scale_note), FadeIn(cap3), run_time=1.8)
        self.play(FadeIn(split_head), FadeIn(split, lag_ratio=0.3), run_time=1.5)
        timing.hold_to_read(self, cap3, split, settle=0.6)
        self.play(FadeIn(origin))
        timing.hold_to_read(self, origin, settle=0.5)
        self.play(FadeIn(pluto))
        timing.hold_to_read(self, pluto, settle=0.9)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, "The archive cannot run the test cleanly, even on Pluto. Rubin can.")
