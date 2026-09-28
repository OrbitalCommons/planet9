"""Geringer-Sameth, Golovich & Iwabuchi (2025) -- multi-year stacking searches
for solar system bodies.

Digital tracking: add every image of a field along a trial orbit, so a body too
faint for any one exposure adds up. The paper's contribution is the metric that
says how finely the trial orbits must be spaced, hence how many are needed, and
what the search then reaches. The crate computes the stacking gain, the signal
lost to a mismatched orbit, the trial count and its look-elsewhere cost, and
the share of the predicted Planet Nines a ZTF stack would newly reach.
Everything is from anim.json -> papers -> p9-2025-stacking.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    DashedVMobject,
    Dot,
    FadeIn,
    FadeOut,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-stacking"


def rows(items, font_size=15, buff=0.14):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


SUPERSCRIPT = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def power_label(k):
    return {0: "1", 1: "10", 2: "100", 3: "1,000"}.get(k, "10" + str(k).translate(SUPERSCRIPT))


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
    yl.next_to(ax, LEFT, buff=0.55)
    g = VGroup(ax, marks, xl, yl)
    g.shift(LEFT * 2.0)
    return ax, g


class Stacking2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        mm = d["mismatch"]
        p9 = d["p9"]

        self.add(paper.scene_header(CRATE))

        # 1. depth bought by stacking
        single, rubin = d["single_depth"], d["rubin_single_depth"]
        ax, furniture = framed_axes(
            [0, 6, 1], [20, 28, 1], "images stacked along the right orbit",
            "limiting magnitude  (fainter ↑)",
            [(k, power_label(k)) for k in range(0, 7)],
            [(m, str(m)) for m in range(20, 29, 2)])
        line = widgets.curve(ax, np.log10(d["depth"]["frames"]), d["depth"]["mag"],
                             color=P.TEAL)
        n_stack = d["stack_frames"]
        here = Dot(ax.c2p(np.log10(n_stack), d["stacked_depth"]), radius=0.08, color=P.TEAL)
        drop = DashedLine(ax.c2p(np.log10(n_stack), 20), ax.c2p(np.log10(n_stack),
                                                                 d["stacked_depth"]),
                          color=P.TEAL, stroke_width=2)
        rubin_line = DashedLine(ax.c2p(0, rubin), ax.c2p(6, rubin), color=P.PURPLE,
                                stroke_width=2)
        rubin_lab = layout.label(f"one Rubin exposure  {rubin:.1f}", font_size=13,
                                 color=P.PURPLE)
        rubin_lab.next_to(ax.c2p(6, rubin), UP, buff=0.06).shift(LEFT * 1.1)
        facts = rows([
            (f"one ZTF exposure: {single:.1f}", P.FG),
            (f"{n_stack:,} exposures: {d['stacked_depth']:.1f}", P.TEAL),
            (f"a gain of {d['gain_mag']:.1f} magnitudes", P.TEAL),
            (f"magnitude {pub['stacked_depth']:.0f} needs "
             f"{d['frames_to_27'] / 1000:.0f},000", P.FG),
        ])
        facts.move_to([4.7, 1.5, 0])
        cap = layout.caption("Depth grows by 1.25 magnitudes for every tenfold more images",
                             font_size=21)
        self.play(FadeIn(furniture), run_time=0.9)
        self.play(Create(line), FadeIn(cap), FadeIn(facts[0]), run_time=1.2)
        self.play(Create(rubin_line), FadeIn(rubin_lab), Create(drop), FadeIn(here),
                  FadeIn(facts[1:]), run_time=1.2)
        timing.hold_to_read(self, cap, facts, settle=0.7)
        self.play(FadeOut(VGroup(furniture, line, here, drop, rubin_line, rubin_lab, facts,
                                 cap)))

        # 2. but only along the right orbit
        err = 1000.0 * np.array(mm["rate_error"])
        top = float(np.ceil(err[-1]))
        ax2, labels2 = widgets.labeled_axes(
            [0, top, 1], [0, 1, 0.2], x_label="error in the trial sky rate (milliarcsec per day)",
            y_label="fraction of the signal kept", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.0, shift_down=-0.2)
        VGroup(ax2, labels2).shift(LEFT * 2.0)
        exact = widgets.curve(ax2, err, mm["retained_exact"], color=P.ORANGE)
        quad = DashedVMobject(
            widgets.curve(ax2, err, mm["retained_quadratic"], color=P.MUTED, stroke_width=2),
            num_dashes=40)
        cell = mm["cell_mas_day"]
        tol = mm["tolerance"]
        mark = VGroup(
            DashedLine(ax2.c2p(0, tol), ax2.c2p(cell, tol), color=P.PURPLE, stroke_width=2),
            DashedLine(ax2.c2p(cell, 0), ax2.c2p(cell, tol), color=P.PURPLE, stroke_width=2),
            Dot(ax2.c2p(cell, tol), radius=0.07, color=P.PURPLE))
        key2 = rows([
            ("stacked signal, exact", P.ORANGE),
            ("the metric's quadratic form", P.MUTED),
            (f"keep {100 * tol:.0f}%: trial orbits", P.PURPLE),
            (f"every {cell:.1f} mas per day", P.PURPLE),
        ])
        key2.move_to([4.7, 1.5, 0])
        cap2 = layout.caption(
            f"Over {d['baseline_years']:.0f} years a tiny error in the assumed motion smears "
            f"the {mm['psf_arcsec']:.0f}″ image away", font_size=21)
        self.play(Create(ax2), FadeIn(labels2), run_time=0.9)
        self.play(Create(exact), Create(quad), FadeIn(key2[:2]), FadeIn(cap2), run_time=1.3)
        self.play(Create(mark), FadeIn(key2[2:]))
        timing.hold_to_read(self, cap2, key2, settle=0.7)
        self.play(FadeOut(VGroup(ax2, labels2, exact, quad, mark, key2, cap2)))

        # 3. so count the orbits, and pay for them
        tr = d["trials"]
        # offset by half a decade so the axes cross at the lower-left corner
        log_t = np.log10(tr["baseline_days"]) + 0.5
        log_n = np.log10(tr["n"])
        ax3, furniture3 = framed_axes(
            [0, 4, 10], [0, 10, 2], "time spanned by the images",
            "trial orbits",
            [(0.5, "1 day"), (1.5, "10 days"), (2.5, "100 days"), (3.5, "1,000 days")],
            [(k, power_label(k)) for k in range(0, 11, 2)])
        grow = widgets.curve(ax3, log_t, log_n, color=P.ORANGE)
        end = Dot(ax3.c2p(log_t[-1], log_n[-1]), radius=0.08, color=P.ORANGE)
        cost = rows([
            (f"{d['baseline_years']:.0f} years of ZTF", P.FG),
            (f"{d['n_trials'] / 1e6:.0f} million trial orbits", P.ORANGE),
            ("paper: billions", P.MUTED),
            (f"threshold rises {d['baseline_sigma']:.0f}σ to {d['threshold_sigma']:.1f}σ",
             P.RED),
            (f"costing {d['penalty_mag']:.1f} magnitude", P.RED),
            (f"stacking gains {d['gain_mag']:.1f}", P.TEAL),
            (f"net depth {d['net_depth']:.1f}", P.TEAL),
        ])
        cost[3:].shift(DOWN * 0.2)
        cost[5:].shift(DOWN * 0.2)
        cost.move_to([4.7, 0.9, 0])
        cap3 = layout.caption(
            "The number of orbits to try grows as the square of the time spanned",
            font_size=21)
        self.play(FadeIn(furniture3), run_time=0.9)
        self.play(Create(grow), FadeIn(end), FadeIn(cost[:3]), FadeIn(cap3), run_time=1.3)
        timing.hold_to_read(self, cap3, cost[:3], settle=0.4)
        cap3b = layout.caption(
            "More trials mean more chance alarms, but the price is logarithmic",
            font_size=21)
        self.play(FadeIn(cost[3:], lag_ratio=0.2), FadeOut(cap3), FadeIn(cap3b))
        timing.hold_to_read(self, cap3b, cost[3:], settle=0.8)
        self.play(FadeOut(VGroup(furniture3, grow, end, cost, cap3b)))

        # 4. what it would do for Planet Nine
        v = np.array(p9["v_mag"])
        found = np.array(p9["p_found"])
        reach = np.array(p9["in_reach"], dtype=bool)
        edges = np.arange(17.0, 26.01, 0.5)
        h_old, _ = np.histogram(v, bins=edges, weights=found)
        h_new, _ = np.histogram(v, bins=edges, weights=np.where(reach, 1.0 - found, 0.0))
        tot, _ = np.histogram(v, bins=edges)
        peak = int(np.ceil(tot.max() / 100.0) * 100)
        ax4, labels4 = widgets.labeled_axes(
            [17, 26, 1], [0, peak, peak // 4], x_label="apparent magnitude V  (fainter →)",
            y_label="synthetic planets", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.0, shift_down=-0.2)
        VGroup(ax4, labels4).shift(LEFT * 2.0)
        bars_old = widgets.histogram(ax4, edges, h_old, color=P.MUTED, opacity=0.6)
        bars_new = widgets.histogram(ax4, edges, h_new, color=P.TEAL, opacity=0.8, base=h_old)
        bars_left = widgets.histogram(ax4, edges, tot - h_old - h_new, color=P.BLUE,
                                      opacity=0.6, base=h_old + h_new)
        mark4 = widgets.marker_line(ax4, p9["depth"], (0, peak),
                                    f"ZTF stack at {p9['sigma']:.0f}σ  V = {p9['depth']:.1f}",
                                    color=P.TEAL, side=UP)
        key4 = rows([
            (f"already ruled out  {100 * p9['already']:.0f}%", P.MUTED),
            (f"paper: {100 * pub['already']:.0f}%", P.MUTED),
            (f"a ZTF stack would reach  {100 * p9['new']:.0f}%", P.TEAL),
            (f"paper: {100 * pub['new']:.0f}%", P.MUTED),
            (f"beyond it  {100 * p9['remaining']:.0f}%", P.BLUE),
        ])
        key4[2:].shift(DOWN * 0.2)
        key4[4:].shift(DOWN * 0.2)
        key4.move_to([4.7, 1.2, 0])
        cap4 = layout.caption(
            f"Teal: not yet searched, but north of {p9['dec_limit_deg']:.0f}° and brighter "
            f"than V = {p9['depth']:.1f}", font_size=21)
        self.play(Create(ax4), FadeIn(labels4), run_time=0.9)
        self.play(FadeIn(bars_old, lag_ratio=0.1), FadeIn(bars_left, lag_ratio=0.1),
                  FadeIn(key4[:2]), FadeIn(key4[4]), run_time=1.2)
        self.play(FadeIn(bars_new, lag_ratio=0.1), Create(mark4), FadeIn(key4[2:4]),
                  FadeIn(cap4), run_time=1.3)
        timing.hold_to_read(self, cap4, key4, settle=0.9)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, "Stacking six years of ZTF reaches most of the Planet Nines still in hiding.")
