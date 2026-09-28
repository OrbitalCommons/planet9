"""Huang & Gladman (2023) -- primordial orbital alignment of sednoids.

The three sednoids are too far out for Neptune to kick, so the only thing that
has happened to their orbits is a slow, steady turning of the apsidal line
forced by the four giant planets, faster for nearer orbits. Run that clock
backwards and today's modest bunching becomes a tight alignment once, about
4.5 Gyr ago. The scene rewinds the crate's precession tracks for Sedna,
2012 VP113 and 2015 TG387 while the alignment statistic is traced, then shows
how special that epoch is and compares a less-detached sample.
Data: anim.json -> papers -> p9-2024-primordial-alignment.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arrow,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    smooth,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2024-primordial-alignment"

SHORT = {
    "Sedna": "Sedna",
    "Alicanto (2012 VP113)": "2012 VP113",
    "Leleakuhonua (2015 TG387)": "2015 TG387",
}


def unit(deg):
    r = np.radians(deg)
    return np.array([np.cos(r), np.sin(r), 0.0])


def interp_angle(t, ts, degs):
    """Angle track sampled at ``ts`` (unwrapped), evaluated at ``t``."""
    un = np.degrees(np.unwrap(np.radians(degs)))
    return float(np.interp(t, ts, un)) % 360.0


class PrimordialAlignment2024(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        ts = np.array(d["tau_gyr"])
        rbar = np.array(d["r_bar"])
        objs = d["objects"]
        best = d["best_epoch_gyr"]
        pub = d["published_epoch_gyr"]

        self.add(paper.scene_header(CRATE))

        # the clock: perihelion directions as hands on a dial
        c = np.array([-4.1, 0.3, 0.0])
        R = 2.2
        dial = Circle(radius=R, color=P.MUTED, stroke_width=1.3).move_to(c)
        marks = VGroup(*[layout.label(f"{x}°", font_size=13, color=P.MUTED)
                         .move_to(c + unit(x) * (R + 0.3)) for x in (0, 90, 180, 270)])
        sun = orbits.sun(radius=0.08).move_to(c)
        lengths = [1.5, 0.95, 1.95]
        # the shortest hand's name sits beside its tip, clear of the longer hands
        side = [0.0, 0.8, 0.0]
        tau = ValueTracker(0.0)

        def hand(k):
            o = objs[k]
            ang = interp_angle(tau.get_value(), ts, o["varpi_track_deg"])
            tip = c + unit(ang) * lengths[k]
            arr = Arrow(c, tip, buff=0, color=P.GREEN, stroke_width=3.2,
                        max_tip_length_to_length_ratio=0.1)
            lab = layout.label(SHORT[o["name"]], font_size=14, color=P.GREEN)
            lab.move_to(tip + unit(ang) * 0.12 + unit(ang - 90) * side[k])
            lab.shift(unit(ang) * (lab.width / 2 * abs(np.cos(np.radians(ang)))
                                   + lab.height / 2 * abs(np.sin(np.radians(ang)))))
            return VGroup(arr, lab)

        hands = VGroup(*[always_redraw(lambda k=k: hand(k)) for k in range(len(objs))])
        clock = always_redraw(lambda: layout.label(
            "today" if tau.get_value() < 0.005 else f"{tau.get_value():.2f} billion years ago",
            font_size=18, color=P.FG).move_to(c + DOWN * (R + 0.55)))

        rates = VGroup(
            layout.label("one full turn of the apse takes", font_size=16, color=P.MUTED),
            *[layout.label(f"{SHORT[o['name']]}  (a = {o['a_au']:.0f} AU):  "
                           f"{o['period_gyr']:.1f} Gyr", font_size=16, color=P.GREEN)
              for o in sorted(objs, key=lambda o: o["a_au"])],
            layout.label("nearer orbits turn faster", font_size=16, color=P.ORANGE),
        ).arrange(DOWN, buff=0.17, aligned_edge=LEFT).move_to([2.9, 0.9, 0])
        cap = layout.caption("Three sednoids today: arrows point to each perihelion", font_size=22)
        self.play(FadeIn(dial), FadeIn(marks), FadeIn(sun), FadeIn(hands), FadeIn(clock),
                  FadeIn(cap), run_time=1.0)
        timing.hold_to_read(self, cap, settle=0.6)
        cap2 = layout.caption("Beyond Neptune's reach, only the giant planets' slow turning acts",
                              font_size=22)
        self.play(FadeIn(rates), FadeOut(cap), FadeIn(cap2), run_time=0.8)
        timing.hold_to_read(self, cap2, rates, settle=0.8)
        self.play(FadeOut(rates), run_time=0.5)

        # rewind: hands turn back while the alignment is traced
        ax, labs = widgets.labeled_axes(
            [0, 5, 1], [0, 1, 0.25], x_label="billions of years ago",
            y_label="alignment  R̄", y_rotate=True, numbers=True,
            x_length=6.0, y_length=3.7, shift_down=0)
        ax.move_to([2.95, 0.45, 0])
        labs[0].next_to(ax, DOWN, buff=0.5)
        labs[1].next_to(ax, LEFT, buff=0.55)
        key = VGroup(
            layout.label("R̄ = 1: all three point the same way", font_size=14, color=P.MUTED),
            layout.label("R̄ = 0: evenly spread", font_size=14, color=P.MUTED),
        ).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        key.next_to(ax.c2p(0, 1), UP + RIGHT, buff=0.1).shift(UP * 0.05)

        def traced():
            t = tau.get_value()
            m = ts <= t + 1e-9
            if m.sum() < 2:
                return VGroup()
            return widgets.curve(ax, ts[m], rbar[m], color=P.GREEN)

        trace = always_redraw(traced)
        rider = always_redraw(lambda: Dot(
            ax.c2p(tau.get_value(), float(np.interp(tau.get_value(), ts, rbar))),
            radius=0.07, color=P.GREEN))
        cap3 = layout.caption("Run the clock backwards", font_size=22)
        self.play(FadeIn(ax), FadeIn(labs), FadeIn(key), FadeIn(rider), FadeOut(cap2),
                  FadeIn(cap3), run_time=0.8)
        self.add(trace)
        self.play(tau.animate.set_value(best), run_time=9.0, rate_func=smooth)
        pub_dir = DashedLine(c, c + unit(d["published_varpi_deg"]) * R, color=P.FG,
                             stroke_width=1.6).set_stroke(opacity=0.7)
        cap4 = layout.caption(
            f"{best:.2f} Gyr ago all three point to {d['best_varpi_deg']:.0f}°, R̄ = "
            f"{d['best_r_bar']:.2f}   (paper: {d['published_varpi_deg']:.0f}°, dashed)",
            font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), Create(pub_dir), run_time=0.8)
        timing.hold_to_read(self, cap4, settle=1.0)

        # how special that moment is
        self.play(tau.animate.set_value(ts[-1]), run_time=1.6)
        self.remove(trace)
        full = widgets.curve(ax, ts, rbar, color=P.GREEN)
        self.add(full)
        thr = d["tight_r_bar"]
        tight = DashedLine(ax.c2p(0, thr), ax.c2p(5, thr), color=P.ORANGE, stroke_width=1.8)
        tight_lab = layout.label(f"tight: R̄ ≥ {thr:.2f}", font_size=14, color=P.ORANGE)
        tight_lab.next_to(ax.c2p(2.6, thr), UP, buff=0.06)
        pub_line = DashedLine(ax.c2p(pub, 0), ax.c2p(pub, 1), color=P.FG, stroke_width=1.6)
        pub_tag = layout.label(f"paper: {pub:.1f} Gyr", font_size=14, color=P.FG)
        pub_tag.next_to(ax.c2p(pub, 0.08), LEFT, buff=0.1)
        cap5 = layout.caption(
            f"Tight alignment only {100 * d['fraction_tight']:.0f}% of the time, "
            "and only at the Solar System's birth", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), Create(tight), FadeIn(tight_lab), Create(pub_line),
                  FadeIn(pub_tag), tau.animate.set_value(best), run_time=1.2)
        readout = paper.result_readout("tightest alignment", f"{best:.2f} Gyr ago",
                                       color=P.GREEN).scale(0.6)
        readout.move_to(key, aligned_edge=LEFT)
        self.play(FadeOut(key), FadeIn(readout), run_time=0.5)
        timing.hold_to_read(self, cap5, settle=1.0)

        # contrast: a less-detached sample never re-aligns
        b = d["brown2017"]
        other = widgets.curve(ax, ts, b["r_bar"], color=P.MUTED, stroke_width=2.4)
        other_lab = layout.label(f"the {b['n']} distant objects of Brown (2017)", font_size=14,
                                 color=P.MUTED)
        other_lab.move_to(ax.c2p(2.55, 0.8))
        cap6 = layout.caption("Rewinding Brown's ten, most of which feel Neptune, finds no birth alignment",
                              font_size=22)
        self.play(FadeOut(cap5), FadeIn(cap6), Create(other), FadeIn(other_lab), run_time=1.2)
        timing.hold_to_read(self, cap6, settle=1.0)
        self.play(FadeOut(cap6), run_time=0.4)

        layout.show_takeaway(
            self, "Imprinted at birth and only turned since: a planet still out there would erase it.")
