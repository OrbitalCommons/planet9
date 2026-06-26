"""Eldadi & Loeb (2026) -- the Loeb-Turner alpha-slope test on TNO photometry.

Reflected sunlight falls off as d^-4 (out and back); self-luminous emission only
as d^-2. The slope alpha of log(flux) vs log(distance) tells them apart. Run
across the MPC archive through a six-criterion pipeline, the "self-luminous"
detections all trace to one instrument -- a calibration tell, not new physics.
Reproduced in p9-2026-alpha-slope.
"""
import numpy as np
from manim import (
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    LEFT,
    RIGHT,
    Scene,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2026-alpha-slope"


class AlphaSlope2026(Scene):
    def construct(self):
        tb = paper.title_block(
            CRATE, "Reflected or self-luminous?", cite="Eldadi & Loeb (2026)")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # --- log(flux) vs log(distance): two diagnostic slopes -----------------
        ax = widgets.axes([0.0, 2.0, 0.5], [-9.0, 0.5, 2.0],
                          x_length=7.6, y_length=3.9, font_size=16, shift_down=0.4)
        xlab = layout.label("log heliocentric distance", font_size=18, color=P.FG)
        xlab.next_to(ax, DOWN, buff=0.25)
        ylab = layout.label("log flux", font_size=18, color=P.FG)
        ylab.rotate(np.pi / 2).next_to(ax, LEFT, buff=0.15)
        self.play(Create(ax), FadeIn(xlab), FadeIn(ylab))

        xs = np.array([0.1, 1.9])
        # Both lines share a left anchor so the fan-out reads as "same body".
        anchor = -0.4
        reflected = anchor - 4.0 * xs   # alpha = -4
        selflum = anchor - 2.0 * xs     # alpha = -2

        line_ref = ax.plot_line_graph(xs, reflected, line_color=P.TEAL,
                                      add_vertex_dots=False, stroke_width=4)
        line_self = ax.plot_line_graph(xs, selflum, line_color=P.RED,
                                       add_vertex_dots=False, stroke_width=4)
        self.play(Create(line_self))
        l_self = layout.label("alpha = -2  self-luminous", font_size=17, color=P.RED)
        l_self.next_to(ax.c2p(1.9, selflum[1]), RIGHT, buff=0.12)
        self.play(FadeIn(l_self))
        self.play(Create(line_ref))
        l_ref = layout.label("alpha = -4  reflected sunlight", font_size=17, color=P.TEAL)
        l_ref.next_to(ax.c2p(1.9, reflected[1]), RIGHT, buff=0.12)
        self.play(FadeIn(l_ref))
        timing.hold_to_read(self, l_ref, l_self, settle=0.7)
        self.play(FadeOut(line_ref), FadeOut(line_self), FadeOut(l_ref),
                  FadeOut(l_self), FadeOut(ax), FadeOut(xlab), FadeOut(ylab))

        eq = layout.explain_equation(
            self,
            [r"F", r"\propto", r"d^{\alpha}", r"\quad",
             r"\alpha=-4\ (\text{reflected}),\ -2\ (\text{self-lum.})"],
            [
                (0, "the measured flux of the body"),
                (2, "falls off as distance to the power alpha"),
                (2, "reflected light loses d-squared out and d-squared back"),
                (4, "thermal/internal light only loses d-squared -- the slope tells them apart"),
            ],
            color=P.FG, scale=0.7, where=UP * 0.4)
        self.play(FadeOut(eq))

        # --- the eligibility funnel -------------------------------------------
        funnel = widgets.value_rows(
            [
                "8,557  candidate (KBO x observatory x band) bins",
                "1,089  pass Q1-Q3",
                "186  also pass Q4-Q6",
            ],
            color=P.FG, font_size=22, line_buff=0.30)
        funnel.move_to(UP * 0.9)
        for row in funnel:
            self.play(FadeIn(row, shift=DOWN * 0.1), run_time=0.5)
        timing.hold_to_read(self, funnel, settle=0.6)

        split = widgets.value_rows(
            [
                "53  reflected   (alpha = -4)",
                "24  self-luminous   (alpha = -2)",
                "109  anomalous",
            ],
            color=P.FG, font_size=22, line_buff=0.26)
        split[0].set_color(P.TEAL)
        split[1].set_color(P.RED)
        split[2].set_color(P.ORANGE)
        split.next_to(funnel, DOWN, buff=0.55)
        self.play(FadeIn(split, lag_ratio=0.3))
        timing.hold_to_read(self, split, settle=0.8)
        self.play(FadeOut(funnel), FadeOut(split))

        # --- the punchline: the 24 all came from one instrument ---------------
        ro = paper.result_readout(
            "self-luminous detections -- all from one instrument",
            "Pan-STARRS only", color=P.RED).scale(0.8)
        ro.move_to(UP * 0.7)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro, settle=0.8)

        pluto = layout.label(
            "Even Pluto (22 bins) cannot cleanly recover alpha = -4: 0 bins do.",
            font_size=20, color=P.MUTED)
        pluto.next_to(ro, DOWN, buff=0.55)
        self.play(FadeIn(pluto))
        timing.hold_to_read(self, pluto, settle=0.8)
        self.play(FadeOut(ro), FadeOut(pluto))

        layout.show_takeaway(
            self,
            "One instrument can't separate the slopes -- "
            "the alpha-test needs cross-calibrated data.")
