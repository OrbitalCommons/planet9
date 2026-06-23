"""Lai (2016) -- analytic theory of the solar obliquity from Planet Nine.

An analytic resonance-overlap / secular treatment of how P9 tilts the Sun,
giving the obliquity as a function of P9's mass, distance, and inclination.
Reproduced in p9-2016-obliquity-lai.
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2016-obliquity-lai"


class ObliquityLai2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The obliquity, from first principles")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 30, 10], y_range=[0, 12, 3], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        xl = Text("Planet Nine inclination (deg)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("induced solar obliquity (deg)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        curve = ax.plot(lambda i: 6.0 * np.sin(np.deg2rad(2 * i)) / np.sin(np.deg2rad(40)), x_range=[0, 30, 0.5], color=P.ORANGE)
        self.play(Create(curve), run_time=1.4)
        target = ax.get_horizontal_line(ax.c2p(30, 6.0), color=P.TEAL, stroke_width=2)
        obs = Text("observed 6°", font_size=15, color=P.TEAL).next_to(ax.c2p(2, 6), UP, buff=0.05)
        self.play(Create(target), FadeIn(obs))
        timing.hold_to_read(self, obs)

        eq = layout.explain_equation(
            self,
            [r"\varepsilon_\odot", r"\propto", r"\sin 2i_9"],
            [
                (0, "the Sun's induced obliquity"),
                (2, "strongest for P9 inclinations near 45°"),
            ],
            scale=0.95,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Analytic theory reproduces the 6° tilt for plausible P9 inclinations.")
