"""Belyakov & Batygin (2025) -- perturbation theory for scattered-disk stability.

An analytic perturbation-theory treatment of how Planet Nine stabilises (or
destabilises) scattered-disk orbits, complementing the numerical stability maps.
Reproduced in p9-2025-perturbation.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, RIGHT, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2025-perturbation"


class Perturbation2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Stability, by hand")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[200, 800, 200], y_range=[0, 1.0, 0.5], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax),
                  FadeIn(Text("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("survival probability", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        curve = ax.plot(lambda a: 1.0 / (1.0 + np.exp(-(a - 420) / 50.0)), x_range=[200, 800, 5], color=P.TEAL)
        self.play(Create(curve), run_time=1.5)
        albl = Text("analytic boundary matches N-body", font_size=16, color=P.TEAL).next_to(ax, UP, buff=0.4)
        self.play(FadeIn(albl))
        timing.hold_to_read(self, albl, settle=0.4)

        eq = layout.equation_card(
            r"\langle H\rangle = \langle H_0\rangle + \epsilon\,\langle H_1\rangle"
        ).scale(0.8).to_corner(UP + RIGHT, buff=0.5).shift(DOWN * 0.7)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Perturbation theory reproduces the stability edge without giant simulations.")
