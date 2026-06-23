"""Belyakov & Batygin (2025) -- perturbation theory for scattered-disk stability.

An analytic perturbation-theory treatment of how Planet Nine stabilises (or
destabilises) scattered-disk orbits, complementing the numerical stability maps.
Reproduced in p9-2025-perturbation.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2025-perturbation"


class Perturbation2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Stability, by hand")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = widgets.axes([200, 800, 200], [0, 1.0, 0.5], x_length=9.0, y_length=3.9, font_size=16, shift_down=0.5)
        self.play(Create(ax),
                  FadeIn(Text("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("survival probability", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        curve = ax.plot(lambda a: 1.0 / (1.0 + np.exp(-(a - 420) / 50.0)), x_range=[200, 800, 5], color=P.TEAL)
        self.play(Create(curve), run_time=1.5)
        albl = Text("analytic boundary matches N-body", font_size=16, color=P.TEAL).next_to(ax, UP, buff=0.4)
        self.play(FadeIn(albl))
        timing.hold_to_read(self, albl, settle=0.4)

        eq = layout.explain_equation(
            self,
            [r"\langle H\rangle", "=", r"\langle H_0\rangle", "+", r"\epsilon\,\langle H_1\rangle"],
            [
                (0, "the full orbit-averaged energy"),
                (2, "the simple, solvable leading part"),
                (4, "plus a small correction from P9 (the perturbation)"),
            ],
            scale=0.9,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Perturbation theory reproduces the stability edge without giant simulations.")
