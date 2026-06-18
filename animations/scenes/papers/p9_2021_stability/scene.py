"""Batygin, Mardling & Nesvorny (2021) -- the stability boundary of the disk.

Distant scattered-disk orbits survive the age of the solar system only outside a
perihelion/semimajor-axis boundary; inside it they are cleared. This boundary
shapes which P9 orbits are dynamically allowed. Reproduced in p9-2021-stability.
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2021-stability"


class Stability2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The edge of survival")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[200, 1000, 200], y_range=[30, 120, 30], x_length=9.5, y_length=4.0,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = Text("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
        yl = Text("perihelion q (AU)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        boundary = ax.plot(lambda a: 40 + 0.045 * (a - 200), x_range=[200, 1000, 5], color=P.RED)
        self.play(Create(boundary), run_time=1.4)
        unstable = Text("cleared (unstable)", font_size=18, color=P.RED).move_to(ax.c2p(420, 45))
        stable = Text("survives", font_size=20, color=P.TEAL, weight="BOLD").move_to(ax.c2p(700, 100))
        self.play(FadeIn(unstable), FadeIn(stable))
        self.play(FadeIn(layout.takeaway(
            "Only detached-enough orbits last 4.5 Gyr -- which bounds the allowed P9 family.")))
        self.wait(0.9)
