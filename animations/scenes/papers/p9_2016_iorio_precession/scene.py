"""Iorio (2016) -- constraints from anomalous planetary perihelion precession.

A distant planet would add an anomalous precession to the inner planets' orbits;
the lack of a measured excess (esp. for Saturn) constrains P9's mass and
distance. Reproduced in p9-2016-iorio-precession.
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2016-iorio-precession"


class Iorio2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Reading the planets' wobble")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[200, 1000, 200], y_range=[0, 5, 1], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        xl = Text("Planet Nine distance (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("predicted excess precession", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        pred = ax.plot(lambda d: 60.0 / (d / 250.0) ** 3, x_range=[230, 1000, 5], color=P.ORANGE)
        bound = ax.plot(lambda d: 1.0, x_range=[200, 1000, 5], color=P.RED)
        self.play(Create(pred))
        self.play(Create(bound))
        self.add(Text("observed upper bound", font_size=15, color=P.RED).next_to(ax.c2p(800, 1), UP, buff=0.05))
        self.play(FadeIn(layout.takeaway(
            "No measured excess precession -> nearby/massive solutions are excluded.")))
        self.wait(0.9)
