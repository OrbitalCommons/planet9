"""Caceres & Gomes (2018) -- TNO evolution under a low-perihelion Planet Nine.

If P9 itself has a relatively low perihelion, it scatters and reshapes the distant
population differently -- a variant that changes which orbits get confined.
Reproduced in p9-2018-low-perihelion.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, RIGHT, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2018-low-perihelion"


class LowPerihelion2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A planet that dips closer in")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        # high-q P9 vs low-q P9 (more eccentric, perihelion nearer the planets)
        high = orbits.ellipse_orbit(3.4, 0.3, color=P.MUTED, varpi=np.deg2rad(20)).set_stroke(opacity=0.5)
        low = orbits.ellipse_orbit(3.4, 0.7, color=P.BLUE, varpi=np.deg2rad(20))
        self.play(FadeIn(sun), Create(high))
        self.add(Text("high-perihelion P9", font_size=14, color=P.MUTED).next_to(high, UP, buff=0.05))
        self.play(Create(low))
        llbl = Text("low-perihelion P9 (more eccentric)", font_size=15, color=P.BLUE).to_edge(DOWN, buff=1.5)
        self.play(FadeIn(llbl))
        timing.hold_to_read(self, llbl, settle=0.4)

        eq = layout.equation_card(r"q_9 = a_9\,(1-e_9)").scale(0.8).to_corner(UP + RIGHT, buff=0.5).shift(DOWN * 0.6)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "P9's own perihelion is a free parameter -- and it changes how it sculpts the belt.")
