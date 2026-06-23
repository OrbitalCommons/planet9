"""Becker et al. (2018) -- discovery & analysis of 2015 BP519.

A real, extreme high-inclination (54 deg) ETNO -- exactly the kind of object P9
can produce, and hard to make otherwise. Reproduced in p9-2018-bp519
(pumping history, inclination).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, RIGHT, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2018-bp519"


class Bp519Discovery2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A real object that shouldn't exist")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        ecl = orbits.ellipse_orbit(2.8, 0.6, color=P.MUTED, varpi=0).set_stroke(opacity=0.4)
        # the high-inclination orbit (drawn tilted/edge-on via vertical squash)
        bp = orbits.ellipse_orbit(2.8, 0.6, color=P.GREEN, varpi=np.deg2rad(20))
        bp.apply_matrix(np.array([[1, 0, 0], [0, 0.35, 0], [0, 0, 1]]))
        self.play(FadeIn(sun), Create(ecl))
        self.add(Text("ecliptic", font_size=14, color=P.MUTED).next_to(ecl, DOWN, buff=0.05))
        self.play(Create(bp))
        clab = Text("2015 BP519: i ≈ 54°", font_size=18, color=P.GREEN).to_edge(UP, buff=1.5)
        self.play(FadeIn(clab))
        timing.hold_to_read(self, clab, settle=0.6)
        # the extreme inclination that marks it as a P9 product
        eq = layout.equation(r"i \approx 54^\circ", color=P.GREEN).scale(0.9).to_corner(UP + RIGHT, buff=0.5)
        layout.show_equation(self, eq, settle=1.0)
        ro = paper.result_readout("inclination (a P9 fingerprint)", "i ≈ 54°", color=P.GREEN)
        ro.scale(0.8).to_edge(DOWN, buff=1.3)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro, settle=1.2)
        layout.show_takeaway(
            self, "Steeply-inclined ETNOs are hard to make WITHOUT a sculptor like Planet Nine.")
