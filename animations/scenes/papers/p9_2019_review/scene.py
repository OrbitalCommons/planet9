"""Batygin et al. (2019) -- the Planet Nine review.

The synthesis paper that revised the favoured parameters downward (about 5 Earth
masses, ~500 AU) and laid out the full evidentiary case. Reproduced in
p9-2019-review (revised parameters, detection prospects).
"""
import numpy as np
from manim import Create, DOWN, FadeIn, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import dataio, layout, orbits, paper

CRATE = "p9-2019-review"


class Review2019(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The case, consolidated")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        sun = orbits.sun()
        p9 = orbits.ellipse_orbit(3.2, 0.25, color=P.BLUE, varpi=np.deg2rad(20))
        swarm = orbits.etno_swarm([(2.3, 0.7, np.deg2rad(200 + 8 * k)) for k in range(5)], color=P.GREEN)
        self.play(FadeIn(sun), Create(p9), Create(swarm, lag_ratio=0.1))
        noms = dataio.section("nominal_orbits")
        nom = None
        if noms:
            nom = next((o for o in noms if "review" in o["name"] or "2019" in o["name"]), None)
        if nom:
            value = f"≈ {nom['mass_earth']:.0f} M⊕,  a ≈ {nom['a']:.0f} AU,  e ≈ {nom['e']:.2g}"
        else:
            value = "≈ 5 M⊕,  a ≈ 500 AU,  e ≈ 0.25"
        ro = paper.result_readout("revised best-fit", value, color=P.BLUE)
        ro.scale(0.8).to_edge(DOWN, buff=1.2)
        self.play(FadeIn(ro))
        self.play(FadeIn(layout.takeaway("By 2019 the favoured planet got smaller and closer -- and the case, broader.")))
        self.wait(0.9)
