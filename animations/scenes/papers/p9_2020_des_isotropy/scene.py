"""(2020) -- testing isotropy of detached objects with DES.

A two-angle correlation test on the DES distant-object sample: are the orbital
angles consistent with isotropy, or clustered? Reproduced in p9-2020-des-isotropy.
"""
import numpy as np
from manim import Arrow, Circle, Create, FadeIn, MathTex, Scene, UP, UR, VGroup, Write
import p9_manim as P
from p9_manim import dataio, layout, orbits, paper

CRATE = "p9-2020-des-isotropy"


class DesIsotropy2020(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Is the DES sample isotropic?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        sun = orbits.sun(radius=0.1)
        circ = Circle(radius=2.6, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        self.play(FadeIn(sun), Create(circ))
        # REAL ETNO apsidal directions, framed as the isotropy test
        objs, r_bar = dataio.real_etnos("etno")
        if objs:
            varpis = np.array([o["varpi_rad"] for o in objs])
        else:
            rng = np.random.default_rng(2020)
            varpis = rng.uniform(0, 2 * np.pi, 9)
        arr = VGroup(*[Arrow(np.zeros(3), 2.4 * np.array([np.cos(v), np.sin(v), 0]), color=P.PURPLE, buff=0, stroke_width=3)
                       for v in varpis])
        self.play(Create(arr, lag_ratio=0.07))
        rbar = r_bar if r_bar is not None else orbits.mean_resultant_length(varpis)
        self.play(Write(MathTex(rf"\bar{{R}} = {rbar:.2f}", color=P.RED).scale(0.6).to_corner(UR, buff=0.4)))
        self.play(FadeIn(layout.takeaway("DES's own footprint sees its objects fairly isotropic -- coverage matters.")))
        self.wait(0.9)
