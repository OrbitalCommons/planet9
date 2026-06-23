"""(2019) clustering analysis -- quantifying the alignment with circular stats.

Uses Poincare-style variables and the mean resultant length to put a number and a
significance on the ETNO clustering. Reproduced in p9-2019-clustering.
"""
import numpy as np
from manim import Arrow, Circle, Create, DOWN, FadeIn, FadeOut, Scene, UP, UR, VGroup, Write
import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2019-clustering"


class Clustering2019(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Putting a number on it")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        sun = orbits.sun(radius=0.1)
        circ = Circle(radius=2.7, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        self.play(FadeIn(sun), Create(circ))
        # REAL apsidal directions of the Brown 2017 ETNO sample
        objs, r_bar = dataio.real_etnos("etno")
        if objs:
            varpis = np.array([o["varpi_rad"] for o in objs])
        else:
            rng = np.random.default_rng(2019)
            varpis = np.deg2rad(255) + rng.normal(0, np.deg2rad(16), 10)
        arr = VGroup(*[Arrow(np.zeros(3), 2.5 * np.array([np.cos(v), np.sin(v), 0]), color=P.GREEN, buff=0, stroke_width=3)
                       for v in varpis])
        self.play(Create(arr, lag_ratio=0.06))
        rbar = r_bar if r_bar is not None else orbits.mean_resultant_length(varpis)
        # the main circular-statistics estimator, walked term by term
        eq = layout.explain_equation(
            self,
            [r"\bar R", "=", r"\left|", r"\tfrac1N\sum e^{\,i\varpi_k}", r"\right|"],
            [
                (0, "mean resultant length (0=random, 1=aligned)"),
                (3, "average of the unit apsidal vectors"),
            ],
            scale=0.9,
        )
        self.play(FadeOut(eq))
        # its companion significance test, and the measured value
        pray = layout.equation(r"p_{\rm Rayleigh} = e^{-N\bar R^2}", color=P.GREEN).scale(0.6)
        stat = layout.equation(rf"\bar{{R}} = {rbar:.2f}", color=P.TEAL).scale(0.9)
        VGroup(pray, stat).arrange(DOWN, buff=0.35).to_corner(UR, buff=0.5)
        self.play(Write(pray))
        timing.hold_to_read(self, pray, settle=0.8)
        self.play(Write(stat))
        timing.hold_to_read(self, stat, settle=1.2)
        layout.show_takeaway(self, "Circular statistics convert 'looks aligned' into a defensible significance.")
