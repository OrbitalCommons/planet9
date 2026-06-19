"""(2019) clustering analysis -- quantifying the alignment with circular stats.

Uses Poincare-style variables and the mean resultant length to put a number and a
significance on the ETNO clustering. Reproduced in p9-2019-clustering.
"""
import numpy as np
from manim import Arrow, Circle, Create, FadeIn, MathTex, Scene, UP, UR, VGroup, Write
import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2019-clustering"


class Clustering2019(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Putting a number on it")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        sun = orbits.sun(radius=0.1)
        circ = Circle(radius=2.7, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        self.play(FadeIn(sun), Create(circ))
        rng = np.random.default_rng(2019)
        varpis = np.deg2rad(255) + rng.normal(0, np.deg2rad(16), 10)
        arr = VGroup(*[Arrow(np.zeros(3), 2.5 * np.array([np.cos(v), np.sin(v), 0]), color=P.GREEN, buff=0, stroke_width=3)
                       for v in varpis])
        self.play(Create(arr, lag_ratio=0.06))
        rbar = orbits.mean_resultant_length(varpis)
        self.play(Write(MathTex(rf"\bar{{R}} = {rbar:.2f}", color=P.TEAL).scale(0.9).to_corner(UR, buff=0.5)))
        self.play(FadeIn(layout.takeaway("Circular statistics convert 'looks aligned' into a defensible significance.")))
        self.wait(0.9)
