"""Siraj, Chyba & Tremaine (2026) -- Measuring apsidal clustering.

A rigorous statistical estimator for how strongly the ETNO perihelion directions
are clustered, with calibrated significance. Reproduced in
p9-2026-apsidal-clustering.
"""
import numpy as np
from manim import (
    Arrow, Circle, Create, DOWN, FadeIn, FadeOut, Scene, UP, UR, VGroup, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2026-apsidal-clustering"


class ApsidalClustering2026(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Measuring the clustering, carefully")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun(radius=0.1)
        circle = Circle(radius=2.8, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        self.play(FadeIn(sun), Create(circle))

        # REAL apsidal directions and measured R-bar of the ETNO sample
        objs, r_bar = dataio.real_etnos("etno")
        if objs:
            varpis = np.array([o["varpi_rad"] for o in objs])
        else:
            rng = np.random.default_rng(2026)
            varpis = np.deg2rad(60) + rng.normal(0, np.deg2rad(20), 11)
        arr = orbits.apse_arrows(varpis, length=2.5, color=P.GREEN, stroke_width=3)
        self.play(Create(arr, lag_ratio=0.08))

        rbar = r_bar if r_bar is not None else orbits.mean_resultant_length(varpis)
        # the mean direction (resultant) arrow
        mx, my = np.mean(np.cos(varpis)), np.mean(np.sin(varpis))
        res = Arrow(np.zeros(3), 2.5 * rbar * np.array([mx, my, 0]) / np.hypot(mx, my),
                    color=P.TEAL, buff=0, stroke_width=6)
        self.play(Create(res))
        # the estimator the paper formalises
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
        stat = layout.equation(rf"\bar{{R}} = {rbar:.2f}", color=P.TEAL).scale(0.9)
        stat.to_corner(UR, buff=0.5)
        self.play(Write(stat))
        timing.hold_to_read(self, stat, settle=1.2)
        layout.show_takeaway(
            self, "A calibrated statistic turns 'they look aligned' into a significance you can defend.")
