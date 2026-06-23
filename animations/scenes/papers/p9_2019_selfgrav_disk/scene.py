"""Sefilian & Touma (2019) -- shepherding by a self-gravitating disk.

A massive disk of trans-Neptunian debris can collectively shepherd ETNO apsides
into alignment -- producing the clustering with NO ninth planet. Reproduced in
p9-2019-selfgrav-disk.
"""
import numpy as np
from manim import (
    Annulus, Create, DOWN, FadeIn, RIGHT, Scene, Text, UP, VGroup, Write, rate_functions,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2019-selfgrav-disk"


class SelfgravDisk2019(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Maybe it's the disk itself")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        disk = Annulus(inner_radius=1.6, outer_radius=3.4, color=P.PURPLE, fill_opacity=0.12, stroke_width=1)
        disk_lbl = Text("massive self-gravitating debris disk", font_size=16, color=P.PURPLE).to_edge(UP, buff=1.5)
        self.play(FadeIn(sun), FadeIn(disk), FadeIn(disk_lbl))

        rng = np.random.default_rng(2019)
        start = rng.uniform(0, 2 * np.pi, 6)
        target = np.deg2rad(70) + rng.normal(0, np.deg2rad(15), 6)
        swarm = VGroup(*[orbits.ellipse_orbit(2.5, 0.65, color=P.GREEN, varpi=v, stroke_width=1.6, opacity=0.85)
                         for v in start])
        self.play(Create(swarm, lag_ratio=0.1), run_time=1.6)
        targ = [orbits.ellipse_orbit(2.5, 0.65, color=P.GREEN, varpi=v, stroke_width=1.6, opacity=0.85) for v in target]
        self.play(*[s.animate.become(t) for s, t in zip(swarm, targ)], run_time=3.0, rate_func=rate_functions.smooth)

        eq = layout.equation_card(
            r"\dot\varpi_{\rm disk} \propto \dfrac{M_{\rm disk}}{M_\odot}\,n"
        ).scale(0.8).to_corner(UP + RIGHT, buff=0.5).shift(DOWN * 0.7)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Collective disk gravity can mimic a planet -- if enough mass survives out there.")
