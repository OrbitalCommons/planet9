"""(2025) -- new distant-TNO discoveries that stress the clustering picture.

Objects like 2017 OF201 and Ammonite (2023 KQ14) extend the sample; some fit the
P9 picture and some sit awkwardly, tightening the test. Reproduced in
p9-2025-new-discoveries.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, Indicate, Scene, Text, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2025-new-discoveries"


class NewDiscoveries2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "New objects, sharper test")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        rng = np.random.default_rng(2025)
        known = np.deg2rad(250) + rng.normal(0, np.deg2rad(15), 7)
        swarm = orbits.etno_swarm([(2.6, 0.72, v) for v in known], color=P.GREEN, opacity=0.7)
        self.play(FadeIn(sun), Create(swarm, lag_ratio=0.08), run_time=1.6)
        self.add(Text("known clustered ETNOs", font_size=15, color=P.GREEN).to_edge(UP, buff=1.5))

        # two new objects: one aligned, one off
        new_aligned = orbits.ellipse_orbit(2.7, 0.7, color=P.TEAL, varpi=np.deg2rad(248), stroke_width=2.5)
        new_off = orbits.ellipse_orbit(2.7, 0.7, color=P.ORANGE, varpi=np.deg2rad(120), stroke_width=2.5)
        self.play(Create(new_aligned))
        self.add(Text("2017 OF201 (fits)", font_size=14, color=P.TEAL).to_edge(DOWN, buff=1.7))
        self.play(Create(new_off))
        self.add(Text("Ammonite (off the cluster)", font_size=14, color=P.ORANGE).to_edge(DOWN, buff=1.3))
        self.play(Indicate(new_off, color=P.ORANGE))
        self.play(FadeIn(layout.takeaway(
            "Each new distant world is a fresh, unforgiving test of the hypothesis.")))
        self.wait(0.9)
