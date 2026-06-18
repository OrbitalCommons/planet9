"""Siraj et al. (2025) -- "Planet Y".

A separate, smaller (Mercury-to-Earth-mass) perturber, closer in, would warp the
mean plane of the Kuiper Belt -- a distinct signal from the distant Planet Nine.
Reproduced in p9-2025-planet-y (Laplace-plane warp).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, ParametricFunction, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2025-planet-y"


class PlanetY2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A different, closer culprit")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)
        # the Kuiper belt mean plane, warped by an inclined inner perturber
        flat = ParametricFunction(lambda t: np.array([t, 0, 0]), t_range=[-6, 6], color=P.MUTED, stroke_width=2)
        flat.set_stroke(opacity=0.4)
        warp = ParametricFunction(lambda t: np.array([t, 0.5 * np.sin(t * 0.6), 0]),
                                  t_range=[-6, 6], color=P.GREEN, stroke_width=3)
        flat_lbl = Text("unperturbed mean plane", font_size=15, color=P.MUTED).to_edge(DOWN, buff=1.6)
        self.play(Create(flat), FadeIn(flat_lbl))
        self.play(Create(warp), run_time=1.6)

        py = orbits.body_dot(1.4, 0.2, np.deg2rad(40), color=P.BLUE, radius=0.1)
        py_lbl = Text("'Planet Y' (smaller, closer, inclined)", font_size=15, color=P.BLUE).next_to(py, UP, buff=0.1)
        self.play(FadeIn(py), FadeIn(py_lbl))
        self.play(FadeIn(layout.takeaway(
            "Not every anomaly needs Planet Nine -- a nearer small world can warp the belt's plane.")))
        self.wait(0.9)
