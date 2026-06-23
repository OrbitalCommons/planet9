"""Preface D -- How a hidden planet sculpts the belt.

Scene: P08Sculpting.
"""
import numpy as np
from manim import (
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    Indicate,
    Scene,
    Text,
    UP,
    UR,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, orbits, timing


class P08Sculpting(Scene):
    def construct(self):
        self.add(layout.concept_badge("PREFACE 08"))
        title = Text("How one planet herds a thousand orbits",
                     color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.6)
        self.play(Write(title))

        sun = orbits.sun()
        varpi_p9 = np.deg2rad(20)
        p9 = orbits.ellipse_orbit(3.7, 0.45, color=P.BLUE, varpi=varpi_p9)
        p9_apse = orbits.apse_arrow(3.7, 0.45, varpi_p9, color=P.BLUE)
        p9_lbl = Text("Planet Nine", color=P.BLUE, font_size=18).next_to(p9_apse.get_end(), UP, buff=0.1)
        self.play(FadeIn(sun), Create(p9), FadeIn(p9_apse), FadeIn(p9_lbl))

        rng = np.random.default_rng(3)
        n = 7
        start = rng.uniform(0, 2 * np.pi, n)
        # secular outcome: apsides confined ANTI-aligned with P9 (varpi_p9 + pi)
        target = varpi_p9 + np.pi + rng.normal(0, np.deg2rad(16), n)

        swarm = VGroup(*[orbits.ellipse_orbit(2.3, 0.7, color=P.GREEN, varpi=v,
                                              stroke_width=1.6, opacity=0.8) for v in start])
        self.play(Create(swarm, lag_ratio=0.1), run_time=1.8)
        cap = layout.caption("secular torque slowly rotates each orbit...")
        self.play(FadeIn(cap))
        self.wait(0.3)

        target_swarm = [orbits.ellipse_orbit(2.3, 0.7, color=P.GREEN, varpi=v,
                                             stroke_width=1.6, opacity=0.8) for v in target]
        self.play(*[s.animate.become(t) for s, t in zip(swarm, target_swarm)], run_time=3.0)
        self.play(FadeOut(cap))
        self.play(Indicate(swarm, color=P.TEAL, scale_factor=1.03))

        eq = layout.equation(r"\Delta\varpi \approx 180^\circ", color=P.TEAL).scale(0.8)
        eq.to_corner(UR, buff=0.6)
        self.play(Write(eq))
        timing.hold_to_read(self, eq, settle=0.4)

        # the secular Hamiltonian governs the apsidal confinement
        ham = layout.explain_equation(
            self,
            [r"\langle H\rangle", "=", r"H(e,\,\Delta\varpi)"],
            [
                (0, "orbit-averaged (secular) energy"),
                (2, "depends only on eccentricity e and apsidal offset Δϖ"),
            ],
            scale=0.9,
            where=DOWN * 1.4,
        )
        self.play(FadeOut(ham))

        layout.show_takeaway(
            self, "Confined apsides, lifted perihelia, tilted orbits -- one planet explains them all.")
