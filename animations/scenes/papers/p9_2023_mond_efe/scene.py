"""Brown & Mathur (2023) -- MOND as an alternative to Planet Nine.

In Modified Newtonian Dynamics, the Milky Way's external gravitational field
imposes a weak, fixed-direction tide on very distant orbits, which can align
their apsides over Gyr -- producing clustering with NO ninth planet. Reproduced
in p9-2023-mond-efe.
"""
import numpy as np
from manim import (
    Arrow,
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    RIGHT,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2023-mond-efe"


class MondEfe2023(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Clustering without a planet?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)
        # the galactic external field direction (fixed)
        efe = Arrow(np.zeros(3), RIGHT * 3.0 + UP * 0.6, color=P.ORANGE, buff=0, stroke_width=4)
        efe_lbl = Text("Milky Way external field", font_size=16, color=P.ORANGE).next_to(efe.get_end(), UP, buff=0.1)
        self.play(Create(efe), FadeIn(efe_lbl))

        rng = np.random.default_rng(2023)
        start = rng.uniform(0, 2 * np.pi, 6)
        field_dir = np.arctan2(0.6, 3.0)
        target = field_dir + np.pi + rng.normal(0, np.deg2rad(15), 6)  # apsides align to the EFE tide
        swarm = VGroup(*[orbits.ellipse_orbit(2.4, 0.7, color=P.GREEN, varpi=v, stroke_width=1.6, opacity=0.8)
                         for v in start])
        self.play(Create(swarm, lag_ratio=0.1), run_time=1.6)
        self.play(*orbits.precess(swarm, start, target), run_time=3.0)

        eq = layout.explain_equation(
            self,
            [r"a_0", r"\approx", r"1.2\times10^{-10}\ \mathrm{m\,s^{-2}}", r";\quad",
             r"g_{\rm ext}\ \text{(galactic)}"],
            [
                (0, "the MOND acceleration scale, where gravity goes weird"),
                (2, "tiny — comparable to deep-space accelerations"),
                (4, "the Milky Way's external field tugs the whole Solar System"),
            ],
            scale=0.78,
            where=DOWN * 1.2)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "A serious alternative: the clustering may need modified gravity, not a planet.")
