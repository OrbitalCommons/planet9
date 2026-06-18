"""Brown & Batygin (2016) -- Observational Constraints on Planet Nine.

Quantifies the clustering: surviving distant objects in the right semimajor-axis
band should show confined apsides AND a high fraction of detached (high-q)
orbits. Reproduced in p9-2016-constraints (clustering_metric, parameter_grid).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, MathTex, Scene, Text, UP, UR, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2016-constraints"


class Constraints2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Turning a hunch into a number")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        rng = np.random.default_rng(5)
        varpis = np.deg2rad(250) + rng.normal(0, np.deg2rad(15), 8)
        swarm = orbits.etno_swarm([(2.6, 0.74, v) for v in varpis], color=P.GREEN)
        self.play(FadeIn(sun), Create(swarm, lag_ratio=0.1), run_time=2.0)

        rbar = orbits.mean_resultant_length(varpis)
        m1 = MathTex(rf"\bar{{R}}_{{\Delta\varpi}} \approx {rbar:.2f}", color=P.TEAL).scale(0.7)
        m2 = MathTex(r"f_{q>60\,\mathrm{AU}}\ \text{high}", color=P.GREEN).scale(0.7)
        VGroup(m1, m2).arrange(DOWN, buff=0.35).to_corner(UR, buff=0.5)
        self.play(Write(m1)); self.play(Write(m2))
        self.play(FadeIn(layout.takeaway(
            "Two testable metrics: apsidal confinement and a surplus of detached orbits.")))
        self.wait(0.9)
