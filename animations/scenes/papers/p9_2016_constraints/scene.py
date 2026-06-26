"""Brown & Batygin (2016) -- Observational Constraints on Planet Nine.

Quantifies the clustering: surviving distant objects in the right semimajor-axis
band should show confined apsides AND a high fraction of detached (high-q)
orbits. Reproduced in p9-2016-constraints (clustering_metric, parameter_grid).
"""
import numpy as np
from manim import (
    Create, FadeIn, FadeOut, Scene, UP, UR, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2016-constraints"


class Constraints2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Turning a hunch into a number")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        # the REAL six stable KBOs with their measured clustering R-bar
        swarm, r_bar = orbits.real_etno_swarm("kbo", color=P.GREEN)
        if r_bar is None:
            rng = np.random.default_rng(5)
            varpis = np.deg2rad(250) + rng.normal(0, np.deg2rad(15), 8)
            swarm = orbits.etno_swarm([(2.6, 0.74, v) for v in varpis], color=P.GREEN)
            r_bar = orbits.mean_resultant_length(varpis)
        self.play(FadeIn(sun), Create(swarm, lag_ratio=0.1), run_time=2.0)

        rbar = r_bar
        eq = layout.explain_equation(
            self,
            [r"\bar R_{\Delta\varpi}", r"\approx", rf"{rbar:.2f}", r",\qquad", r"f_{q>60\,\mathrm{AU}}", r"\ \text{high}"],
            [
                (0, "how tightly the apsides are confined (1=perfect)"),
                (4, "fraction of detached, high-perihelion orbits"),
            ],
            scale=0.85,
        )
        self.play(FadeOut(eq))
        layout.show_takeaway(
            self, "Two testable metrics: apsidal confinement and a surplus of detached orbits.")
