"""Batygin & Brown (2016) -- "Evidence for a Distant Giant Planet".

The founding paper: six stable distant KBOs share both their longitudes of
perihelion ϖ and their orbital poles, at a joint false-alarm probability of
~0.007% -- the clue that launched the search. Reproduced in p9-2016-evidence
(`joint_clustering_significance`, `observed_clustering_stats`).
"""
import numpy as np
from manim import (
    Arrow,
    Circle,
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    Indicate,
    MathTex,
    Scene,
    Text,
    UP,
    UR,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2016-evidence"


class Evidence2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Six orbits that line up")
        self.play(Write(tb[0]), FadeIn(tb[1], shift=UP * 0.15))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        # the six REAL stable KBOs (Batygin & Brown 2016): real e and apsidal
        # directions, with the measured clustering R-bar.
        swarm, r_bar = orbits.real_etno_swarm("kbo", color=P.GREEN)
        objs, _ = dataio.real_etnos("kbo")
        if objs:
            varpis = np.array([o["varpi_rad"] for o in objs])
        else:
            rng = np.random.default_rng(2016)
            varpis = np.deg2rad(255) + rng.normal(0, np.deg2rad(13), 6)
            swarm = orbits.etno_swarm([(2.7, 0.72, v) for v in varpis], color=P.GREEN)
        apses = VGroup(*[orbits.apse_arrow(2.7, 0.72, v, color=P.GREEN, scale=1.0) for v in varpis])
        self.play(FadeIn(sun))
        self.play(Create(swarm, lag_ratio=0.15), run_time=2.0)
        self.play(*[Create(a) for a in apses], run_time=1.0)
        self.play(Indicate(apses, color=P.TEAL, scale_factor=1.05))

        rbar = r_bar if r_bar is not None else orbits.mean_resultant_length(varpis)
        # the mean-resultant-length estimator for apsidal clustering
        eq = layout.explain_equation(
            self,
            [r"\bar R_\varpi", "=", r"\left|", r"\frac{1}{N}\sum_{k=1}^{N} e^{\,i\varpi_k}", r"\right|"],
            [
                (0, "mean resultant length (0=random, 1=aligned)"),
                (3, "average of the N unit apsidal vectors"),
            ],
            scale=0.85,
        )
        self.play(FadeOut(eq))
        stat = MathTex(rf"\bar{{R}}_\varpi \approx {rbar:.2f}", color=P.TEAL).scale(0.8)
        stat.to_corner(UR, buff=0.5)
        self.play(Write(stat))
        timing.hold_to_read(self, stat, settle=1.2)

        # the headline significance the crate reproduces
        readout = paper.result_readout("joint false-alarm probability", "P ≈ 0.007%", color=P.GREEN)
        readout.scale(0.9).to_edge(DOWN, buff=0.5)
        self.play(FadeIn(readout))
        pjoint = layout.equation(r"P_{\rm joint} \approx 7\times10^{-5}", color=P.GREEN).scale(0.7)
        pjoint.next_to(readout, UP, buff=0.3)
        self.play(Write(pjoint))
        timing.hold_to_read(self, readout, pjoint, settle=1.0)
        self.play(FadeOut(readout), FadeOut(pjoint))

        layout.show_takeaway(
            self, "Random orbits would not do this -- a massive shepherd is implied.")
