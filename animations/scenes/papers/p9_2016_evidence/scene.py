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
from p9_manim import layout, orbits, paper

CRATE = "p9-2016-evidence"


class Evidence2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Six orbits that line up")
        self.play(Write(tb[0]), FadeIn(tb[1], shift=UP * 0.15))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        # the six stable KBOs: long, eccentric, apsides clustered to one side
        rng = np.random.default_rng(2016)
        varpis = np.deg2rad(255) + rng.normal(0, np.deg2rad(13), 6)
        swarm = orbits.etno_swarm([(2.7, 0.72, v) for v in varpis], color=P.GREEN)
        apses = VGroup(*[orbits.apse_arrow(2.7, 0.72, v, color=P.GREEN, scale=1.0) for v in varpis])
        self.play(FadeIn(sun))
        self.play(Create(swarm, lag_ratio=0.15), run_time=2.0)
        self.play(*[Create(a) for a in apses], run_time=1.0)
        self.play(Indicate(apses, color=P.TEAL, scale_factor=1.05))

        rbar = orbits.mean_resultant_length(varpis)
        stat = MathTex(rf"\bar{{R}}_\varpi \approx {rbar:.2f}", color=P.TEAL).scale(0.8)
        stat.to_corner(UR, buff=0.6).shift(DOWN * 0.6)
        self.play(Write(stat))
        self.wait(0.3)

        # the headline significance the crate reproduces
        readout = paper.result_readout("joint false-alarm probability", "P ≈ 0.007%", color=P.GREEN)
        readout.scale(0.9).to_edge(DOWN, buff=0.5)
        self.play(FadeIn(readout))
        self.wait(0.6)
        self.play(FadeOut(readout))

        self.play(FadeIn(layout.takeaway(
            "Random orbits would not do this -- a massive shepherd is implied.")))
        self.wait(0.9)
