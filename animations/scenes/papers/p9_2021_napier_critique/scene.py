"""Napier et al. (2021) -- No evidence for orbital clustering in the ETNOs.

Modelling the pointed surveys' selection functions, the apparent clustering is
consistent with detection biases: the data do not require a planet. Reproduced
in p9-2021-napier-critique.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, MathTex, Scene, Text, UP, UR, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper

CRATE = "p9-2021-napier-critique"


class NapierCritique2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The skeptics' rebuttal")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        # the SAME real ETNO objects -- but the critique argues they are
        # consistent with a uniform parent once each survey's bias is folded in
        swarm, _ = orbits.real_etno_swarm("etno", color=P.GREEN)
        if not swarm.submobjects:
            rng = np.random.default_rng(2021)
            varpis = np.deg2rad(250) + rng.normal(0, np.deg2rad(18), 7)
            swarm = orbits.etno_swarm([(2.6, 0.72, v) for v in varpis], color=P.GREEN)
        self.play(FadeIn(sun), Create(swarm, lag_ratio=0.1), run_time=1.8)

        p = MathTex(r"p\text{-value consistent with uniform}", color=P.RED).scale(0.6).to_corner(UR, buff=0.4)
        self.play(Write(p))
        self.play(FadeIn(layout.takeaway(
            "Account for where the surveys looked, and the 'clustering' may melt away.")))
        self.wait(0.9)
