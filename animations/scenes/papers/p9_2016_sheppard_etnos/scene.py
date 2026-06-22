"""Sheppard & Trujillo -- new extreme TNO discoveries.

New ETNOs found in dedicated deep surveys extend the clustered sample that
motivates Planet Nine. Reproduced in p9-2016-sheppard-etnos.
"""
import numpy as np
from manim import Create, FadeIn, Indicate, Scene, UP, Write
import p9_manim as P
from p9_manim import dataio, layout, orbits, paper

CRATE = "p9-2016-sheppard-etnos"


class SheppardEtnos2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "More objects, same direction")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        sun = orbits.sun()
        rng = np.random.default_rng(11)
        # the REAL ETNO sample as the existing clustered base
        base, _ = orbits.real_etno_swarm("etno", color=P.GREEN)
        if not base.submobjects:
            base = orbits.etno_swarm([(2.6, 0.72, np.deg2rad(250) + rng.normal(0, np.deg2rad(14)))
                                      for _ in range(5)], color=P.GREEN)
        self.play(FadeIn(sun), Create(base, lag_ratio=0.1), run_time=1.6)
        new = orbits.etno_swarm([(2.7, 0.74, np.deg2rad(250) + rng.normal(0, np.deg2rad(14)))
                                 for _ in range(3)], color=P.TEAL, stroke_width=2.5)
        self.play(Create(new, lag_ratio=0.15))
        self.play(Indicate(new, color=P.TEAL))
        self.play(FadeIn(layout.takeaway("Each newly-found ETNO has landed in the same apsidal cluster.")))
        self.wait(0.9)
