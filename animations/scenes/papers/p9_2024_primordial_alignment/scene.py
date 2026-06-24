"""Huang & Gladman (2024) -- primordial orbital alignment of sednoids.

The sednoids' apsides may have been aligned at birth (e.g. by a stellar flyby or
the Sun's birth cluster) and simply preserved -- alignment without an ongoing
planet. Reproduced in p9-2024-primordial-alignment (lookback convergence).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, Scene, UP, VGroup, Write, rate_functions,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2024-primordial-alignment"


class PrimordialAlignment2024(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Aligned at birth?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)
        rng = np.random.default_rng(2024)
        # run TIME BACKWARD: today's spread converges to a tight primordial cluster
        today = np.deg2rad(250) + rng.normal(0, np.deg2rad(28), 7)
        birth = np.deg2rad(250) + rng.normal(0, np.deg2rad(6), 7)
        swarm = VGroup(*[orbits.ellipse_orbit(2.6, 0.72, color=P.GREEN, varpi=v, stroke_width=1.6, opacity=0.8)
                         for v in today])
        self.play(Create(swarm, lag_ratio=0.1), run_time=1.6)
        rwlab = layout.label("rewind to the Sun's birth cluster...", font_size=18, color=P.MUTED).to_edge(UP, buff=1.5)
        self.play(FadeIn(rwlab))
        timing.hold_to_read(self, rwlab, settle=0.5)
        # integrating apsidal precession backward to the primordial direction
        eq = layout.explain_equation(
            self,
            [r"\varpi(t)", "=", r"\varpi_0", "-", r"\dot\varpi", r"\,(t_0 - t)"],
            [
                (0, "apsidal direction at an earlier time t"),
                (2, "today's measured direction"),
                (4, "the precession rate that we rewind by"),
            ],
            scale=0.9,
        )
        self.play(FadeOut(eq))
        self.play(*orbits.precess(swarm, today, birth), run_time=3.0,
                  rate_func=rate_functions.smooth)
        layout.show_takeaway(
            self, "Maybe they were aligned at birth and froze -- no present-day planet required.")
