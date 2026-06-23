"""Clement & Sheppard (2020) -- precession / apsidal behaviour of distant TNOs.

Tracks the apsidal precession of the distant objects and how a Planet Nine sets
the precession pattern that maintains their alignment. Reproduced in
p9-2020-clement-precession.
"""
import numpy as np
from manim import Create, FadeIn, FadeOut, Scene, Text, UP, UR, Write, rate_functions, ValueTracker, always_redraw
import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2020-clement-precession"


class ClementPrecession2020(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Precessing in step")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        sun = orbits.sun()
        self.add(sun)
        # two orbits precessing together (coupled) -> alignment maintained
        d = ValueTracker(0.0)
        live = always_redraw(lambda: orbits.etno_swarm(
            [(2.5, 0.68, np.deg2rad(40) + d.get_value()),
             (2.7, 0.72, np.deg2rad(55) + d.get_value()),
             (2.3, 0.66, np.deg2rad(48) + d.get_value())], color=P.GREEN))
        self.add(live)
        self.add(Text("apsides drift together → stay clustered", font_size=17, color=P.TEAL).to_edge(UP, buff=1.5))
        self.play(d.animate.set_value(np.pi), run_time=4.0, rate_func=rate_functions.linear)

        eq = layout.explain_equation(
            self,
            [r"\dot\varpi_k", r"\approx", r"\text{const}", r"\Rightarrow", r"\Delta\varpi\ \text{preserved}"],
            [
                (0, "each object's apsides precess at nearly the same rate"),
                (4, "so the spread between them stays fixed — the cluster holds"),
            ],
            scale=0.78,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Coupled precession keeps the cluster intact as it slowly rotates.")
