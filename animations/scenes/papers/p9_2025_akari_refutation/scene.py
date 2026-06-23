"""(2025) -- refutation of the IRAS/AKARI Planet Nine candidate.

Follow-up shows the far-IR candidate's implied orbit is inconsistent with the
dynamical constraints (and/or the motion does not check out): a cautionary
counter-result. Reproduced in p9-2025-akari-refutation.
"""
import numpy as np
from manim import (
    Cross, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2025-akari-refutation"


class AkariRefutation2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A candidate that didn't survive")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        # the candidate's implied orbit, then crossed out as inconsistent
        cand = orbits.ellipse_orbit(3.0, 0.5, color=P.ORANGE, varpi=np.deg2rad(30))
        self.play(FadeIn(sun), Create(cand))
        self.add(Text("IRAS/AKARI candidate's implied orbit", font_size=16, color=P.ORANGE).to_edge(UP, buff=1.5))
        self.play(FadeIn(orbits.body_dot(3.0, 0.5, np.deg2rad(30), color=P.ORANGE, radius=0.1)))

        x = Cross(scale_factor=1.6, stroke_color=P.RED, stroke_width=8).move_to(cand.get_center())
        self.play(Create(x))
        note = Text("inconsistent with the orbit / motion constraints", font_size=16, color=P.RED).to_edge(DOWN, buff=1.9)
        self.play(FadeIn(note))
        timing.hold_to_read(self, note, settle=0.5)

        eq = layout.explain_equation(
            self,
            [r"\text{orbit}", r"\nRightarrow", r"\text{consistency}"],
            [(0, "the candidate's implied orbit"),
             (2, "fails the dynamical / motion constraints")],
            color=P.RED, scale=0.85, where=UP * 1.7)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Extraordinary claims need follow-up -- this candidate did not hold up.")
