"""Holman et al. (2025) -- a PS1 shift-and-stack search.

Rather than linking catalogue detections, this co-adds Pan-STARRS pixels along
millions of trial orbits, reaching ~1 mag deeper than the catalogue search.
Reproduced in p9-2025-ps1-holman (stack depth, reach).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, Line, Scene, Text, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2025-ps1-holman"


class Ps1Holman2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Stacking pixels along a guess")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # several faint frames that, summed along the right track, reveal a source
        frames = VGroup()
        for k in range(5):
            f = Text("·", font_size=40, color=P.MUTED).move_to([-4 + k * 0.5, 1.2 - k * 0.2, 0]).set_opacity(0.4)
            frames.add(f)
        self.play(FadeIn(frames, lag_ratio=0.1))
        track = Line([-4, 1.2, 0], [-2, 0.4, 0], color=P.TEAL, stroke_width=2).set_stroke(opacity=0.6)
        self.play(Create(track))
        stacked = Text("●", font_size=44, color=P.GREEN).move_to([2.5, 0.6, 0])
        arrow = Line([-1.5, 0.6, 0], [2.0, 0.6, 0], color=P.MUTED, stroke_width=2)
        self.play(Create(arrow), FadeIn(stacked))
        note = Text("co-add along trial orbit → +1 mag depth (r ≈ 22.5)", font_size=17, color=P.TEAL).to_edge(DOWN, buff=1.9)
        self.play(FadeIn(note))
        timing.hold_to_read(self, note, settle=0.6)

        eq = layout.explain_equation(
            self,
            [r"\Delta m", r"\approx", "1", r"\ \text{mag}", r"\ (\text{shift--stack})"],
            [(0, "depth gained over a single image"),
             (2, "about one magnitude deeper"),
             (4, "from co-adding pixels along a trial orbit")],
            color=P.TEAL, scale=0.8, where=UP * 1.7)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Digital tracking digs below the single-image limit -- at the cost of many trials.")
