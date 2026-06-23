"""Hadden, Li, Payne & Holman (2018) -- Chaotic dynamics from Planet Nine.

Each mean-motion resonance with P9 has a finite width in semimajor axis; when
neighbouring resonances overlap (Chirikov's criterion), motion becomes chaotic
and objects diffuse. Reproduced in p9-2018-chaotic-dynamics (libration widths).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, Line, Rectangle, Scene, Text, UP, VGroup, Write,
    rate_functions, ValueTracker, always_redraw,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2018-chaotic-dynamics"


class ChaoticDynamics2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "When resonances collide")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        axis = Line([-6, -1.5, 0], [6, -1.5, 0], color=P.FG, stroke_width=2)
        axlbl = Text("semi-major axis", font_size=18, color=P.MUTED).next_to(axis, DOWN, buff=0.2)
        self.play(Create(axis), FadeIn(axlbl))

        centers = np.linspace(-4.5, 4.5, 7)
        labels = ["3:2", "2:1", "5:2", "3:1", "7:2", "4:1", "9:2"]
        w = ValueTracker(0.25)
        bands = VGroup(*[
            always_redraw(lambda c=c: Rectangle(width=w.get_value(), height=2.4, stroke_width=0,
                                                color=P.PURPLE).set_fill(P.PURPLE, opacity=0.35)
                          .move_to([c, -0.3, 0]))
            for c in centers])
        ticks = VGroup(*[Text(t, font_size=14, color=P.PURPLE).move_to([c, 1.1, 0]) for c, t in zip(centers, labels)])
        self.add(bands); self.play(FadeIn(ticks))
        msg = Text("each resonance has a width", font_size=18, color=P.PURPLE).to_edge(UP, buff=1.6)
        self.play(FadeIn(msg))
        timing.hold_to_read(self, msg, settle=0.4)

        # widen until they overlap -> chaos
        self.play(w.animate.set_value(1.6), run_time=3.0, rate_func=rate_functions.smooth)

        eq = layout.explain_equation(
            self,
            [r"K", "=", r"\dfrac{\delta a_{\rm res}}{\Delta a_{\rm sep}}", r"\gtrsim", "1"],
            [
                (0, "Chirikov resonance-overlap parameter"),
                (2, "resonance width divided by spacing to its neighbour"),
                (4, "once they touch (K>1), motion turns chaotic"),
            ],
            scale=0.9,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Overlapping resonances -> chaos: distant orbits diffuse over billions of years.")
