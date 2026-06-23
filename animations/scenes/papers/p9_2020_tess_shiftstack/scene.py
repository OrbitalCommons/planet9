"""Rice & Laughlin (2020) -- targeted shift-stacking of TESS for distant bodies.

Co-adding TESS frames along trial tracks (the "shift-stack") recovers known
distant TNOs and sets the method's reach for Planet Nine. Reproduced in
p9-2020-tess-shiftstack (depth gain, trial-track grid).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, Line, Scene, Text, UP, VGroup, Write,
    ValueTracker, always_redraw, rate_functions,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2020-tess-shiftstack"


class TessShiftstack2020(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Sweeping trial tracks across TESS")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # a sweeping trial-track angle; the correct angle stacks the signal
        ang = ValueTracker(-0.6)
        track = always_redraw(lambda: Line(
            [-4, -0.5, 0], [-4 + 6 * np.cos(ang.get_value()), -0.5 + 6 * np.sin(ang.get_value()), 0],
            color=P.TEAL, stroke_width=3))
        self.add(track)
        self.add(Text("trial track direction (∝ assumed orbit)", font_size=16, color=P.TEAL).to_edge(UP, buff=1.5))
        self.play(ang.animate.set_value(0.25), run_time=3.0, rate_func=rate_functions.smooth)
        found = Text("● recovered source", font_size=22, color=P.GREEN).move_to([1.5, 0.6, 0])
        self.play(FadeIn(found))
        timing.hold_to_read(self, found, settle=0.6)

        eq = layout.equation(r"\Delta m = 1.25\,\log_{10} N", color=P.GREEN, scale=0.85)
        eq.to_edge(DOWN, buff=1.5)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Pick the right track and the moving signal adds up; pick wrong and it smears.")
