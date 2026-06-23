"""Beust (2016) -- clustering by a secular resonance.

A distant KBO can be captured in a secular resonance with Planet Nine, where its
apsidal angle Δϖ librates about a fixed value instead of circulating -- locking
in the alignment. Reproduced in p9-2016-secular-resonance (libration center/amp).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write, always_redraw,
    rate_functions, ValueTracker, Dot,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2016-secular-resonance"


class SecularResonance2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Locked in by a slow resonance")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 10, 2], y_range=[-120, 120, 60], x_length=9.5, y_length=4.0,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = Text("time (arb.)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
        yl = Text("Δϖ − Δϖ₀  (deg)", font_size=16, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        lib = ax.plot(lambda t: 80 * np.sin(1.3 * t), x_range=[0, 10, 0.05], color=P.TEAL)
        self.play(Create(lib), run_time=2.0)
        t = ValueTracker(0.0)
        dot = always_redraw(lambda: Dot(ax.c2p(t.get_value(), 80 * np.sin(1.3 * t.get_value())),
                                        color=P.GREEN, radius=0.07))
        self.add(dot)
        self.play(t.animate.set_value(10.0), run_time=3.0, rate_func=rate_functions.linear)
        msg = Text("Δϖ oscillates -- it never runs away", font_size=18, color=P.TEAL).to_edge(UP, buff=1.5)
        self.play(FadeIn(msg))
        timing.hold_to_read(self, msg, settle=0.6)
        self.play(FadeOut(msg), FadeOut(yl))

        eq = layout.explain_equation(
            self,
            [r"\ddot{\Delta\varpi}", "=", "-", r"\,\omega_{\rm lib}^2", r"\,\Delta\varpi"],
            [
                (0, "how the apsidal angle accelerates"),
                (3, "the libration frequency, squared"),
                (4, "restoring force pulls it back -- it oscillates"),
            ],
            scale=0.95,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "A secular resonance holds the apsidal angle fixed -- alignment without luck.")
