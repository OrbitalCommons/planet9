"""Holman & Payne (2016) -- Cassini range observations of Saturn.

An independent analysis of the Earth-Saturn range: a too-massive or too-close P9
would have left a residual Cassini did not see, setting an upper bound on mass at
each distance. Reproduced in p9-2016-holman-payne (range residual model).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2016-holman-payne"


class HolmanPayne2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "What Saturn does NOT feel")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[200, 1000, 200], y_range=[0, 20, 5], x_length=9.3, y_length=4.0,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = Text("Planet Nine distance (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
        yl = Text("max allowed mass (M⊕)", font_size=16, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        # allowed mass grows steeply with distance (tidal residual ~ M / d^3)
        bound = ax.plot(lambda d: 0.6 * (d / 250.0) ** 3, x_range=[200, 1000, 5], color=P.ORANGE)
        self.play(Create(bound), run_time=1.5)
        excl = Text("excluded by Cassini ranging", font_size=18, color=P.RED).move_to(ax.c2p(420, 14))
        ok = Text("allowed", font_size=20, color=P.TEAL, weight="BOLD").move_to(ax.c2p(820, 4))
        self.play(FadeIn(excl), FadeIn(ok))
        timing.hold_to_read(self, excl, ok)

        eq = layout.equation_card(r"m_{\max}(d) \propto d^{3}", scale=1.0)
        eq.next_to(ax, UP, buff=0.15)
        layout.show_equation(self, eq)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self,
            "Ranging caps the mass at each distance -- closer planets must be lighter.")
