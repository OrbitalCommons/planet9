"""WISE co-add search.

Stacking N WISE exposures deepens the limit by ~2.5·log10(√N), and the W2 band
(4.6 µm) catches a cold body better than W1 -- extending the thermal reach.
Reproduced in p9-2016-wise-coadd (coadd depth, band advantage).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, LEFT, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2016-wise-coadd"


class WiseCoadd2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Deeper by stacking")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[1, 1000, 200], y_range=[15, 19, 1], x_length=9.0, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        xl = Text("frames stacked (N)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("limiting magnitude", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        depth = ax.plot(lambda n: 16.5 + 1.25 * np.log10(max(n, 1)), x_range=[1, 1000, 5], color=P.ORANGE)
        self.play(Create(depth), run_time=1.6)

        eq = layout.equation_card(
            r"m_{\rm lim} \to m_0 + 1.25\,\log_{10}\!\sqrt{N}", color=P.ORANGE, scale=0.85)
        eq.next_to(ax.c2p(550, 16.5 + 1.25 * np.log10(550)), UP, buff=0.15)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Stacking + the W2 band stretch WISE's cold-body reach -- but only so far.")
