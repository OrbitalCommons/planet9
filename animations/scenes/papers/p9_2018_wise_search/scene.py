"""Meisner et al. / WISE W1 shift-and-stack search.

WISE scanned the whole sky in the thermal infrared, but at 3.4 µm (W1) a ~40 K
body is deep on the Wien tail -- so WISE constrains only relatively nearby /
massive Planet Nine solutions. Reproduced in p9-2018-wise-search.
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, LEFT, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2018-wise-search"


class WiseSearch2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "All-sky, but only so deep")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[2, 20, 4], y_range=[0, 400, 100], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.5)
        xl = Text("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("WISE W1 reach (AU)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        # shallow reach ~ grows slowly with mass (cold body faint at 3.4 um)
        reach = ax.plot(lambda m: 120 + 9.0 * m, x_range=[2, 20, 0.5], color=P.ORANGE)
        self.play(Create(reach), run_time=1.4)
        note = Text("a 40 K world barely emits at 3.4 µm", font_size=16, color=P.ORANGE).next_to(ax, UP, buff=0.1)
        self.play(FadeIn(note))
        self.play(FadeIn(layout.takeaway(
            "Whole-sky thermal coverage -- but shallow for the coldest, most distant cases.")))
        self.wait(0.9)
