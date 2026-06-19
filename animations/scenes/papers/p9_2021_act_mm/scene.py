"""(2021) -- ACT millimetre search for Planet Nine.

The Atacama Cosmology Telescope's mm maps can be searched for a moving
point-source; non-detection constrains warm/near solutions. Reproduced in
p9-2021-act-mm (mm detectability).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2021-act-mm"


class ActMm2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Hunting in the CMB maps")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[2, 20, 4], y_range=[0, 800, 200], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        xl = Text("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("ACT reach (AU)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        reach = ax.plot(lambda m: 430 + 24 * m, x_range=[2, 20, 0.5], color=P.PURPLE)
        self.play(Create(reach), run_time=1.3)
        self.add(Text("229 GHz, ~8 mJy point-source limit", font_size=16, color=P.PURPLE).next_to(ax, UP, buff=0.4))
        self.play(FadeIn(layout.takeaway(
            "CMB telescopes double as outer-solar-system surveys at millimetre wavelengths.")))
        self.wait(0.9)
