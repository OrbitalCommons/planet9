"""Fortney et al. (2016) -- thermal evolution / emergent flux of Planet Nine.

A several-Earth-mass body retains internal heat, so it is warmer (and brighter in
the mid/far-IR) than pure solar equilibrium would predict -- which sets how
detectable it is. Reproduced in p9-2016-fortney-thermal.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, LEFT, Scene, UP, UR, Write,
)

import p9_manim as P
from p9_manim import layout, paper, widgets

CRATE = "p9-2016-fortney-thermal"


class FortneyThermal2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Still warm from birth")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = widgets.axes([2, 20, 4], [0, 60, 20], x_length=9.0, y_length=3.9, font_size=18, shift_down=0.5)
        xl = layout.label("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = layout.label("effective temperature (K)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        eq = ax.plot(lambda m: 10.0, x_range=[2, 20, 0.5], color=P.MUTED)
        warm = ax.plot(lambda m: 10.0 + 2.0 * m, x_range=[2, 20, 0.5], color=P.ORANGE)
        self.play(Create(eq))
        self.add(layout.label("solar equilibrium only", font_size=14, color=P.MUTED).next_to(ax.c2p(14, 10), UP, buff=0.05))
        self.play(Create(warm))
        self.add(layout.label("with internal heat", font_size=14, color=P.ORANGE).next_to(ax.c2p(14, 38), UP, buff=0.05))

        eq = layout.explain_equation(
            self,
            [r"T_{\rm eff}", r"=", r"\max", r"\left(T_{\rm eq},\ T_{\rm int}\right)"],
            [
                (0, "the temperature P9 actually radiates at"),
                (2, "whichever source dominates"),
                (3, "sunlight-only equilibrium vs leftover internal heat"),
            ],
            scale=0.8,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Leftover formation heat makes P9 warmer -- and easier for IR/mm surveys.")
