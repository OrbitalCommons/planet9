"""Linder & Mordasini (2016) -- thermal evolution & brightness of Planet Nine.

Cooling/evolution models give P9's luminosity and temperature as it ages, setting
its present-day brightness in reflected and thermal light. Reproduced in
p9-2016-linder-evolution.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, Scene, UP, UR, Write
import p9_manim as P
from p9_manim import layout, paper, widgets

CRATE = "p9-2016-linder-evolution"


class LinderEvolution2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Cooling for 4.5 billion years")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = widgets.axes([0, 4.5, 1], [0, 120, 30], x_length=9.0, y_length=3.8, font_size=16, shift_down=0.5)
        self.play(Create(ax),
                  FadeIn(layout.label("age (Gyr)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(layout.label("effective temperature (K)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        for m, col in [(10, P.BLUE), (5, P.TEAL)]:
            curve = ax.plot(lambda t, m=m: 30 + (m * 9) / (1 + t) ** 0.6, x_range=[0.05, 4.5, 0.05], color=col)
            self.play(Create(curve), run_time=1.0)
            self.add(layout.label(f"{m} M⊕", font_size=15, color=col).next_to(curve.get_start(), UP, buff=0.05))
        self.add(ax.get_vertical_line(ax.c2p(4.5, 40), color=P.MUTED, stroke_width=1.5))

        eq = layout.explain_equation(
            self,
            [r"T_{\rm eff}(t)\downarrow", r",\qquad", r"L", r"\propto", r"R^{2}\,T^{4}"],
            [
                (0, "the planet cools as it ages"),
                (2, "its luminosity (how bright it glows)"),
                (4, "falls with surface area times temperature⁴"),
            ],
            scale=0.8,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Today it sits near ~40 K -- warm enough to glow faintly in the far-IR.")
