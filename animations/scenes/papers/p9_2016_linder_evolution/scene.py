"""Linder & Mordasini (2016) -- thermal evolution & brightness of Planet Nine.

Cooling/evolution models give P9's luminosity and temperature as it ages, setting
its present-day brightness in reflected and thermal light. Reproduced in
p9-2016-linder-evolution.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2016-linder-evolution"


class LinderEvolution2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Cooling for 4.5 billion years")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 4.5, 1], y_range=[0, 120, 30], x_length=9.0, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax),
                  FadeIn(Text("age (Gyr)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("effective temperature (K)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        for m, col in [(10, P.BLUE), (5, P.TEAL)]:
            curve = ax.plot(lambda t, m=m: 30 + (m * 9) / (1 + t) ** 0.6, x_range=[0.05, 4.5, 0.05], color=col)
            self.play(Create(curve), run_time=1.0)
            self.add(Text(f"{m} M⊕", font_size=15, color=col).next_to(curve.get_start(), UP, buff=0.05))
        self.add(ax.get_vertical_line(ax.c2p(4.5, 40), color=P.MUTED, stroke_width=1.5))
        self.play(FadeIn(layout.takeaway("Today it sits near ~40 K -- warm enough to glow faintly in the far-IR.")))
        self.wait(0.9)
