"""(2021) -- the perihelion gap.

The distant TNOs show a deficit of objects at intermediate perihelion (the "gap")
separating the scattering disk from the detached population -- a feature P9 can
carve. Reproduced in p9-2021-perihelion-gap.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, UR, Write
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2021-perihelion-gap"


class PerihelionGap2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The gap in between")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[33, 80, 10], y_range=[0, 1.05, 0.5], x_length=9.3, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax), FadeIn(Text("perihelion q (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)))

        def dist(q):
            return np.exp(-0.5 * ((q - 37) / 2.5) ** 2) + 0.8 * np.exp(-0.5 * ((q - 65) / 8) ** 2)
        curve = ax.plot(dist, x_range=[33, 80, 0.4], color=P.GREEN)
        self.play(Create(curve), run_time=1.5)
        gapline = ax.get_vertical_line(ax.c2p(48, dist(48)), color=P.RED, stroke_width=2)
        self.add(gapline, Text("the gap (~42–55 AU)", font_size=15, color=P.RED).next_to(ax.c2p(48, 0.5), UP, buff=0.1))

        eq = layout.equation_card(
            r"N(q):\ \text{deficit at}\ 42 \lesssim q \lesssim 55\ \mathrm{AU}",
            scale=0.74).to_corner(UR, buff=0.5)
        layout.show_equation(self, eq)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "A scarcity at middling perihelia separates scattered from detached -- a P9 carve-out.")
