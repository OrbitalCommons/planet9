"""Khain et al. (2018) -- the Kuiper belt's response to Planet Nine.

P9 broadens the perihelion distribution of the distant scattered disk -- pumping
some objects to large q ("broad" vs "narrow" outcomes). Reproduced in
p9-2018-kuiper-belt (population analysis).
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2018-kuiper-belt"


class KuiperBelt2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Spreading the perihelia")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[30, 100, 20], y_range=[0, 1.05, 0.5], x_length=9.0, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax), FadeIn(Text("perihelion q (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)))

        def g(x, mu, s):
            return np.exp(-0.5 * ((x - mu) / s) ** 2)
        narrow = ax.plot(lambda q: g(q, 38, 5), x_range=[30, 100, 0.5], color=P.MUTED)
        broad = ax.plot(lambda q: g(q, 60, 18), x_range=[30, 100, 0.5], color=P.GREEN)
        self.play(Create(narrow))
        self.add(Text("no P9 (narrow)", font_size=14, color=P.MUTED).next_to(ax.c2p(38, 1.0), UP, buff=0.05))
        self.play(Create(broad))
        self.add(Text("with P9 (broad, detached)", font_size=14, color=P.GREEN).next_to(ax.c2p(66, g(66, 60, 18)), UP, buff=0.1))
        self.play(FadeIn(layout.takeaway("P9 lifts a tail of objects to high perihelion -- detaching them from Neptune.")))
        self.wait(0.9)
