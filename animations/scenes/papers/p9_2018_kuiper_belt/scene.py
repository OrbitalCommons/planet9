"""Khain et al. (2018) -- the Kuiper belt's response to Planet Nine.

P9 broadens the perihelion distribution of the distant scattered disk -- pumping
some objects to large q ("broad" vs "narrow" outcomes). Reproduced in
p9-2018-kuiper-belt (population analysis).
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, Scene, UP, UR, Write
import p9_manim as P
from p9_manim import layout, paper, widgets

CRATE = "p9-2018-kuiper-belt"


class KuiperBelt2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Spreading the perihelia")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = widgets.axes([30, 100, 20], [0, 1.05, 0.5], x_length=9.0, y_length=3.8, font_size=16, shift_down=0.5)
        self.play(Create(ax), FadeIn(layout.label("perihelion q (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)))

        def g(x, mu, s):
            return np.exp(-0.5 * ((x - mu) / s) ** 2)
        narrow = ax.plot(lambda q: g(q, 38, 5), x_range=[30, 100, 0.5], color=P.MUTED)
        broad = ax.plot(lambda q: g(q, 60, 18), x_range=[30, 100, 0.5], color=P.GREEN)
        self.play(Create(narrow))
        self.add(layout.label("no P9 (narrow)", font_size=14, color=P.MUTED).next_to(ax.c2p(38, 1.0), UP, buff=0.05))
        self.play(Create(broad))
        self.add(layout.label("with P9 (broad, detached)", font_size=14, color=P.GREEN).next_to(ax.c2p(66, g(66, 60, 18)), UP, buff=0.1))

        eq = layout.explain_equation(
            self,
            [r"e", r"=", r"1", r"-", r"\dfrac{q}{a}"],
            [
                (0, "eccentricity — how stretched the orbit is"),
                (4, "perihelion over semi-major axis"),
            ],
            scale=0.9,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "P9 lifts a tail of objects to high perihelion -- detaching them from Neptune.")
