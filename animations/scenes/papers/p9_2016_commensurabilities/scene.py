"""(2016) -- period commensurabilities of the distant TNOs with Planet Nine.

If P9 traps objects in resonance, their orbital periods should sit near simple
integer ratios of P9's period -- a measurable commensurability excess. Reproduced
in p9-2016-commensurabilities.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, Line, RIGHT, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2016-commensurabilities"


class Commensurabilities2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Periods that share a beat")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = widgets.axes([0, 1, 0.25], [0, 1.05, 0.5], x_length=9.3, y_length=3.8, font_size=16, shift_down=0.5)
        self.play(Create(ax), FadeIn(Text("distance from nearest simple period ratio", font_size=17, color=P.FG).next_to(ax, DOWN, buff=0.25)))
        # an excess of objects near ratio (small distance) vs uniform
        excess = ax.plot(lambda x: 1.0 - 0.8 * x, x_range=[0, 1, 0.02], color=P.PURPLE)
        unif = ax.plot(lambda x: 0.5, x_range=[0, 1, 0.02], color=P.MUTED).set_stroke(opacity=0.5)
        self.play(Create(unif), Create(excess))
        elbl = Text("excess near commensurability", font_size=15, color=P.PURPLE).next_to(ax.c2p(0.2, 0.84), UP, buff=0.05)
        self.play(FadeIn(elbl))
        timing.hold_to_read(self, elbl, settle=0.4)

        eq = layout.equation_card(
            r"\dfrac{P_{\rm TNO}}{P_9} \approx \dfrac{p}{q}"
        ).scale(0.7).to_edge(RIGHT, buff=0.5).shift(UP * 1.2)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Periods clustering near integer ratios would betray resonant capture by P9.")
