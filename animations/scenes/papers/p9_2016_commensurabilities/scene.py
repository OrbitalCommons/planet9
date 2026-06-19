"""(2016) -- period commensurabilities of the distant TNOs with Planet Nine.

If P9 traps objects in resonance, their orbital periods should sit near simple
integer ratios of P9's period -- a measurable commensurability excess. Reproduced
in p9-2016-commensurabilities.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, Line, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2016-commensurabilities"


class Commensurabilities2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Periods that share a beat")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[0, 1, 0.25], y_range=[0, 1.05, 0.5], x_length=9.3, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax), FadeIn(Text("distance from nearest simple period ratio", font_size=17, color=P.FG).next_to(ax, DOWN, buff=0.25)))
        # an excess of objects near ratio (small distance) vs uniform
        excess = ax.plot(lambda x: 1.0 - 0.8 * x, x_range=[0, 1, 0.02], color=P.PURPLE)
        unif = ax.plot(lambda x: 0.5, x_range=[0, 1, 0.02], color=P.MUTED).set_stroke(opacity=0.5)
        self.play(Create(unif), Create(excess))
        self.add(Text("excess near commensurability", font_size=15, color=P.PURPLE).next_to(ax.c2p(0.2, 0.84), UP, buff=0.05))
        self.play(FadeIn(layout.takeaway("Periods clustering near integer ratios would betray resonant capture by P9.")))
        self.wait(0.9)
