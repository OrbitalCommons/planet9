"""Resonance hopping -- objects jump between resonances.

A distant TNO is not locked forever: it sticks in one resonance, escapes, and is
captured by a neighbour, random-walking in semimajor axis. Reproduced in
p9-2017-resonance-hopping.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write, rate_functions, ValueTracker, always_redraw, Dot
import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2017-resonance-hopping"


class ResonanceHopping2017(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Sticking, then hopping")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[0, 10, 2], y_range=[300, 500, 50], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.4)
        self.play(Create(ax),
                  FadeIn(Text("time", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)),
                  FadeIn(Text("semi-major axis (AU)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        for lvl in (340, 390, 440):
            self.add(ax.plot(lambda x, l=lvl: l, x_range=[0, 10, 5], color=P.MUTED).set_stroke(opacity=0.4))
        # stair-step random walk between resonance levels
        steps = [340, 340, 390, 390, 390, 440, 390, 440, 440]
        def val(t):
            i = min(int(t), len(steps) - 1)
            return steps[i]
        tr = ValueTracker(0.0)
        dot = always_redraw(lambda: Dot(ax.c2p(tr.get_value(), val(tr.get_value())), color=P.GREEN, radius=0.08))
        self.add(dot)
        self.play(tr.animate.set_value(8.9), run_time=5.0, rate_func=rate_functions.linear)

        eq = layout.equation_card(
            r"a:\ \text{random walk among } a_{\rm res}"
        ).scale(0.85).to_edge(DOWN, buff=1.6)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Resonances trap, release, and re-trap -- so a's spread over time.")
