"""Kohne & Batygin (2020) -- secular dynamics to octupole order.

Going beyond the quadrupole, the octupole term in P9's secular potential
modulates the apsidal libration and widens the range of confined orbits.
Reproduced in p9-2020-secular-octupole (octupole modulation amplitude).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, Scene, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2020-secular-octupole"


class SecularOctupole2020(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The next term matters")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = widgets.axes([0, 12, 3], [-100, 100, 50], x_length=9.3, y_length=3.9, font_size=16, shift_down=0.4)
        xl = layout.label("time (arb.)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
        self.play(Create(ax), FadeIn(xl))

        quad = ax.plot(lambda t: 70 * np.sin(t), x_range=[0, 12, 0.05], color=P.MUTED)
        oct_ = ax.plot(lambda t: 70 * np.sin(t) * (1 + 0.3 * np.sin(0.3 * t)), x_range=[0, 12, 0.05], color=P.TEAL)
        self.play(Create(quad))
        self.add(layout.label("quadrupole only", font_size=15, color=P.MUTED).next_to(ax.c2p(10, 70), UP, buff=0.05))
        self.play(Create(oct_))
        olbl = layout.label("+ octupole modulation", font_size=15, color=P.TEAL).next_to(ax.c2p(3, 90), UP, buff=0.05)
        self.play(FadeIn(olbl))
        timing.hold_to_read(self, olbl, settle=0.4)

        eq = layout.explain_equation(
            self,
            [r"\epsilon_{\rm oct}", "=", r"\dfrac{a}{a_9}", r"\dfrac{e_9}{1-e_9^2}"],
            [
                (0, "strength of the octupole (next-order) term"),
                (2, "the TNO's size relative to P9's orbit"),
                (3, "grows with how eccentric Planet Nine is"),
            ],
            scale=0.9,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Higher-order terms reshape the libration -- and enlarge the confined zone.")
