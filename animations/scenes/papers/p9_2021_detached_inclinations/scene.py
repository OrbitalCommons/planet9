"""Anderson & Kaib (2021) -- Planet Nine and the inclinations of detached KBOs.

P9 forces a non-zero inclination on the detached population: their orbital planes
are pumped and broadened in a way a P9-free model does not reproduce. Reproduced
in p9-2021-detached-inclinations (forced inclination).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2021-detached-inclinations"


class DetachedInclinations2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Tilting the detached belt")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 40, 10], y_range=[0, 1.05, 0.5], x_length=9.3, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = Text("inclination (deg)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
        yl = Text("relative number", font_size=16, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        def gauss(x, mu, s):
            return np.exp(-0.5 * ((x - mu) / s) ** 2)

        no_p9 = ax.plot(lambda i: gauss(i, 6, 4), x_range=[0, 40, 0.5], color=P.MUTED)
        with_p9 = ax.plot(lambda i: gauss(i, 18, 9), x_range=[0, 40, 0.5], color=P.GREEN)
        self.play(Create(no_p9))
        self.add(Text("no Planet Nine", font_size=15, color=P.MUTED).next_to(ax.c2p(6, 1.0), UP, buff=0.05))
        self.play(Create(with_p9))
        wlbl = Text("with Planet Nine", font_size=15, color=P.GREEN).next_to(ax.c2p(22, gauss(22, 18, 9)), UP, buff=0.1)
        self.play(FadeIn(wlbl))
        timing.hold_to_read(self, wlbl, settle=0.4)

        eq = layout.equation_card(r"i_{\rm forced}(a)").scale(0.95).to_edge(DOWN, buff=1.6)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "P9 forces broader, higher inclinations in the detached belt than no-planet models.")
