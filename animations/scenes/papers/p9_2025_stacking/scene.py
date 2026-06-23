"""Geringer-Sameth, Golovich & Iwabuchi (2025) -- multi-year stacking searches.

Matched-filter "digital tracking" over many years of images, with a metric on
orbit-parameter space and an honest look-elsewhere (trials) penalty. Reproduced
in p9-2025-stacking.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, Scene, Text, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2025-stacking"


class Stacking2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Searching a haystack of orbits")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # a grid of trial-orbit cells, one lighting up (the matched filter peak)
        cells = VGroup()
        rng = np.random.default_rng(7)
        for i in range(7):
            for j in range(4):
                v = rng.uniform(0.0, 0.25)
                c = Text("·", font_size=28, color=P.MUTED).move_to([-3.6 + i * 1.2, -0.4 + j * 0.9, 0]).set_opacity(0.4 + v)
                cells.add(c)
        self.play(FadeIn(cells, lag_ratio=0.01))
        peak = Text("●", font_size=40, color=P.GREEN).move_to([-3.6 + 4 * 1.2, -0.4 + 2 * 0.9, 0])
        self.play(FadeIn(peak))
        top = Text("matched-filter peak in (a, e, angles)", font_size=16, color=P.GREEN).to_edge(UP, buff=1.5)
        bot = Text("...but huge trials → look-elsewhere penalty", font_size=16, color=P.RED).to_edge(DOWN, buff=1.9)
        self.play(FadeIn(top), FadeIn(bot))
        timing.hold_to_read(self, top, bot, settle=0.5)

        eq = layout.equation(
            r"\mathrm{SNR} \propto \sqrt{N}\,,\quad \text{trials penalty}", color=P.GREEN, scale=0.8)
        eq.next_to(bot, UP, buff=0.3)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Stacking years of data is powerful -- if you count your trials honestly.")
