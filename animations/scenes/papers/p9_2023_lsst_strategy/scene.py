"""Schwamb et al. (2023) -- how LSST's observing strategy drives discovery.

The fraction of a Planet Nine population Rubin/LSST can discover depends on three
levers: per-visit depth, how many nights are required to link a track, and the
footprint's declination reach. Reproduced in p9-2023-lsst-strategy.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, LEFT, Rectangle, Scene, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2023-lsst-strategy"


class LsstStrategy2023(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Designing the survey that could win")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        levers = [
            ("depth per visit  (r ≈ 24.5)", 0.86, P.TEAL),
            ("nights to link a track", 0.74, P.PURPLE),
            ("footprint declination reach", 0.80, P.BLUE),
        ]
        rows = VGroup()
        for name, frac, col in levers:
            track = Rectangle(width=7.0, height=0.6, color=P.MUTED, stroke_width=1.5)
            fill = Rectangle(width=7.0 * frac, height=0.6, stroke_width=0, color=col).set_fill(col, opacity=0.5)
            fill.align_to(track, LEFT).align_to(track, UP)
            lab = layout.label(name, font_size=17, color=P.FG).next_to(track, UP, buff=0.08)
            rows.add(VGroup(track, fill, lab))
        rows.arrange(DOWN, buff=0.7).shift(DOWN * 0.6)
        for r in rows:
            self.play(Create(r[0]), Create(r[1]), FadeIn(r[2]), run_time=0.7)
        timing.hold_to_read(self, rows, settle=0.6)

        eq = layout.explain_equation(
            self,
            [r"f_{\rm disc}", "=", r"f(\text{depth},\,N_{\rm link},\,\delta_{\min})"],
            [(0, "fraction of P9 orbits LSST can discover"),
             (2, "set by depth, nights to link, declination reach")],
            color=P.TEAL, scale=0.78, where=UP * 2.5)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Tune depth, linking, and footprint -- and LSST discovers most of the viable region.")
