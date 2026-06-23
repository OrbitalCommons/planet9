"""Payne, Holman & Pal (2019) -- TESS yield via shift-and-stacking.

Co-adding ~thousands of TESS full-frame images along trial orbital tracks pushes
the depth to I_C ~ 22, deep enough to recover a fraction of Planet Nine
solutions. Reproduced in p9-2019-tess-yield (TessStack, detectable_fraction).
"""
import numpy as np
from manim import (
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    Line,
    RIGHT,
    Scene,
    Text,
    UP,
    Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2019-tess-yield"


def _optical_depth(name, fallback):
    """Real survey limiting magnitude from viability()'s survey_optical_depths."""
    try:
        for n, d in dataio.viability().get("survey_optical_depths", []):
            if n == name:
                return float(d)
    except Exception:
        pass
    return fallback


class TessYield2019(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Stacking thousands of frames")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = widgets.axes([18, 24, 1], [0, 1.05, 0.5], x_length=9.0, y_length=3.8, font_size=18, shift_down=0.5)
        xl = Text("apparent magnitude (I_C)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("detection efficiency", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, [-1, 0, 0], buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        depth = _optical_depth("TESS", 22.0)
        k = 4.0
        eff = ax.plot(lambda m: 1.0 / (1.0 + np.exp(k * (m - depth))), x_range=[18, 24, 0.05], color=P.PURPLE)
        self.play(Create(eff), run_time=1.6)

        dline = ax.get_vertical_line(ax.c2p(depth, 0.5), color=P.TEAL, stroke_width=2)
        dlbl = Text(f"stacked depth ≈ {depth:.1f}", font_size=18, color=P.TEAL).next_to(ax.c2p(depth, 1.0), UP, buff=0.1)
        self.play(Create(dline), FadeIn(dlbl))
        timing.hold_to_read(self, dlbl, settle=0.6)

        eq = layout.explain_equation(
            self,
            [r"m_{\rm stack}", "=", "m_1", "+", "1.25", r"\log_{10} N"],
            [(0, "depth reached after co-adding"),
             (2, "single-frame limiting magnitude"),
             (5, "gain grows with the log of frame count N")],
            color=P.TEAL, scale=0.8, where=UP * 2.45)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Even a survey built for exoplanets can be re-aimed at the outer dark.")
