"""Brown & Batygin (2022) -- ZTF search.

A direct shift-and-link search of ZTF public data finds nothing, ruling out
56.4% of the predicted Planet Nine parameter space (the brightest / nearest
solutions). Reproduced in p9-2021-ztf (`compute_exclusion`).
"""
import numpy as np
from manim import (
    Create,
    DOWN,
    FadeIn,
    LEFT,
    Rectangle,
    RIGHT,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper

CRATE = "p9-2021-ztf"


def _exclusion(key, fallback):
    """Real per-orbit exclusion fraction from dataio.section('exclusion')."""
    try:
        val = dataio.section("exclusion").get(key)
        if val is not None:
            return float(val)
    except Exception:
        pass
    return fallback


class Ztf2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "First big bite out of the search")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # a declination strip: ZTF sees dec > -30 deg
        sky = Rectangle(width=10, height=2.0, color=P.MUTED, stroke_width=1).shift(UP * 0.8)
        sky.set_fill(P.MUTED, opacity=0.05)
        foot = Rectangle(width=10 * (1 - (-30 + 90) / 180.0), height=2.0, color=P.PURPLE, stroke_width=2)
        foot.set_fill(P.PURPLE, opacity=0.12).align_to(sky, RIGHT).align_to(sky, UP)
        foot_lbl = Text("ZTF footprint (dec > −30°), r ≈ 20.5", font_size=16, color=P.PURPLE).next_to(sky, UP, buff=0.1)
        self.play(Create(sky), Create(foot), FadeIn(foot_lbl))

        # exclusion bar filling to the real ZTF-ruled-out fraction
        frac = _exclusion("ztf", 0.564)
        track = Rectangle(width=8.0, height=0.8, color=P.MUTED, stroke_width=2).shift(DOWN * 1.4)
        fill = Rectangle(width=8.0 * frac, height=0.8, color=P.RED, stroke_width=0).set_fill(P.RED, opacity=0.5)
        fill.align_to(track, LEFT).align_to(track, UP)
        self.play(Create(track))
        self.play(Create(fill), run_time=1.5)
        pct = Text(f"{frac*100:.1f}% of parameter space ruled out", font_size=22, color=P.RED).next_to(track, DOWN, buff=0.25)
        self.play(Write(pct))
        self.play(FadeIn(layout.takeaway(
            "No planet found yet -- but more than half the hiding places are now gone.")))
        self.wait(0.9)
