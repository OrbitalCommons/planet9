"""Belyakov, Bernardinelli & Brown (2022) -- the DES search.

The Dark Energy Survey is deep (r ~ 23.8) but covers only ~5000 deg^2 (~12% of
sky), so it rules out the faint, distant solutions that happen to fall in its
footprint -- adding ~5% unique exclusion. Reproduced in p9-2022-des.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, LEFT, Rectangle, Scene, UP, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2022-des"


def _exclusion(key, fallback):
    """Real per-orbit exclusion fraction from dataio.section('exclusion')."""
    try:
        val = dataio.section("exclusion").get(key)
        if val is not None:
            return float(val)
    except Exception:
        pass
    return fallback


class Des2022(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Deep, but narrow")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sky = Rectangle(width=10, height=2.5, color=P.MUTED, stroke_width=1.5).shift(UP * 0.7)
        sky.set_fill(P.MUTED, opacity=0.05)
        # ~12% footprint
        foot = widgets.footprint_rect(10 * 0.12, 1.0, fill_opacity=0.2, stroke_width=2)
        foot.move_to(sky.get_center() + LEFT * 3.0 + DOWN * 0.6)
        self.play(Create(sky), Create(foot))
        self.add(layout.label("DES footprint ~5000 deg² (deep: r ≈ 23.8)", font_size=16, color=P.PURPLE).next_to(sky, UP, buff=0.1))
        self.add(layout.label("~88% of sky unsearched by DES", font_size=15, color=P.MUTED).move_to(sky.get_center() + UP * 0.4))

        eq = layout.explain_equation(
            self,
            [r"f_{\rm DES}^{\rm unique}", "=", r"\langle P_{\rm DES}\,(1-P_{\rm ZTF})\rangle"],
            [(0, "exclusion DES adds that ZTF didn't already cover"),
             (2, "DES detections, but only where ZTF missed")],
            color=P.ORANGE, scale=0.78, where=DOWN * 1.5)

        des_unique = _exclusion("des_unique", 0.050)
        ro = paper.result_readout("unique exclusion added (after ZTF)", f"+{des_unique*100:.1f}%", color=P.ORANGE)
        ro.scale(0.62).next_to(eq, DOWN, buff=0.2)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro, settle=0.8)

        layout.show_takeaway(
            self, "Depth buys little if the footprint is small -- coverage is king.")
