"""Phan et al. (2025) -- IRAS + AKARI far-IR candidate search.

A cold Planet Nine would have moved measurably between the IRAS (1983) and AKARI
(2006) far-IR all-sky surveys. Cross-matching sources with consistent fluxes and
the right ~23-year motion yields a short candidate list. Reproduced in
p9-2025-iras-akari.
"""
from manim import (
    Arrow,
    Create,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    LEFT,
    RIGHT,
    Scene,
    Text,
    UP,
    Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2025-iras-akari"


class IrasAkari2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A 23-year baseline in the far-IR")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        bg = widgets.star_field(n=80, seed=1983, x=(-6, 6), y=(-2.2, 2.0), radius=0.02, opacity=0.5)
        self.play(FadeIn(bg, lag_ratio=0.005))

        # IRAS 1983 detection, AKARI 2006 detection, shifted by proper motion
        p1 = Dot(LEFT * 2.0 + DOWN * 0.2, radius=0.10, color=P.RED)
        l1 = Text("IRAS 1983 (60 µm)", font_size=16, color=P.RED).next_to(p1, DOWN, buff=0.2)
        self.play(FadeIn(p1), FadeIn(l1))
        p2 = Dot(LEFT * 2.0 + RIGHT * 1.6 + UP * 0.25, radius=0.10, color=P.ORANGE)
        l2 = Text("AKARI 2006 (90 µm)", font_size=16, color=P.ORANGE).next_to(p2, UP, buff=0.2)
        arr = Arrow(p1.get_center(), p2.get_center(), color=P.TEAL, buff=0.12, stroke_width=3)
        self.play(FadeIn(p2), FadeIn(l2), Create(arr))
        mv = Text("proper motion over 23 yr matches a ~500–700 AU body",
                  font_size=16, color=P.TEAL).next_to(arr, DOWN, buff=0.3)
        self.play(FadeIn(mv))
        timing.hold_to_read(self, mv, settle=0.6)

        eq = layout.explain_equation(
            self,
            [r"\mu", "=", r"\dfrac{\Delta\theta}{\Delta t}", r"\ (23\ \mathrm{yr})"],
            [(0, "proper motion across the sky"),
             (2, "angular shift divided by elapsed time"),
             (3, "measured over the IRAS-to-AKARI baseline")],
            color=P.TEAL, scale=0.8, where=UP * 1.7)
        self.play(FadeOut(eq))
        ro = paper.result_readout("candidates after flux + motion cuts", "1 good candidate", color=P.ORANGE).scale(0.7)
        ro.next_to(mv, DOWN, buff=0.45)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro, settle=0.8)

        layout.show_takeaway(
            self, "Thermal surveys give an independent channel -- candidates still need follow-up.")
