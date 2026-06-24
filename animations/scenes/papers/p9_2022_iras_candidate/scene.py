"""(2022) -- an IRAS far-IR Planet Nine candidate.

An earlier IRAS-based candidate search: flag moving far-IR sources consistent
with a cold, distant body's flux and motion -- the precursor to the IRAS+AKARI
pair search. Reproduced in p9-2022-iras-candidate.
"""
from manim import Create, DOWN, FadeIn, FadeOut, Scene, UP, Write, Dot
import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2022-iras-candidate"


class IrasCandidate2022(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Combing the IRAS archive")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        bg = widgets.star_field(n=90, seed=1983, x=(-6, 6), y=(-2.2, 2.0), radius=0.02, opacity=0.5)
        self.play(FadeIn(bg, lag_ratio=0.004))
        cand = Dot([1.0, 0.3, 0], radius=0.12, color=P.RED)
        self.play(FadeIn(cand))
        clbl = layout.label("IRAS 60 µm candidate", font_size=16, color=P.RED).next_to(cand, UP, buff=0.15)
        self.play(FadeIn(clbl))
        timing.hold_to_read(self, clbl, settle=0.5)

        eq = layout.explain_equation(
            self,
            [r"F_{60}", r"\Rightarrow", "d", r"\sim", r"\text{few}\times10^2\ \mathrm{AU}"],
            [(0, "measured 60-micron far-IR flux"),
             (2, "implied distance to the body"),
             (4, "a few hundred AU for a cold planet")],
            color=P.ORANGE, scale=0.8, where=UP * 1.7)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "A tantalising single-epoch flag -- needing a second epoch to confirm motion.")
