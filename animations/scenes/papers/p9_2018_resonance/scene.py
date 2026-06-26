"""Resonance capture probability.

Catalogues which P9 resonances most readily capture distant TNOs and with what
probability -- the statistical backbone of the resonance picture. Reproduced in
p9-2018-resonance (resonance catalog, capture probability).
"""
import numpy as np
from manim import Create, DOWN, FadeIn, LEFT, Rectangle, RIGHT, Scene, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2018-resonance"


class Resonance2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Which resonances catch the most?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        data = [("2:1", 0.9), ("3:1", 0.7), ("5:2", 0.55), ("3:2", 0.45), ("4:1", 0.35)]
        rows = VGroup()
        for name, p in data:
            track = Rectangle(width=6.5, height=0.5, color=P.MUTED, stroke_width=1.5)
            fill = Rectangle(width=6.5 * p, height=0.5, stroke_width=0, color=P.PURPLE).set_fill(P.PURPLE, opacity=0.5)
            fill.align_to(track, LEFT).align_to(track, UP)
            lab = layout.label(f"{name}   {p*100:.0f}%", font_size=16, color=P.FG).next_to(track, LEFT, buff=0.2)
            rows.add(VGroup(lab, track, fill))
        rows.arrange(DOWN, buff=0.35).shift(DOWN * 0.3)
        for r in rows:
            self.play(FadeIn(r[0]), Create(r[1]), Create(r[2]), run_time=0.5)

        eq = layout.equation_card(r"P_{\rm capture}(p{:}q)").scale(0.75).to_corner(UP + RIGHT, buff=0.6).shift(DOWN * 0.7)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Capture is most likely into the strong low-order resonances.")
