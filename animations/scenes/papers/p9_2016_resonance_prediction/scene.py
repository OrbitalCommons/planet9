"""Resonance prediction -- where P9 would park TNOs.

If Planet Nine confines objects in mean-motion resonances, they should pile up at
specific semimajor axes (the n:1, n:2 locations). Reproduced in
p9-2016-resonance-prediction.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, Line, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2016-resonance-prediction"


class ResonancePrediction2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Where to expect the pile-ups")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        axis = Line([-6, -0.5, 0], [6, -0.5, 0], color=P.FG, stroke_width=2)
        self.play(Create(axis), FadeIn(Text("semi-major axis", font_size=18, color=P.MUTED).next_to(axis, DOWN, buff=0.2)))
        a9 = 460.0
        res = [("3:2", a9 * (2 / 3) ** (2 / 3)), ("2:1", a9 * 0.5 ** (2 / 3)),
               ("3:1", a9 * (1 / 3) ** (2 / 3)), ("5:2", a9 * (2 / 5) ** (2 / 3))]
        amin, amax = 200, 460
        for name, a in res:
            x = -6 + (a - amin) / (amax - amin) * 11
            tick = Line([x, -0.8, 0], [x, 0.4, 0], color=P.PURPLE, stroke_width=3)
            self.play(Create(tick), FadeIn(Text(name, font_size=15, color=P.PURPLE).next_to(tick, UP, buff=0.05)), run_time=0.5)
        self.play(FadeIn(layout.takeaway("A testable prediction: TNOs should cluster at P9's resonant distances.")))
        self.wait(0.9)
