"""Resonance prediction -- where P9 would park TNOs.

If Planet Nine confines objects in mean-motion resonances, they should pile up at
specific semimajor axes (the n:1, n:2 locations). Reproduced in
p9-2016-resonance-prediction.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, Line, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper, dataio, timing

CRATE = "p9-2016-resonance-prediction"


class ResonancePrediction2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Where to expect the pile-ups")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        axis = Line([-6, -0.5, 0], [6, -0.5, 0], color=P.FG, stroke_width=2)
        self.play(Create(axis), FadeIn(Text("semi-major axis", font_size=18, color=P.MUTED).next_to(axis, DOWN, buff=0.2)))
        data = dataio.section("resonance")
        if data:
            # real mean-motion-resonance semimajor axes interior to P9
            a9 = float(data["a9"])
            res = [(it["label"], float(it["a_res"])) for it in data["items"]]
            a_vals = [a for _, a in res]
            amin = min(a_vals) - 20
            amax = a9
        else:
            a9 = 460.0
            res = [("3:2", a9 * (2 / 3) ** (2 / 3)), ("2:1", a9 * 0.5 ** (2 / 3)),
                   ("3:1", a9 * (1 / 3) ** (2 / 3)), ("5:2", a9 * (2 / 5) ** (2 / 3))]
            amin, amax = 200, 460
        for name, a in res:
            x = -6 + (a - amin) / (amax - amin) * 11
            tick = Line([x, -0.8, 0], [x, 0.4, 0], color=P.PURPLE, stroke_width=3)
            self.play(Create(tick), FadeIn(Text(name, font_size=15, color=P.PURPLE).next_to(tick, UP, buff=0.05)), run_time=0.5)

        eq = layout.explain_equation(
            self,
            [r"a_{\rm res}", "=", "a_9", r"\,(p/q)^{2/3}"],
            [
                (0, "where a resonant TNO parks (its semimajor axis)"),
                (2, "set by Planet Nine's own semimajor axis"),
                (3, "scaled by the period ratio p:q (Kepler's third law)"),
            ],
            scale=1.0,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "A testable prediction: TNOs should cluster at P9's resonant distances.")
