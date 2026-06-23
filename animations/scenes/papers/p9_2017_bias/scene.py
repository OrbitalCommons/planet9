"""Brown (2017) -- Observational bias and the clustering of distant KBOs.

A careful treatment of each survey's detection probability: even after debiasing,
a residual clustering signal remains -- the counter-argument to pure-bias claims.
Reproduced in p9-2017-bias (bias_function).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2017-bias"


class Bias2017(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Debiasing the sky")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 360, 90], y_range=[0, 1.05, 0.5], x_length=9.3, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = Text("longitude of perihelion ϖ (deg)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
        self.play(Create(ax), FadeIn(xl))

        # survey sensitivity vs ϖ (the bias), and the debiased residual excess
        sens = ax.plot(lambda v: 0.5 + 0.35 * np.cos(np.deg2rad(v - 60)), x_range=[0, 360, 2], color=P.PURPLE)
        excess = ax.plot(lambda v: 0.55 + 0.3 * np.cos(np.deg2rad(v - 250)), x_range=[0, 360, 2], color=P.GREEN)
        self.play(Create(sens))
        slab = Text("survey sensitivity (bias)", font_size=15, color=P.PURPLE).next_to(ax.c2p(60, 0.85), UP, buff=0.05)
        self.play(FadeIn(slab))
        timing.hold_to_read(self, slab, settle=0.6)
        self.play(Create(excess))
        elab = Text("debiased object excess", font_size=15, color=P.GREEN).next_to(ax.c2p(250, 0.85), UP, buff=0.05)
        self.play(FadeIn(elab))
        timing.hold_to_read(self, elab, settle=0.6)

        # debiasing: recover the intrinsic distribution by dividing out S(varpi)
        eq = layout.explain_equation(
            self,
            [r"N_{\rm intrinsic}(\varpi)", "=", r"\frac{N_{\rm obs}(\varpi)}{S(\varpi)}"],
            [
                (0, "the true number of objects at each direction"),
                (2, "counts observed, divided by survey sensitivity"),
            ],
            scale=0.8,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "The excess peaks where sensitivity does NOT -- so bias alone cannot explain it.")
