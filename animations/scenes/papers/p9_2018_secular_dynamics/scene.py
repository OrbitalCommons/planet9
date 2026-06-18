"""Li, Hadden, Payne & Holman (2018) -- secular dynamics of TNOs and Planet Nine.

In the eccentricity-vector plane (k, h) = e(cosΔϖ, sinΔϖ), P9's secular forcing
shifts the equilibrium off-centre (a forced eccentricity) and orbits circulate /
librate around it. Reproduced in p9-2018-secular-dynamics (free precession rate).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, Dot, FadeIn, LEFT, ParametricFunction, Scene, Text, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2018-secular-dynamics"


class SecularDynamics2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The eccentricity-vector picture")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[-0.8, 0.8, 0.4], y_range=[-0.8, 0.8, 0.4], x_length=5.6, y_length=5.6,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.3)
        xl = Text("k = e·cos Δϖ", font_size=16, color=P.FG).next_to(ax, DOWN, buff=0.15)
        yl = Text("h = e·sin Δϖ", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        # forced equilibrium offset along -k (anti-aligned), level curves around it
        kf = -0.32
        forced = Dot(ax.c2p(kf, 0), color=P.BLUE, radius=0.07)
        flbl = Text("forced eccentricity", font_size=15, color=P.BLUE).next_to(forced, DOWN, buff=0.1)
        self.play(FadeIn(forced), FadeIn(flbl))
        loops = VGroup(*[
            ParametricFunction(lambda t, r=r: ax.c2p(kf + r * np.cos(t), r * np.sin(t)),
                               t_range=[0, 2 * np.pi], color=P.TEAL).set_stroke(opacity=0.8)
            for r in (0.12, 0.24, 0.36)])
        self.play(Create(loops, lag_ratio=0.2), run_time=1.8)
        self.play(FadeIn(layout.takeaway(
            "P9's torque sets a forced eccentricity; orbits cycle around it, anti-aligned.")))
        self.wait(0.9)
