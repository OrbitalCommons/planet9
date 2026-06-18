"""Batygin & Morbidelli (2017) -- Dynamical Evolution Induced by Planet Nine.

The secular (orbit-averaged) phase space has a stable island in which a distant
KBO's apsidal angle Δϖ relative to P9 librates -- i.e. stays confined -- rather
than circulating freely. Reproduced in p9-2017-dynamics (phase-space/Hamiltonian).
"""
import numpy as np
from manim import (
    Axes,
    Create,
    DOWN,
    Dot,
    FadeIn,
    LEFT,
    MathTex,
    ParametricFunction,
    Scene,
    Text,
    UP,
    UR,
    VGroup,
    Write,
    rate_functions,
    ValueTracker,
    always_redraw,
)

import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2017-dynamics"


class Dynamics2017(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A trap in phase space")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 360, 90], y_range=[0, 0.8, 0.2], x_length=9.5, y_length=4.2,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = MathTex(r"\Delta\varpi\ \text{(deg, relative to P9)}", color=P.FG).scale(0.6).next_to(ax, DOWN, buff=0.2)
        yl = Text("eccentricity", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        # libration loops (closed contours) centred at Delta-varpi = 180
        def loop(scale):
            return ParametricFunction(
                lambda t: ax.c2p(180 + 70 * scale * np.cos(t), 0.42 + 0.26 * scale * np.sin(t)),
                t_range=[0, 2 * np.pi], color=P.TEAL)
        loops = VGroup(*[loop(s) for s in (0.35, 0.65, 1.0)])
        self.play(Create(loops, lag_ratio=0.2), run_time=1.8)

        # circulating trajectories (top): Delta-varpi sweeps all values
        circ = VGroup(*[ax.plot(lambda x, c=c: c, x_range=[0, 360, 5], color=P.MUTED).set_stroke(opacity=0.5)
                        for c in (0.66, 0.74)])
        self.play(Create(circ))

        center = Dot(ax.c2p(180, 0.42), color=P.BLUE, radius=0.06)
        self.play(FadeIn(center))
        # a test particle librating around the centre
        t = ValueTracker(0.0)
        body = always_redraw(lambda: Dot(
            ax.c2p(180 + 60 * np.cos(t.get_value()), 0.42 + 0.22 * np.sin(t.get_value())),
            color=P.GREEN, radius=0.07))
        self.add(body)
        self.play(t.animate.set_value(2.2 * np.pi), run_time=4.0, rate_func=rate_functions.linear)

        lbl = Text("libration: Δϖ stays confined", font_size=18, color=P.TEAL).to_corner(UR, buff=0.5)
        self.play(Write(lbl))
        self.play(FadeIn(layout.takeaway(
            "Inside the island, orbits stay anti-aligned with P9 -- that is the clustering.")))
        self.wait(0.9)
