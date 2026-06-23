"""Iorio (2026) -- a distended Planet Nine and Saturn's perihelion precession.

If Planet Nine is extended/distended rather than a point mass, its perturbation of
Saturn's perihelion precession changes -- tightening the planetary-ephemeris
bound. Reproduced in p9-2026-iorio-precession.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2026-iorio-precession"


class Iorio2026(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Saturn feels the difference")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[300, 800, 100], y_range=[0, 4, 1], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax),
                  FadeIn(Text("Planet Nine distance (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("Saturn precession effect", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        point = ax.plot(lambda d: 40.0 / (d / 250.0) ** 3, x_range=[320, 800, 5], color=P.MUTED)
        disten = ax.plot(lambda d: 52.0 / (d / 250.0) ** 3, x_range=[320, 800, 5], color=P.ORANGE)
        pm_lbl = Text("point mass", font_size=14, color=P.MUTED).next_to(ax.c2p(700, 40 / (700 / 250) ** 3), UP, buff=0.05)
        self.play(Create(point), FadeIn(pm_lbl))
        d_lbl = Text("distended P9", font_size=14, color=P.ORANGE).next_to(ax.c2p(450, 52 / (450 / 250) ** 3), UP, buff=0.05)
        self.play(Create(disten), FadeIn(d_lbl))
        timing.hold_to_read(self, pm_lbl, d_lbl)

        eq = layout.explain_equation(
            self,
            [r"\dot\varpi_{\rm Saturn}^{\rm distended}"],
            [
                (0, "Saturn's precession if P9's mass is spread out, not a point"),
            ],
            scale=0.95,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Even P9's mass distribution leaves a (tiny) fingerprint on Saturn.")
