"""Bailey, Batygin & Brown (2016) -- Solar Obliquity Induced by Planet Nine.

An inclined, massive P9 torques the planets over billions of years, slowly
tilting the Sun's spin axis away from the planets' mean plane -- naturally
producing the observed ~6 deg solar obliquity. Reproduced in p9-2016-obliquity.
"""
import numpy as np
from manim import (
    Angle,
    Arc,
    Create,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    MathTex,
    Scene,
    Text,
    UP,
    Write,
    rate_functions,
    ValueTracker,
    always_redraw,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2016-obliquity"


class Obliquity2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Why is the Sun tilted?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = Dot(np.zeros(3), radius=0.55, color=P.SUN)
        plane = Line(np.array([-4.5, -1.6, 0]), np.array([4.5, -1.6, 0]),
                     color=P.MUTED, stroke_width=2)
        plane_lbl = Text("planets' mean plane", font_size=16, color=P.MUTED).next_to(plane, DOWN, buff=0.1)
        self.play(FadeIn(sun), Create(plane), FadeIn(plane_lbl))

        # spin axis starts vertical, precesses to ~6 deg over the Sun's age
        vert = Line(np.zeros(3), UP * 2.2, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        tilt = ValueTracker(0.0)
        axis = always_redraw(lambda: Line(
            np.zeros(3), 2.2 * np.array([np.sin(np.deg2rad(tilt.get_value())),
                                         np.cos(np.deg2rad(tilt.get_value())), 0]),
            color=P.ORANGE, stroke_width=4))
        self.add(vert, axis)
        deg = always_redraw(lambda: Text(f"{tilt.get_value():.1f}°", font_size=24, color=P.ORANGE)
                            .next_to(axis.get_end(), UP, buff=0.1))
        self.add(deg)

        # a far inclined P9 orbit doing the torquing
        p9 = orbits.ellipse_orbit(3.2, 0.4, color=P.BLUE, varpi=np.deg2rad(30))
        p9 = p9.copy().apply_matrix(np.array([[1, 0, 0], [0, 0.4, 0], [0, 0, 1]])).shift(UP * 1.2)
        p9_lbl = Text("inclined Planet Nine", font_size=15, color=P.BLUE).next_to(p9, UP, buff=0.05)
        self.play(Create(p9), FadeIn(p9_lbl))

        prec = layout.equation(
            r"\dot{\boldsymbol{s}} \propto \boldsymbol{s}\times\boldsymbol{\ell}_9",
            scale=0.85)
        prec.next_to(p9_lbl, UP, buff=0.2)
        layout.show_equation(self, prec)

        self.play(tilt.animate.set_value(6.0), run_time=4.0, rate_func=rate_functions.smooth)
        self.wait(0.3)
        ro = paper.result_readout("solar obliquity produced", "≈ 6°", color=P.ORANGE)
        ro.scale(0.8).to_edge(DOWN, buff=1.4)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro)
        self.play(FadeOut(ro), FadeOut(prec))

        eq = layout.equation_card(r"\varepsilon_\odot \approx 6^\circ", scale=1.0)
        eq.to_edge(DOWN, buff=1.4)
        layout.show_equation(self, eq)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self,
            "The Sun's 6° tilt -- a long-standing puzzle -- falls out for free.")
