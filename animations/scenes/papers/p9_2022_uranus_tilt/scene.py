"""Lu & Laughlin (2022) -- primordial spin-orbit resonance origin of Uranus' tilt.

A slow sweep of a spin-orbit resonance in the early solar system can tip a
planet's spin axis far over -- the same obliquity machinery invoked for the Sun
under Planet Nine. Reproduced in p9-2022-uranus-tilt.
"""
import numpy as np
from manim import (
    Create, DOWN, Dot, FadeIn, FadeOut, Line, Scene, Text, UP, Write, rate_functions,
    ValueTracker, always_redraw,
)

import p9_manim as P
from p9_manim import layout, paper, timing

CRATE = "p9-2022-uranus-tilt"


class UranusTilt2022(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "How to tip a planet on its side")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        planet = Dot(np.zeros(3), radius=0.7, color=P.TEAL)
        self.play(FadeIn(planet))
        vert = Line(DOWN * 1.8, UP * 1.8, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        tilt = ValueTracker(0.0)
        axis = always_redraw(lambda: Line(
            -1.8 * np.array([np.sin(np.deg2rad(tilt.get_value())), np.cos(np.deg2rad(tilt.get_value())), 0]),
            1.8 * np.array([np.sin(np.deg2rad(tilt.get_value())), np.cos(np.deg2rad(tilt.get_value())), 0]),
            color=P.ORANGE, stroke_width=4))
        deg = always_redraw(lambda: Text(f"obliquity {tilt.get_value():.0f}°", font_size=22, color=P.ORANGE)
                            .to_edge(UP, buff=1.6))
        self.add(vert, axis, deg)
        self.play(tilt.animate.set_value(98.0), run_time=4.0, rate_func=rate_functions.smooth)
        self.wait(0.5)

        eq = layout.explain_equation(
            self,
            [r"\varepsilon:", r"0^\circ", r"\to", r"98^\circ"],
            [
                (0, "a planet's axial tilt (obliquity)"),
                (1, "starts upright"),
                (3, "a resonance can topple it sideways, like Uranus"),
            ],
            scale=0.95,
            where=DOWN * 1.2)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self,
            "Resonant obliquity sweeps are real -- the same physics the Sun's 6° tilt invokes.")
