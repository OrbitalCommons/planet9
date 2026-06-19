"""Gomes, Deienno & Morbidelli (2016) -- the planetary system's tilt from P9.

A complementary framing of the obliquity: P9 tilts the giant planets' orbital
plane away from the Sun's equator, accumulating to ~6 deg. Reproduced in
p9-2016-obliquity-gomes.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, Line, Scene, Text, UP, Write, rate_functions, ValueTracker, always_redraw
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2016-obliquity-gomes"


class ObliquityGomes2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Tilting the whole planetary plane")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        eq = Line([-4.5, 0, 0], [4.5, 0, 0], color=P.SUN, stroke_width=3)
        self.play(Create(eq), FadeIn(Text("solar equator", font_size=15, color=P.SUN).next_to(eq, DOWN, buff=0.1)))
        tilt = ValueTracker(0.0)
        plane = always_redraw(lambda: Line(
            -4.0 * np.array([np.cos(np.deg2rad(tilt.get_value())), np.sin(np.deg2rad(tilt.get_value())), 0]),
            4.0 * np.array([np.cos(np.deg2rad(tilt.get_value())), np.sin(np.deg2rad(tilt.get_value())), 0]),
            color=P.BLUE, stroke_width=3))
        deg = always_redraw(lambda: Text(f"{tilt.get_value():.1f}°", font_size=22, color=P.BLUE).to_corner(UP, buff=1.4).shift([3.5, 0, 0]))
        self.add(plane, deg, Text("planets' plane", font_size=15, color=P.BLUE).to_edge(DOWN, buff=1.6))
        self.play(tilt.animate.set_value(6.0), run_time=3.5, rate_func=rate_functions.smooth)
        self.play(FadeIn(layout.takeaway("Same 6° puzzle, same culprit -- whether you tilt the Sun or the planets.")))
        self.wait(0.9)
