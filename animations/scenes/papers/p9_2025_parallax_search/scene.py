"""Socas-Navarro & Trujillo (2025) -- a targeted parallax search.

A nearby body shows a large reflex parallax as Earth orbits the Sun; pairs of
images months apart can isolate that signature even without long arcs. Reproduced
in p9-2025-parallax-search.
"""
import numpy as np
from manim import (
    Arrow, Create, DOWN, Dot, FadeIn, FadeOut, Scene, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2025-parallax-search"


class ParallaxSearch2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Letting Earth do the work")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun(radius=0.18)
        e1 = Dot([-2.0, -0.3, 0], radius=0.09, color=P.GREEN)
        e2 = Dot([2.0, -0.3, 0], radius=0.09, color=P.GREEN)
        self.play(FadeIn(sun), FadeIn(e1), FadeIn(e2))
        self.add(layout.label("Earth, 6 months apart", font_size=15, color=P.GREEN).next_to(sun, DOWN, buff=0.4))

        p9 = Dot([0.6, 2.0, 0], radius=0.11, color=P.BLUE)
        l1 = Arrow(e1.get_center(), p9.get_center(), color=P.MUTED, buff=0.1, stroke_width=2)
        l2 = Arrow(e2.get_center(), p9.get_center(), color=P.MUTED, buff=0.1, stroke_width=2)
        self.play(FadeIn(p9), Create(l1), Create(l2))
        note = layout.label("apparent shift = parallax → distance", font_size=16, color=P.TEAL).to_edge(DOWN, buff=1.9)
        self.play(FadeIn(note))
        timing.hold_to_read(self, note, settle=0.5)

        eq = layout.explain_equation(
            self,
            ["p", "=", r"\dfrac{1\ \mathrm{AU}}{d}"],
            [(0, "the parallax angle from Earth's baseline"),
             (2, "one AU divided by the body's distance d")],
            color=P.TEAL, scale=0.8, where=UP * 1.7)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Two epochs and Earth's baseline can pin a nearby planet's distance directly.")
