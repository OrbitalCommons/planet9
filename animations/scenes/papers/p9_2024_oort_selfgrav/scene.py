"""Batygin & Nesvorny (2024) -- self-gravitational dynamics in the inner Oort cloud.

The inner Oort cloud's own self-gravity drives collective apsidal precession,
shaping the detached population independent of (or together with) Planet Nine.
Reproduced in p9-2024-oort-selfgrav.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, RIGHT, Scene, Text, UP, VGroup, Write, rate_functions, ValueTracker, always_redraw
import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2024-oort-selfgrav"


class OortSelfgrav2024(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A self-stirring cloud")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)
        varpis0 = np.linspace(0, 2 * np.pi, 7, endpoint=False)
        drift = ValueTracker(0.0)
        live = always_redraw(lambda: VGroup(*[
            orbits.ellipse_orbit(2.6, 0.7, color=P.GREEN, varpi=v + drift.get_value(), stroke_width=1.6, opacity=0.8)
            for v in varpis0]))
        self.add(live)
        clbl = Text("collective self-gravity → coherent precession", font_size=17, color=P.TEAL).to_edge(UP, buff=1.5)
        self.play(FadeIn(clbl))
        timing.hold_to_read(self, clbl, settle=0.4)
        self.play(drift.animate.set_value(np.pi), run_time=4.0, rate_func=rate_functions.linear)

        eq = layout.equation_card(
            r"\dot\varpi_{\rm self} \propto \dfrac{M_{\rm IOC}}{M_\odot}\,n"
        ).scale(0.8).to_corner(UP + RIGHT, buff=0.5).shift(DOWN * 0.7)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "The disk can stir itself -- another way to get coherent apsidal motion.")
