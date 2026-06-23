"""Preface C -- The clue and the dynamics: clustering, precession, resonance.

Scenes: P05Clustering, P06Precession, P07Resonance.
"""
import numpy as np
from manim import (
    Circle,
    Create,
    DOWN,
    DecimalNumber,
    Dot,
    FadeIn,
    FadeOut,
    LEFT,
    Line,
    RIGHT,
    Rotate,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
    always_redraw,
    rate_functions,
    ValueTracker,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, timing


class P05Clustering(Scene):
    """Random vs aligned longitudes of perihelion; the mean resultant length."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 05"))
        title = Text("The clue: the orbits point the same way",
                     color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.6)
        self.play(Write(title))

        # The ALIGNED state is the REAL six stable KBOs (Batygin & Brown 2016);
        # the random state stays illustrative as the contrast.
        objs, r_bar = dataio.real_etnos("kbo")
        if objs:
            clustered_varpi = np.array([o["varpi_rad"] for o in objs])
            n = len(objs)
        else:
            n = 9
            clustered_varpi = np.deg2rad(60) + np.random.default_rng(7).normal(0, np.deg2rad(14), n)
        random_varpi = np.random.default_rng(7).uniform(0, 2 * np.pi, n)

        sun = orbits.sun(radius=0.1)
        circle = Circle(radius=3.0, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.4)
        self.play(FadeIn(sun), Create(circle))

        def arrows(varpis):
            return orbits.apse_arrows(varpis, length=2.7, color=P.GREEN, stroke_width=3)

        # the mean resultant length R-bar: how confined the directions are
        rbar_def = layout.explain_equation(
            self,
            [r"\bar R", "=", r"\left|\dfrac{1}{N}\sum_{k=1}^{N} e^{\,i\varpi_k}\right|"],
            [
                (0, "mean resultant length: 0 = scattered, 1 = perfectly aligned"),
                (2, "average each orbit's direction ϖ as a unit vector, take the length"),
            ],
            scale=0.85,
            where=UP * 1.3,
        )
        self.play(FadeOut(rbar_def))

        arr = arrows(random_varpi)
        rbar0 = orbits.mean_resultant_length(random_varpi)
        lbl0 = layout.equation(rf"\bar{{R}} = {rbar0:.2f}\ \text{{(random)}}", color=P.RED).scale(0.7)
        lbl0.to_edge(DOWN, buff=0.9)
        self.play(Create(arr, lag_ratio=0.1), FadeIn(lbl0))
        self.wait(0.6)

        # rotate each arrow to its real (or illustrative) clustered direction
        target = arrows(clustered_varpi)
        rbar1 = r_bar if r_bar is not None else orbits.mean_resultant_length(clustered_varpi)
        lbl1 = layout.equation(rf"\bar{{R}} = {rbar1:.2f}\ \text{{(real KBOs)}}", color=P.TEAL).scale(0.7)
        lbl1.to_edge(DOWN, buff=0.9)
        self.play(*[a.animate.become(t) for a, t in zip(arr, target)], run_time=2.2)
        self.play(FadeOut(lbl0), FadeIn(lbl1))
        timing.hold_to_read(self, lbl1, settle=0.5)
        layout.show_takeaway(
            self, "R-bar near 1 means the orbits are confined -- something is herding them.")


class P06Precession(Scene):
    """Apsidal precession: the ellipse slowly rotates, tracing a rosette."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 06"))
        title = Text("Orbits are not frozen: they precess",
                     color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.6)
        self.play(Write(title))

        sun = orbits.sun()
        self.add(sun)
        a, e = 2.6, 0.6

        # faint trail of past orientations (the rosette)
        trail = VGroup(*[
            orbits.ellipse_orbit(a, e, color=P.BLUE, varpi=v, stroke_width=1, opacity=0.18)
            for v in np.linspace(0, 1.6 * np.pi, 9)
        ])
        self.play(Create(trail, lag_ratio=0.15), run_time=2.0)

        vt = ValueTracker(0.0)
        live = always_redraw(lambda: orbits.ellipse_orbit(a, e, color=P.BLUE, varpi=vt.get_value()))
        apse = always_redraw(lambda: orbits.apse_arrow(a, e, vt.get_value(), color=P.ORANGE))
        self.add(live, apse)
        self.play(vt.animate.set_value(1.6 * np.pi), run_time=4.0, rate_func=rate_functions.linear)

        eq = layout.equation(r"\dot{\varpi}\ \neq\ 0", color=P.ORANGE).scale(0.9).to_corner(UP + RIGHT, buff=0.6)
        self.play(Write(eq))
        timing.hold_to_read(self, eq, settle=0.4)

        rate_eq = layout.explain_equation(
            self,
            [r"\dot\varpi", r"\approx", r"\tfrac34", "n", r"\dfrac{m_p}{M_\odot}",
             r"\left(\dfrac{a_p}{a}\right)^{2}", r"\dfrac{1}{(1-e^2)^2}"],
            [
                (0, "how fast the orbit's apse direction precesses"),
                (3, "mean motion: faster orbits precess faster"),
                (4, "perturber's mass relative to the Sun"),
                (5, "(perturber's distance / this orbit's), squared"),
                (6, "eccentric orbits precess faster"),
            ],
            scale=0.72,
            where=DOWN * 1.4,
        )
        self.play(FadeOut(rate_eq))

        layout.show_takeaway(
            self, "Differential precession scrambles alignment -- so confinement needs a cause.")


class P07Resonance(Scene):
    """A 2:1 mean-motion resonance: repeated kicks at the same phase."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 07"))
        title = Text("Mean-motion resonance: kicks that add up",
                     color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.6)
        self.play(Write(title))

        sun = orbits.sun()
        a1 = 1.6
        a2 = a1 * 2 ** (2 / 3)  # P2 = 2 P1  ->  a2 = a1 * 2^(2/3)
        o1 = Circle(radius=a1, color=P.GREEN, stroke_width=2)
        o2 = Circle(radius=a2, color=P.BLUE, stroke_width=2)
        self.play(FadeIn(sun), Create(o1), Create(o2))
        lbl = layout.equation(r"P_{\rm out} : P_{\rm in} = 2 : 1", color=P.FG).scale(0.7).to_corner(UP + RIGHT, buff=0.6)
        self.play(Write(lbl))

        # resonant semi-major axis: where periods lock into a p:q ratio
        res_eq = layout.explain_equation(
            self,
            [r"a_{\rm res}", "=", r"a_p", r"\left(\dfrac{p}{q}\right)^{2/3}"],
            [
                (0, "the orbit size that locks into resonance"),
                (2, "perturber's semi-major axis"),
                (3, "integer period ratio p:q, raised to the 2/3 power"),
            ],
            scale=0.85,
            where=UP * 1.7,
        )
        self.play(FadeOut(res_eq))

        vt = ValueTracker(0.0)
        d1 = always_redraw(lambda: Dot(a1 * np.array([np.cos(2 * vt.get_value()), np.sin(2 * vt.get_value()), 0]),
                                       color=P.GREEN, radius=0.09))
        d2 = always_redraw(lambda: Dot(a2 * np.array([np.cos(vt.get_value()), np.sin(vt.get_value()), 0]),
                                       color=P.BLUE, radius=0.09))
        self.add(d1, d2)
        # conjunction direction (they line up at the same longitude each outer orbit)
        conj = Line(np.zeros(3), RIGHT * (a2 + 0.4), color=P.ORANGE, stroke_width=2).set_stroke(opacity=0.6)
        self.add(conj)
        self.play(vt.animate.set_value(4 * np.pi), run_time=6.0, rate_func=rate_functions.linear)
        layout.show_takeaway(
            self, "When periods form a simple ratio, perturbations repeat and accumulate.")
