"""Preface A -- Foundations: the hook, Newton's ellipse, and Kepler's laws.

Scenes: P00Hook, P01Ellipse, P02KeplerLaws.
"""
import numpy as np
from manim import (
    Arc,
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    GrowArrow,
    Indicate,
    LEFT,
    Line,
    MathTex,
    Polygon,
    RIGHT,
    Scene,
    UP,
    Unwrite,
    ValueTracker,
    VGroup,
    Write,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import layout, orbits, timing


class P00Hook(Scene):
    def construct(self):
        title = layout.title_card(
            "Planet Nine",
            "the case for a world we have never seen",
        )
        self.play(Write(title[0]), run_time=1.2)
        self.play(FadeIn(title[1], shift=UP * 0.2))
        self.wait(0.4)
        self.play(title.animate.scale(0.5).to_edge(UP))

        sun = orbits.sun()
        # a faint, huge Planet Nine orbit + a few clustered ETNOs
        p9 = orbits.ellipse_orbit(3.6, 0.4, color=P.BLUE, varpi=np.deg2rad(20))
        swarm = orbits.etno_swarm(
            [(2.4, 0.7, np.deg2rad(a)) for a in (180, 196, 168, 205, 190)],
            color=P.GREEN,
        )
        self.play(FadeIn(sun))
        self.play(Create(swarm, lag_ratio=0.2), run_time=2.0)
        self.play(Create(p9), run_time=1.5)
        self.play(Indicate(p9, color=P.BLUE, scale_factor=1.05))
        self.wait(0.5)
        self.play(*[FadeOut(m) for m in [title, sun, p9, swarm]])


class P01Ellipse(Scene):
    def construct(self):
        badge = layout.concept_badge("PREFACE 01")
        self.add(badge)

        law = layout.equation_card(r"\vec{F} = -\,\frac{G M m}{r^2}\,\hat{r}").scale(1.1)
        layout.show_equation(self, law)
        sub = layout.caption("one inverse-square pull -> the bound orbit is an ELLIPSE")
        self.play(law.animate.scale(0.6).to_edge(UP, buff=0.9), FadeIn(sub))
        timing.hold_to_read(self, sub, settle=0.6)
        self.play(FadeOut(sub))

        a, e = 3.3, 0.6
        sun = orbits.sun()
        ell = orbits.ellipse_orbit(a, e, color=P.BLUE)
        self.play(Create(ell), run_time=1.6)
        self.play(FadeIn(sun))
        sun_lbl = MathTex(r"\odot", color=P.SUN).scale(0.7).next_to(sun, DOWN, buff=0.12)
        self.play(FadeIn(sun_lbl))

        # major axis + semi-major axis a
        peri = orbits.orbit_point(a, e, 0.0)
        apo = orbits.orbit_point(a, e, np.pi)
        center = (peri + apo) / 2.0
        major = Line(apo, peri, color=P.MUTED, stroke_width=2)
        a_line = Line(center, peri, color=P.TEAL, stroke_width=4)
        a_lbl = layout.equation("a", color=P.TEAL).scale(0.8).next_to(a_line, UP, buff=0.12)
        self.play(Create(major))
        self.play(Create(a_line), FadeIn(a_lbl))

        # perihelion / aphelion markers
        q_dot = orbits.body_dot(a, e, 0.0, color=P.GREEN)
        Q_dot = orbits.body_dot(a, e, np.pi, color=P.RED)
        q_lbl = layout.equation(r"q = a(1-e)", color=P.GREEN).scale(0.6).next_to(q_dot, RIGHT, buff=0.15)
        Q_lbl = layout.equation(r"Q = a(1+e)", color=P.RED).scale(0.6).next_to(Q_dot, LEFT, buff=0.15)
        self.play(FadeIn(q_dot, q_lbl), FadeIn(Q_dot, Q_lbl))
        timing.hold_to_read(self, q_lbl, Q_lbl, settle=0.6)

        # the orbit equation r(nu) as a brief hero beat
        r_eq = layout.equation_card(r"r(\nu) = \dfrac{a(1-e^2)}{1 + e\cos\nu}").scale(0.8)
        r_eq.to_edge(DOWN, buff=0.9)
        layout.show_equation(self, r_eq)
        self.play(FadeOut(r_eq))

        e_lbl = layout.equation(r"e:\ \text{eccentricity (how stretched)}", color=P.FG).scale(0.6)
        e_lbl.to_edge(DOWN, buff=0.6)
        self.play(FadeIn(e_lbl))
        timing.hold_to_read(self, e_lbl, settle=0.4)
        self.play(FadeOut(e_lbl))

        layout.show_takeaway(self, "A bound orbit is an ellipse; the Sun sits at one focus.")


class P02KeplerLaws(Scene):
    def construct(self):
        self.add(layout.concept_badge("PREFACE 02"))
        a, e = 3.2, 0.65
        sun = orbits.sun()
        ell = orbits.ellipse_orbit(a, e, color=P.BLUE)
        self.add(ell, sun)

        # --- Kepler II: equal areas in equal times ---
        k2 = layout.caption("Kepler II: the Sun-body line sweeps equal AREAS in equal TIMES")
        self.play(FadeIn(k2))

        # equal-area sectors: narrow near perihelion, wide near aphelion
        def sector(nu0, nu1):
            pts = [np.zeros(3)]
            for nu in np.linspace(nu0, nu1, 30):
                pts.append(orbits.orbit_point(a, e, nu))
            return Polygon(*pts, color=P.GREEN, fill_opacity=0.35, stroke_width=1)

        peri_sec = sector(np.deg2rad(-26), np.deg2rad(26))
        apo_sec = sector(np.deg2rad(150), np.deg2rad(210))
        self.play(FadeIn(peri_sec), FadeIn(apo_sec))

        area_eq = layout.equation(
            r"\dfrac{dA}{dt} = \tfrac12 r^2\dot\nu = \text{const}", color=P.GREEN
        ).scale(0.7).to_corner(UP + RIGHT, buff=0.5)
        self.play(Write(area_eq))
        timing.hold_to_read(self, area_eq, settle=0.4)

        M = ValueTracker(0.0)
        body = always_redraw(
            lambda: orbits.body_dot(a, e, orbits.nu_from_mean_anomaly(e, M.get_value()), color=P.GREEN)
        )
        radius = always_redraw(
            lambda: Line(np.zeros(3), orbits.orbit_point(a, e, orbits.nu_from_mean_anomaly(e, M.get_value())),
                         color=P.MUTED, stroke_width=2)
        )
        self.add(radius, body)
        self.play(M.animate.set_value(2 * np.pi), run_time=5.0, rate_func=rate_functions.linear)
        self.wait(0.2)
        self.play(FadeOut(peri_sec), FadeOut(apo_sec), FadeOut(radius), FadeOut(body),
                  FadeOut(k2), FadeOut(area_eq))

        # --- vis-viva: speed at any point on the orbit ---
        vv = layout.equation(r"v^2 = GM\!\left(\dfrac{2}{r}-\dfrac{1}{a}\right)", color=P.TEAL)
        vv.to_edge(UP, buff=0.9)
        self.play(Write(vv))
        timing.hold_to_read(self, vv, settle=0.5)
        self.play(FadeOut(vv))

        # --- Kepler III ---
        k3 = layout.equation_card(r"P^2 \;=\; a^3 \quad (\text{years, AU})")
        k3.to_edge(UP, buff=1.0)
        layout.show_equation(self, k3)
        layout.show_takeaway(
            self, "Farther = slower. Planet Nine's orbit takes ~10,000 yr -- we cannot wait one.")
