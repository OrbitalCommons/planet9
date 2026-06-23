"""Preface B -- Orienting an orbit, and the geography of the outer system.

Scenes: P03Elements, P04Geography.
"""
import numpy as np
from manim import (
    Arc,
    Arrow,
    Circle,
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    GrowArrow,
    LEFT,
    Line,
    RIGHT,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, orbits, timing


class P03Elements(Scene):
    """The orbital elements -- especially longitude of perihelion (the clue)."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 03"))
        intro = layout.caption("Six numbers fix an orbit in space.")
        self.play(FadeIn(intro))
        self.wait(0.4)
        self.play(FadeOut(intro))

        sun = orbits.sun()
        # the ecliptic plane, seen face-on, as a faint reference disk
        ecl = Circle(radius=3.4, color=P.MUTED, stroke_width=1.5).set_stroke(opacity=0.5)
        ecl_lbl = Text("ecliptic plane", color=P.MUTED, font_size=18).next_to(ecl, DOWN, buff=0.05)
        ref = Arrow(np.zeros(3), RIGHT * 3.4, color=P.FG, buff=0, stroke_width=2)
        ref_lbl = Text("reference direction", color=P.FG, font_size=16).next_to(ref, RIGHT, buff=0.05)
        self.play(Create(ecl), FadeIn(ecl_lbl), GrowArrow(ref), FadeIn(ref_lbl), FadeIn(sun))

        Omega = np.deg2rad(55)   # longitude of ascending node
        omega = np.deg2rad(70)   # argument of perihelion
        varpi = Omega + omega    # longitude of perihelion
        a, e = 2.6, 0.55

        # ascending node direction
        node = Line(np.zeros(3), 3.4 * np.array([np.cos(Omega), np.sin(Omega), 0]),
                    color=P.PURPLE, stroke_width=2)
        node_lbl = Text("ascending node", color=P.PURPLE, font_size=15)
        node_lbl.next_to(node.get_end(), UP, buff=0.05)
        Om_arc = Arc(radius=1.0, start_angle=0, angle=Omega, color=P.PURPLE)
        Om_lbl = layout.equation(r"\Omega", color=P.PURPLE).scale(0.7).move_to(
            1.3 * np.array([np.cos(Omega / 2), np.sin(Omega / 2), 0]))
        self.play(Create(node), FadeIn(node_lbl), Create(Om_arc), FadeIn(Om_lbl))
        self.wait(0.3)

        # the orbit, perihelion along varpi
        orbit = orbits.ellipse_orbit(a, e, color=P.BLUE, varpi=varpi)
        apse = orbits.apse_arrow(a, e, varpi, color=P.ORANGE)
        peri_lbl = Text("perihelion", color=P.ORANGE, font_size=15).next_to(apse.get_end(), UP, buff=0.05)
        self.play(Create(orbit))
        self.play(GrowArrow(apse), FadeIn(peri_lbl))

        # argument of perihelion omega (node -> perihelion)
        om_arc = Arc(radius=1.6, start_angle=Omega, angle=omega, color=P.GREEN)
        om_lbl = layout.equation(r"\omega", color=P.GREEN).scale(0.7).move_to(
            1.9 * np.array([np.cos(Omega + omega / 2), np.sin(Omega + omega / 2), 0]))
        self.play(Create(om_arc), FadeIn(om_lbl))
        self.wait(0.3)

        # longitude of perihelion varpi = Omega + omega
        var_eq = layout.equation(
            r"\varpi = \Omega + \omega", t2c={r"\Omega": P.PURPLE, r"\omega": P.GREEN}
        ).scale(0.8).to_corner(UP + RIGHT, buff=0.5)
        self.play(Write(var_eq))
        self.play(apse.animate.set_color(P.TEAL))
        timing.hold_to_read(self, var_eq, settle=0.5)

        # Kepler's equation: where the body is along the orbit at a given time
        kep = layout.equation_card(r"M = E - e\sin E").scale(0.8)
        kep.to_edge(DOWN, buff=0.9)
        layout.show_equation(self, kep)
        self.play(FadeOut(kep))

        # inclination note
        inc = layout.equation(r"i:\ \text{tilt of the orbit out of this plane}", color=P.FG).scale(0.55)
        inc.to_edge(DOWN, buff=0.9)
        self.play(FadeIn(inc))
        timing.hold_to_read(self, inc, settle=0.4)
        self.play(FadeOut(inc))

        layout.show_takeaway(
            self, "Clustering is about ϖ (and the tilt direction): which way orbits point.")


class P04Geography(Scene):
    """Scale of the outer solar system on a log distance axis."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 04"))
        title = Text("How far is 'far'?", color=P.FG, font_size=34, weight="BOLD").to_edge(UP, buff=0.7)
        self.play(Write(title))

        # log distance axis from 1 to 2000 AU
        x0, x1 = -6.0, 6.0
        lo, hi = np.log10(1.0), np.log10(2000.0)

        def xf(au):
            return x0 + (np.log10(au) - lo) / (hi - lo) * (x1 - x0)

        axis = Line(LEFT * 6, RIGHT * 6, color=P.FG, stroke_width=2).shift(DOWN * 0.5)
        ax_lbl = Text("distance from the Sun (AU, log scale)", color=P.MUTED, font_size=18)
        ax_lbl.next_to(axis, DOWN, buff=0.4)
        self.play(Create(axis), FadeIn(ax_lbl))

        marks = [
            (1, "Earth", P.GREEN), (5.2, "Jupiter", P.FG), (30, "Neptune", P.BLUE),
            (45, "Kuiper Belt", P.MUTED), (250, "Sednoids /\nETNOs", P.GREEN),
            (600, "Planet Nine?", P.TEAL),
        ]
        anims = []
        for au, name, col in marks:
            x = xf(au)
            tick = Line(UP * 0.12, DOWN * 0.12, color=col, stroke_width=3).move_to([x, -0.5, 0])
            lbl = Text(name, color=col, font_size=16).next_to(tick, UP, buff=0.15)
            val = Text(f"{au:g}", color=P.MUTED, font_size=13).next_to(tick, DOWN, buff=0.12)
            self.play(Create(tick), FadeIn(lbl), FadeIn(val), run_time=0.55)
        self.wait(0.4)

        per_eq = layout.equation_card(r"P \approx a^{3/2}\ \text{(yr, AU)}").scale(0.8)
        per_eq.next_to(title, DOWN, buff=0.3)
        layout.show_equation(self, per_eq)
        self.play(FadeOut(per_eq))

        layout.show_takeaway(
            self, "Beyond Neptune lies a vast, barely-mapped realm -- where Planet Nine would hide.")
