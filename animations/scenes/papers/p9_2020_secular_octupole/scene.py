"""Köhne & Batygin (2020) -- retrograde Jupiter Trojans and high-inclination TNOs.

Asteroid 514107 Ka'epaoka'awela shares Jupiter's orbit but travels it
backwards (i = 163°). Köhne & Batygin run its clones back 100 Myr and find them
on polar orbits beyond Neptune, in the reservoir of high-inclination objects
that Planet Nine feeds. The reproduction crate (p9-2020-secular-octupole) does
not integrate the clones; it computes the octupole-order secular Hamiltonian
through which an eccentric Planet Nine reshapes distant orbits, and that phase
portrait is what the last beat shows (anim.json -> papers ->
p9-2020-secular-octupole). The Trojan's orbit and the clone result are the
published values, labelled as such.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arc,
    Axes,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    MoveAlongPath,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2020-secular-octupole"


def dial(centre, radius, color=P.MUTED):
    """A 0-180° inclination dial: prograde on the right, retrograde on the left."""
    arc = Arc(radius=radius, start_angle=0, angle=np.pi, color=color, stroke_width=2)
    arc.shift(centre)
    ticks = VGroup()
    for deg in (0, 45, 90, 135, 180):
        th = np.deg2rad(deg)
        u = np.array([np.cos(th), np.sin(th), 0.0])
        ticks.add(Line(centre + 0.93 * radius * u, centre + 1.07 * radius * u,
                       color=color, stroke_width=2))
        ticks.add(layout.label(f"{deg}°", font_size=16, color=color)
                  .move_to(centre + 1.22 * radius * u))
    base = Line(centre + LEFT * radius, centre + RIGHT * radius, color=color, stroke_width=1)
    return VGroup(arc, ticks, base)


def needle(centre, radius, deg, color):
    th = np.deg2rad(deg)
    tip = centre + radius * np.array([np.cos(th), np.sin(th), 0.0])
    return VGroup(Line(centre, tip, color=color, stroke_width=5), Dot(tip, radius=0.07, color=color))


def band(centre, radius, lo, hi, color):
    """Filled wedge between inclinations lo and hi (degrees)."""
    from manim import AnnularSector

    return AnnularSector(inner_radius=0, outer_radius=radius, angle=np.deg2rad(hi - lo),
                         start_angle=np.deg2rad(lo), color=color, fill_opacity=0.25,
                         stroke_width=0).shift(centre)


class EdgeAxes(Axes):
    """Axes that meet at the lower-left corner even when a range spans zero."""

    @staticmethod
    def _origin_shift(axis_range):
        return axis_range[0]


def edge_axes(x_range, y_range, x_label, y_label, x_length, y_length):
    ax = EdgeAxes(x_range=x_range, y_range=y_range, x_length=x_length, y_length=y_length,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16,
                               "numbers_to_exclude": []})
    ax.add_coordinates()
    xl = layout.label(x_label, font_size=18).next_to(ax, DOWN, buff=0.12)
    yl = layout.label(y_label, font_size=15).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.12)
    return ax, VGroup(xl, yl)


class SecularOctupole2020(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        tro = d["trojan"]
        e_jup = dataio.body("Jupiter")["e"]
        self.add(paper.scene_header(CRATE))

        # 1. a Trojan going the wrong way round
        scale = 2.3 / tro["jupiter_a_au"]
        sun = orbits.sun(radius=0.13)
        jup_orbit = orbits.ellipse_orbit(tro["jupiter_a_au"] * scale, e_jup, color=P.FG,
                                         stroke_width=2.5, varpi=0.3)
        tro_orbit = orbits.ellipse_orbit(tro["a_au"] * scale, tro["e"], color=P.GREEN,
                                         stroke_width=2.5, varpi=np.deg2rad(200))
        VGroup(sun, jup_orbit, tro_orbit).shift(UP * 0.1)
        phase = ValueTracker(0.0)

        def body(a_au, e, varpi, sign, m0, color, r):
            def build():
                nu = orbits.nu_from_mean_anomaly(e, m0 + sign * phase.get_value())
                return Dot(orbits.orbit_point(a_au * scale, e, nu, varpi) + UP * 0.1,
                           radius=r, color=color).set_z_index(4)
            return always_redraw(build)

        jup = body(tro["jupiter_a_au"], e_jup, 0.3, 1.0, 0.0, P.FG, 0.13)
        kae = body(tro["a_au"], tro["e"], np.deg2rad(200), -1.0, 2.2, P.GREEN, 0.08)
        jl = layout.label(f"Jupiter, a = {tro['jupiter_a_au']:.1f} AU", font_size=18, color=P.FG)
        jl.move_to([4.6, 1.2, 0])
        kl = layout.label(f"{tro['name']}\na = {tro['a_au']:.2f} AU, i = {tro['i_deg']:.0f}°",
                          font_size=18, color=P.GREEN)
        kl.move_to([-4.9, -1.0, 0])
        arrow_j = layout.label("↺ forwards", font_size=20, color=P.FG).next_to(jl, DOWN, buff=0.15)
        arrow_k = layout.label("↻ backwards", font_size=20, color=P.GREEN).next_to(kl, UP,
                                                                                buff=0.15)
        self.play(FadeIn(sun), Create(jup_orbit), Create(tro_orbit), run_time=1.2)
        self.add(jup, kae)
        self.play(FadeIn(jl), FadeIn(kl), FadeIn(arrow_j), FadeIn(arrow_k))
        cap = layout.caption("It shares Jupiter's orbit (the 1:−1 resonance) "
                             "but travels it the wrong way round", font_size=22)
        self.play(FadeIn(cap), phase.animate.set_value(2 * np.pi), run_time=6.0, rate_func=linear)
        timing.hold_to_read(self, cap, settle=0.2)
        for m in (jup, kae):
            m.clear_updaters()
        self.play(FadeOut(VGroup(sun, jup_orbit, tro_orbit, jup, kae, jl, kl, arrow_j,
                                 arrow_k, cap)))

        # 2. where it came from: the paper's backward integration
        c = np.array([-2.6, -1.6, 0.0])
        r = 2.6
        dl = dial(c, r)
        now = needle(c, r * 0.95, tro["i_deg"], P.GREEN)
        now_l = layout.label("today\n163°", font_size=18, color=P.GREEN)
        now_l.move_to(c + 1.45 * r * np.array([np.cos(np.deg2rad(166)), np.sin(np.deg2rad(166)), 0]))
        pro = layout.label("prograde", font_size=16, color=P.MUTED).move_to(c + [r * 0.62, -0.25, 0])
        ret = layout.label("retrograde", font_size=16, color=P.MUTED).move_to(c + [-r * 0.62, -0.25, 0])
        self.play(Create(dl), FadeIn(pro), FadeIn(ret))
        self.play(Create(now), FadeIn(now_l))
        wedge = band(c, r * 0.95, 90, 135, P.ORANGE)
        text = VGroup(
            layout.label("Paper: its clones, run back 100 Myr,", font_size=20),
            layout.label("tilt down to i = 90°–135°", font_size=20, color=P.ORANGE),
            layout.label("while a and e grow: polar orbits", font_size=20),
            layout.label("beyond Neptune", font_size=20),
        ).arrange(DOWN, buff=0.16, aligned_edge=LEFT).move_to([3.6, 0.6, 0])
        self.play(FadeIn(text[:2]))
        incl = ValueTracker(tro["i_deg"])
        past = always_redraw(lambda: needle(c, r * 0.95, incl.get_value(), P.ORANGE))
        self.add(past)
        self.play(FadeIn(wedge), incl.animate.set_value(112.5), run_time=2.5)
        past.clear_updaters()
        self.play(FadeIn(text[2:]))
        cap2 = layout.caption("The same reservoir of high-inclination orbits "
                              "that Planet Nine is predicted to fill", font_size=22)
        self.play(FadeIn(cap2))
        timing.hold_to_read(self, text, cap2, settle=0.6)
        self.play(FadeOut(VGroup(dl, now, now_l, pro, ret, wedge, past, text, cap2)))

        # 3. the lever Planet Nine pulls: the octupole
        ax, labels = edge_axes(
            [-180, 180, 90], [0, 0.9, 0.3], "apse angle to Planet Nine  Δϖ (deg)",
            "eccentricity e", x_length=7.4, y_length=4.4)
        VGroup(ax, labels).move_to([-1.9, 0.15, 0])
        start = d["start"]
        a0 = d["a_prototype_au"]
        quad = DashedLine(ax.c2p(-180, start["e"]), ax.c2p(180, start["e"]),
                          color=P.MUTED, stroke_width=3)
        quad_l = layout.label("quadrupole only: e stays put,\nthe apse just circulates",
                              font_size=16, color=P.MUTED)
        quad_l.next_to(ax.c2p(-180, start["e"]), UP + RIGHT, buff=0.12)
        side = VGroup(
            layout.label(f"orbit at a = {a0:.0f} AU", font_size=18),
            layout.label(f"Planet Nine: {d['p9']['mass_earth']:.0f} M⊕,", font_size=18,
                         color=P.BLUE),
            layout.label(f"a = {d['p9']['a_au']:.0f} AU, e = {d['p9']['e']:.1f}", font_size=18,
                         color=P.BLUE),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([4.6, 1.9, 0])
        self.play(Create(ax), FadeIn(labels), FadeIn(side))
        self.play(Create(quad), FadeIn(quad_l))
        cap3 = layout.caption("To lowest order, Planet Nine's pull cannot change "
                              "a distant orbit's shape", font_size=22)
        self.play(FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.4)

        curves = VGroup()
        lib = None
        for cv in d["curves"]:
            if not cv["e"]:
                continue
            is_start = abs(cv["e0"] - start["e"]) < 1e-9
            col = P.ORANGE if cv["motion"] == "libration" else P.MUTED
            m = widgets.curve(ax, cv["dvarpi_deg"], cv["e"], color=col,
                              stroke_width=3.5 if is_start else 1.6)
            if cv["motion"] == "libration":
                m.add_line_to(m.get_start())
            if not is_start:
                m.set_stroke(opacity=0.55)
            curves.add(m)
            if is_start:
                lib = m
        eps = d["epsilon_oct"]
        cap4 = layout.caption(f"Add the octupole (strength ε = {eps:.2f} at {a0:.0f} AU): "
                              "orbits now trace closed loops", font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), FadeOut(quad_l), run_time=0.6)
        self.play(quad.animate.set_opacity(0.35), Create(curves, lag_ratio=0.1), run_time=2.2)
        dot = Dot(lib.get_start(), radius=0.09, color=P.GREEN).set_z_index(5)
        self.add(dot)
        self.play(MoveAlongPath(dot, lib), run_time=4.0, rate_func=linear)
        loop_e = next(cv["e"] for cv in d["curves"] if abs(cv["e0"] - start["e"]) < 1e-9)
        e_lo, e_hi = min(loop_e), max(loop_e)
        res = VGroup(
            layout.label(f"apse held within ±{d['libration_max_deg']:.0f}°", font_size=18,
                         color=P.ORANGE),
            layout.label(f"e swings {e_lo:.2f} → {e_hi:.2f}", font_size=18, color=P.ORANGE),
            layout.label(f"perihelion {a0 * (1 - e_lo):.0f} → {a0 * (1 - e_hi):.0f} AU",
                         font_size=18, color=P.ORANGE),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([4.6, -0.4, 0])
        self.play(FadeIn(res))
        cap5 = layout.caption("The octupole locks the apse and pumps e: "
                              "the lever that moves orbits in and out", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), run_time=0.6)
        self.play(MoveAlongPath(dot, lib), run_time=4.0, rate_func=linear)
        timing.hold_to_read(self, cap5, res, settle=0.8)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Planet Nine reshapes distant orbits; the Trojan traces back to them.")
