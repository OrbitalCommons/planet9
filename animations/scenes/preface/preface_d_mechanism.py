"""Preface D -- How one eccentric planet sculpts the distant belt.

Scenes: P08Rings, P08Island, P08Herding.

All three use the coplanar, orbit-averaged (secular) model of Batygin & Brown
(2016) computed in crates/p9-anim-data/src/preface/d_mechanism.rs: Planet Nine
(10 Earth masses, a = 700 AU, e = 0.6) and the giant planets smeared into rings,
Hamilton's equations integrated for belt orbits. Planet Nine's perihelion points
along +x throughout, so an orbit's orientation on screen *is* its Δϖ.
"""
import numpy as np
from manim import (
    Arc,
    Arrow,
    Circle,
    Create,
    DOWN,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    GrowArrow,
    LEFT,
    LaggedStart,
    Line,
    MathTex,
    Polygon,
    RIGHT,
    Scene,
    Text,
    TracedPath,
    UP,
    VGroup,
    VMobject,
    ValueTracker,
    Write,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, timing, widgets


def _data():
    return dataio.section("preface")["d_mechanism"]


def _header(scene, badge, title):
    scene.add(layout.concept_badge(badge))
    t = Text(title, color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.55)
    scene.play(Write(t), run_time=1.0)
    return t


def _say(scene, text, hold=True, settle=0.6):
    cap = layout.caption(text)
    scene.play(FadeIn(cap, shift=UP * 0.1), run_time=0.5)
    if hold:
        timing.hold_to_read(scene, cap, settle=settle)
    return cap


def _unsay(scene, *caps):
    scene.play(*[FadeOut(c) for c in caps], run_time=0.4)


def _u(angle):
    return np.array([np.cos(angle), np.sin(angle), 0.0])


def _orbit(a, e, varpi, scale, sun, color, stroke_width=2.0, opacity=1.0):
    """A Kepler ellipse at true (linear) scale with the Sun-focus at ``sun``."""
    return orbits.ellipse_orbit(a * scale, e, color=color, varpi=varpi, stroke_width=stroke_width,
                                opacity=opacity).shift(sun)


def _plot_axes(x_range, y_range, x_len, y_len, center, xticks, yticks, xlabel, ylabel):
    """Film axes with legible hand-placed tick labels. Returns (ax, labels)."""
    ax = widgets.axes(x_range, y_range, x_length=x_len, y_length=y_len, shift_down=0)
    ax.move_to(center)
    xt = VGroup(*[layout.label(t, font_size=18, color=P.FG)
                  .next_to(ax.c2p(v, y_range[0]), DOWN, buff=0.12) for v, t in xticks])
    yt = VGroup(*[layout.label(t, font_size=18, color=P.FG)
                  .next_to(ax.c2p(x_range[0], v), LEFT, buff=0.12) for v, t in yticks])
    parts = [xt, yt]
    if xlabel:
        parts.append(layout.label(xlabel, font_size=18, color=P.FG).next_to(xt, DOWN, buff=0.12))
    if ylabel:
        parts.append(layout.label(ylabel, font_size=18, color=P.FG).rotate(np.pi / 2)
                     .next_to(yt, LEFT, buff=0.15))
    return ax, VGroup(*parts)


def _neptune(sun, scale, d):
    return Circle(radius=d["neptune_a"] * scale, color=P.FG, stroke_width=1.5).move_to(sun)


class P08Rings(Scene):
    """Learning goal: averaged over its ~18,500-yr orbit, Planet Nine acts like a
    lopsided ring of mass; its pull on a belt orbit depends on the angle Δϖ
    between the two perihelia, and the slope of that dependence is a torque
    that changes the belt orbit's eccentricity."""

    def construct(self):
        d = _data()
        p9 = d["p9"]
        E = d["energy"]
        _header(self, "PREFACE 08", "Planet Nine's pull, averaged: a lopsided ring")

        # ---- beat 1: the planet sweeps its orbit; average it into a ring -------
        sc = 0.0037
        sun_p = np.array([-2.55, -0.05, 0.0])
        sun = orbits.sun(radius=0.08).move_to(sun_p)
        nep = _neptune(sun_p, sc, d)
        p9_orb = _orbit(p9["a"], p9["e"], 0.0, sc, sun_p, P.BLUE, 2.4)
        p9_lbl = VGroup(
            layout.label("Planet Nine", font_size=20, color=P.BLUE, weight="BOLD"),
            layout.label(f"{p9['mass_earth']:.0f} Earth masses, a = {p9['a']:.0f} AU, e = {p9['e']:.1f}",
                         font_size=16, color=P.BLUE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        p9_lbl.move_to(np.array([-6.8, 2.75, 0]), aligned_edge=LEFT)
        nep_lbl = VGroup(layout.label("Sun and Neptune's orbit", font_size=16, color=P.FG),
                         layout.label("(tiny at this scale)", font_size=16, color=P.FG)
                         ).arrange(DOWN, aligned_edge=RIGHT, buff=0.08)
        nep_lbl.next_to(sun, DOWN, buff=0.25).align_to(sun, RIGHT)
        self.play(FadeIn(sun), Create(nep), FadeIn(nep_lbl), Create(p9_orb), FadeIn(p9_lbl), run_time=1.5)

        xs, ys = np.array(d["ring"]["x"]), np.array(d["ring"]["y"])
        ring_pts = [sun_p + sc * np.array([x, y, 0]) for x, y in zip(xs, ys)]
        ring_pts = ring_pts[::2]
        step = ValueTracker(0.0)
        n_ring = len(ring_pts)
        planet = always_redraw(lambda: Dot(ring_pts[int(step.get_value()) % n_ring], radius=0.09,
                                           color=P.BLUE))
        beads = always_redraw(lambda: VGroup(*[
            Dot(ring_pts[j], radius=0.06, color=P.BLUE).set_opacity(0.85)
            for j in range(min(int(step.get_value()) + 1, n_ring))]))
        cap = _say(self, f"One lap takes ~{p9['period_yr']:,.0f} years: fast near the Sun, slow far out",
                   hold=False)
        self.add(beads, planet)
        self.play(step.animate.set_value(n_ring - 1), run_time=5.0, rate_func=rate_functions.linear)
        _unsay(self, cap)
        cap = _say(self, f"A dot every 1/{n_ring} of a lap: they crowd far out, where the planet lingers",
                   hold=False)
        self.play(FadeOut(planet), run_time=0.4)
        timing.hold_to_read(self, cap, settle=0.4)
        _unsay(self, cap)
        cap = _say(self, "Over millions of years it acts like a lopsided ring of mass, heavy far out")
        _unsay(self, cap)

        # ---- beat 2: a belt orbit, turned through every orientation ----------
        a_b, e_b = E["a"], E["e"]
        dw = ValueTracker(0.0)
        belt = always_redraw(lambda: _orbit(a_b, e_b, np.deg2rad(dw.get_value()), sc, sun_p, P.GREEN, 2.4))
        belt_apse = always_redraw(lambda: Arrow(
            sun_p, sun_p + sc * a_b * (1 - e_b) * 3.2 * _u(np.deg2rad(dw.get_value())), buff=0,
            color=P.GREEN, stroke_width=4, tip_length=0.14))
        p9_apse = Arrow(sun_p, sun_p + 0.95 * RIGHT, buff=0, color=P.BLUE, stroke_width=4, tip_length=0.14)
        ang = always_redraw(lambda: Arc(radius=0.55, start_angle=0,
                                        angle=max(np.deg2rad(dw.get_value()), 1e-3),
                                        color=P.ORANGE, stroke_width=3).shift(sun_p))
        ang_lbl = always_redraw(lambda: MathTex(r"\Delta\varpi", color=P.ORANGE).scale(0.7).move_to(
            sun_p + 0.85 * _u(np.deg2rad(max(dw.get_value(), 40.0)) / 2)))
        belt_lbl = layout.label(f"a belt orbit: a = {a_b:.0f} AU, e = {e_b:.1f}", font_size=16,
                                color=P.GREEN)
        belt_lbl.next_to(p9_lbl, DOWN, buff=0.2, aligned_edge=LEFT)

        dws = np.array(E["dvarpi_deg"])
        hs = np.array(E["h_norm"])
        tq = np.array(E["torque_norm"])
        ax1, l1 = _plot_axes([0, 360, 90], [-1, 1, 1], 5.2, 1.9, np.array([3.9, 1.45, 0]),
                             [], [], None, None)
        ax2, l2 = _plot_axes([0, 360, 90], [-1, 1, 1], 5.2, 1.9, np.array([3.9, -1.3, 0]),
                             [(v, f"{v}°") for v in (0, 90, 180, 270, 360)], [], "Δϖ", None)
        h1 = VGroup(layout.label("interaction energy", font_size=18, color=P.FG),
                    layout.label("low", font_size=15, color=P.MUTED).next_to(ax1.c2p(0, -1), LEFT, 0.1),
                    layout.label("high", font_size=15, color=P.MUTED).next_to(ax1.c2p(0, 1), LEFT, 0.1))
        h1[0].next_to(ax1, UP, buff=0.12)
        h2 = VGroup(layout.label("torque: changes e", font_size=18, color=P.FG),
                    layout.label("0", font_size=15, color=P.MUTED).next_to(ax2.c2p(0, 0), LEFT, 0.1))
        h2[0].next_to(ax2, UP, buff=0.12)
        zero2 = DashedLine(ax2.c2p(0, 0), ax2.c2p(360, 0), color=P.MUTED, stroke_width=1.2)

        def partial(ax, ys, color):
            v = dw.get_value()
            m = dws <= v
            xx = np.append(dws[m], v)
            yy = np.append(ys[m], np.interp(v, dws, ys))
            return widgets.curve(ax, xx, yy, color=color, stroke_width=3.2)

        e_curve = always_redraw(lambda: partial(ax1, hs, P.ORANGE))
        e_dot = always_redraw(lambda: Dot(ax1.c2p(dw.get_value(), np.interp(dw.get_value(), dws, hs)),
                                          radius=0.07, color=P.ORANGE))

        cap = _say(self, "Now a belt orbit, smeared the same way. Δϖ is the angle between the "
                         "two perihelia", hold=False)
        self.play(Create(belt), GrowArrow(belt_apse), GrowArrow(p9_apse), FadeIn(belt_lbl),
                  FadeIn(ang), FadeIn(ang_lbl), FadeOut(nep_lbl))
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        cap = _say(self, "Turn it all the way round and add up every pull between the two rings",
                   hold=False)
        self.play(Create(ax1), FadeIn(h1))
        self.add(e_curve, e_dot)
        self.play(dw.animate.set_value(360.0), run_time=6.0, rate_func=rate_functions.linear)
        _unsay(self, cap)
        anti = DashedLine(ax1.c2p(180, -1), ax1.c2p(180, 1), color=P.BLUE, stroke_width=1.6)
        anti_lbl = layout.label("anti-aligned", font_size=15, color=P.BLUE)
        anti_lbl.next_to(ax1.c2p(180, -0.7), RIGHT, buff=0.1)
        cap = _say(self, "The energy depends on orientation, with a turning point at 180°: "
                         "anti-aligned", hold=False)
        self.play(Create(anti), FadeIn(anti_lbl))
        timing.hold_to_read(self, cap, settle=0.4)
        _unsay(self, cap)

        # the torque is minus the slope of that curve
        t_curve = always_redraw(lambda: partial(ax2, tq, P.ORANGE))
        t_dot = always_redraw(lambda: Dot(ax2.c2p(dw.get_value(), np.interp(dw.get_value(), dws, tq)),
                                          radius=0.07, color=P.ORANGE))
        self.play(FadeOut(e_dot), run_time=0.3)
        self.remove(e_curve)
        self.add(partial(ax1, hs, P.ORANGE))
        dw.set_value(0.0)
        cap = _say(self, "Its slope is a torque: it pumps the belt orbit's eccentricity up or down",
                   hold=False)
        self.play(Create(ax2), FadeIn(l2), FadeIn(h2), Create(zero2))
        self.add(t_curve, t_dot)
        self.play(dw.animate.set_value(360.0), run_time=5.0, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)
        z0 = Dot(ax2.c2p(180, 0), radius=0.09, color=P.BLUE)
        cap = _say(self, "At 180° the torque vanishes: an anti-aligned orbit feels no net twist")
        self.play(FadeIn(z0))
        _unsay(self, cap)
        self.play(*[FadeOut(m) for m in self.mobjects[2:]])

        # ---- beat 3: Hamilton's equations, term by term -----------------------
        eq = layout.explain_equation(
            self,
            [r"\frac{d\,\Delta\varpi}{dt}", "=", r"\frac{\partial H}{\partial G}", r"\qquad",
             r"\frac{dG}{dt}", "=", r"-\frac{\partial H}{\partial\,\Delta\varpi}"],
            [
                (2, "turning rate, from the averaged energy H (giants' rings + Planet Nine)"),
                (4, "G = √(GMa(1−e²)): the orbit's angular momentum, i.e. its eccentricity"),
                (6, "changed by the torque: minus the slope you just saw"),
            ],
            scale=1.05,
            where=UP * 0.4,
        )
        self.play(FadeOut(eq))
        layout.show_takeaway(
            self, "Averaged, Planet Nine turns belt orbits and reshapes them, depending on Δϖ.")


class P08Island(Scene):
    """Learning goal: follow the computed secular paths in the (Δϖ, e) plane:
    around Δϖ = 180° they close into an island, so a belt orbit there stays
    anti-aligned forever while its perihelion rises and falls."""

    def construct(self):
        d = _data()
        p9 = d["p9"]
        port = d["portrait"]
        a_b = port["a"]
        _header(self, "PREFACE 08", "The anti-aligned island")

        # ---- right: the (Δϖ, e) map --------------------------------------------
        ax, ax_l = _plot_axes(
            [0, 360, 90], [0.3, 1, 0.1], 5.4, 4.3, np.array([3.75, 0.25, 0]),
            [(v, f"{v}°") for v in (0, 90, 180, 270, 360)],
            [(v, f"{v:.1f}") for v in (0.4, 0.6, 0.8, 1.0)],
            "orientation Δϖ", "eccentricity e")
        q_ticks = VGroup(*[
            layout.label(f"{q:.0f}", font_size=15, color=P.MUTED).next_to(ax.c2p(360, 1 - q / a_b), RIGHT,
                                                                         buff=0.1)
            for q in (30, 100, 200)])
        q_head = layout.label("q (AU)", font_size=15, color=P.MUTED).next_to(ax.c2p(360, 1.0), UP, buff=0.1)
        e_nep = port["e_neptune"]
        nep_band = Polygon(ax.c2p(0, e_nep), ax.c2p(360, e_nep), ax.c2p(360, 1), ax.c2p(0, 1),
                           stroke_width=0, color=P.RED).set_fill(P.RED, opacity=0.25)
        tag = layout.label(f"all belt orbits here: a = {a_b:.0f} AU", font_size=16, color=P.FG)
        tag.next_to(ax, UP, buff=0.12)

        # ---- left: real space, true scale --------------------------------------
        sc = 0.0034
        sun_p = np.array([-3.05, 0.2, 0.0])
        sun = orbits.sun(radius=0.08).move_to(sun_p)
        nep = _neptune(sun_p, sc, d)
        p9_orb = _orbit(p9["a"], p9["e"], 0.0, sc, sun_p, P.BLUE, 2.2)
        b9 = p9["a"] * np.sqrt(1 - p9["e"] ** 2)
        p9_lbl = layout.label("Planet Nine", font_size=18, color=P.BLUE)
        p9_lbl.next_to(sun_p + sc * np.array([-p9["a"] * p9["e"], b9, 0]), UP, buff=0.1)

        cap = _say(self, "Map every belt orbit by two numbers: its orientation Δϖ and its "
                         "eccentricity e", hold=False)
        self.play(FadeIn(sun), Create(nep), Create(p9_orb), FadeIn(p9_lbl),
                  Create(ax), FadeIn(ax_l), FadeIn(tag), FadeIn(q_ticks), FadeIn(q_head))
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        cap = _say(self, "High e means a low perihelion q: the red band is inside Neptune's orbit",
                   hold=False)
        self.play(FadeIn(nep_band))
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)

        # ---- beat 2: follow one orbit ------------------------------------------
        tracks = port["tracks"]
        feat = next(t for t in tracks if t["kind"] == "libration" and t["start_deg"] == 180.0)
        tw = np.array(feat["dvarpi_deg"])
        te = np.array(feat["e"])
        tt = np.array(feat["t_myr"])
        k = ValueTracker(0.0)

        def state():
            x = k.get_value()
            return np.interp(x, np.arange(len(tw)), tw), np.interp(x, np.arange(len(te)), te)

        dot = always_redraw(lambda: Dot(ax.c2p(state()[0] % 360, state()[1]), radius=0.08, color=P.ORANGE))
        trail = TracedPath(dot.get_center, stroke_color=P.ORANGE, stroke_width=3)
        belt = always_redraw(lambda: _orbit(a_b, state()[1], np.deg2rad(state()[0]), sc, sun_p,
                                            P.GREEN, 2.4))
        peri = always_redraw(lambda: Dot(sun_p + sc * a_b * (1 - state()[1]) * _u(np.deg2rad(state()[0])),
                                         radius=0.07, color=P.ORANGE))
        read = always_redraw(lambda: VGroup(
            layout.label(f"t = {np.interp(k.get_value(), np.arange(len(tt)), tt):,.0f} Myr",
                         font_size=18, color=P.FG),
            layout.label(f"Δϖ = {state()[0] % 360:.0f}°", font_size=18, color=P.FG),
            layout.label(f"perihelion q = {a_b * (1 - state()[1]):.0f} AU", font_size=18, color=P.ORANGE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1).move_to(np.array([-6.8, -2.3, 0]), aligned_edge=LEFT))
        cap = _say(self, "Follow one orbit, computed from Hamilton's equations...", hold=False)
        self.play(Create(belt), FadeIn(peri), FadeIn(dot), FadeIn(read))
        self.add(trail)
        timing.hold_to_read(self, cap, settle=0.1)
        _unsay(self, cap)
        cap = _say(self, "...its direction rocks back and forth about 180°, never getting away",
                   hold=False)
        self.play(k.animate.set_value(len(tw) - 1), run_time=10.0, rate_func=rate_functions.linear)
        _unsay(self, cap)
        orbits_per = feat["period_myr"] * 1e6 / p9["period_yr"]
        cap = _say(self, f"One swing: {feat['period_myr']:.0f} Myr, about {round(orbits_per, -3):,.0f} "
                         "laps of Planet Nine", hold=False)
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        cap = _say(self, f"And its perihelion breathes: {feat['q_min']:.0f} AU up to "
                         f"{feat['q_max']:.0f} AU, lifted far from Neptune")
        _unsay(self, cap)

        # ---- beat 3: the whole family of paths ----------------------------------
        def track_curve(t, color, width):
            m = VMobject(color=color, stroke_width=width)
            pts = [ax.c2p(w % 360, e) for w, e in zip(t["dvarpi_deg"], t["e"])]
            # split where Δϖ wraps so circulating paths do not jump across the plot
            segs = [[pts[0]]]
            for p0, p1, w0, w1 in zip(pts[:-1], pts[1:], t["dvarpi_deg"][:-1], t["dvarpi_deg"][1:]):
                if int(w0 // 360) != int(w1 // 360):
                    segs.append([])
                segs[-1].append(p1)
            g = VGroup()
            for s in segs:
                if len(s) > 1:
                    seg = VMobject(color=color, stroke_width=width)
                    seg.set_points_as_corners(s)
                    g.add(seg)
            return g

        island = VGroup(*[track_curve(t, P.ORANGE, 2.4) for t in tracks
                          if t["kind"] == "libration" and t["start_deg"] == 180.0
                          and t["e0"] in (0.7, 0.76, 0.9)])
        circ = VGroup(*[track_curve(t, P.FG, 2.0) for t in tracks
                        if t["kind"] == "circulation" and t["e0"] >= 0.7])
        isl_lbl = layout.label("island: trapped", font_size=18, color=P.ORANGE)
        isl_lbl.move_to(ax.c2p(180, 0.52))
        circ_lbl = layout.label("circulating", font_size=18, color=P.FG)
        circ_lbl.move_to(ax.c2p(255, 0.37))
        cap = _say(self, "Every curve is one computed path. Around 180° they close into an island",
                   hold=False)
        self.play(LaggedStart(*[Create(c) for c in island], lag_ratio=0.3), FadeIn(isl_lbl), run_time=2.5)
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        cap = _say(self, "Outside it, orbits circulate, and their perihelia sink into Neptune's zone",
                   hold=False)
        self.play(LaggedStart(*[Create(c) for c in circ], lag_ratio=0.3), FadeIn(circ_lbl), run_time=2.5)
        timing.hold_to_read(self, cap, settle=0.4)
        _unsay(self, cap)
        cap = _say(self, "Neptune then scatters them away. The trapped ones survive")
        _unsay(self, cap, read)
        layout.show_takeaway(
            self, "Orbits trapped around anti-alignment stay clustered, with perihelia lifted.")


class P08Herding(Scene):
    """Learning goal: start a scattered disk pointing every which way and run
    the computed secular evolution for 4 Gyr: the survivors are the orbits
    trapped anti-aligned with Planet Nine, perihelia lifted; the real distant
    objects show the same pattern."""

    def construct(self):
        d = _data()
        p9 = d["p9"]
        sw = d["swarm"]
        parts = sw["particles"]
        tgrid = np.array(sw["t_myr"])
        n_t = len(tgrid)
        _header(self, "PREFACE 08", "Four billion years of herding")

        # ---- left: real space ----------------------------------------------------
        sc = 0.0027
        sun_p = np.array([-3.0, 0.05, 0.0])
        sun = orbits.sun(radius=0.07).move_to(sun_p)
        nep = _neptune(sun_p, sc, d)
        p9_orb = _orbit(p9["a"], p9["e"], 0.0, sc, sun_p, P.BLUE, 2.6)
        p9_lbl = layout.label("Planet Nine", font_size=18, color=P.BLUE)
        b9 = p9["a"] * np.sqrt(1 - p9["e"] ** 2)
        p9_lbl.next_to(sun_p + sc * np.array([-p9["a"] * p9["e"], b9, 0]), UP, buff=0.1)

        clock = ValueTracker(0.0)

        def idx():
            return np.interp(clock.get_value(), tgrid, np.arange(n_t))

        def st(p):
            x = idx()
            return (np.interp(x, np.arange(n_t), p["dvarpi_deg"]), np.interp(x, np.arange(n_t), p["e"]))

        def fate(p):
            """alive / just removed (flash red) / gone."""
            r = p["removed_at"]
            if r is None or idx() < r:
                return "alive"
            return "flash" if idx() < r + 1.5 else "gone"

        def swarm():
            g = VGroup()
            for p in parts:
                f = fate(p)
                if f == "gone":
                    continue
                w, e = st(p) if f == "alive" else (p["dvarpi_deg"][p["removed_at"]], p["e"][p["removed_at"]])
                col, op = (P.GREEN, 0.55) if f == "alive" else (P.RED, 0.9)
                g.add(_orbit(p["a"], e, np.deg2rad(w), sc, sun_p, col, 1.3, op))
            return g

        swarm_m = always_redraw(swarm)

        # ---- right: every orbit as a dot in (Δϖ, q) -----------------------------
        ax, ax_l = _plot_axes(
            [0, 360, 90], [0, 250, 50], 5.0, 3.5, np.array([3.9, 0.8, 0]),
            [(v, f"{v}°") for v in (0, 90, 180, 270, 360)],
            [(v, f"{v}") for v in (0, 50, 100, 150, 200, 250)],
            "orientation Δϖ", "perihelion q (AU)")
        band = Polygon(ax.c2p(0, 0), ax.c2p(360, 0), ax.c2p(360, sw["q_neptune"]), ax.c2p(0, sw["q_neptune"]),
                       stroke_width=0).set_fill(P.RED, opacity=0.3)
        band_lbl = layout.label("inside Neptune", font_size=15, color=P.RED)
        band_lbl.next_to(ax.c2p(360, sw["q_neptune"] / 2), LEFT, buff=0.1)

        def dots():
            g = VGroup()
            for p in parts:
                f = fate(p)
                if f == "gone":
                    continue
                if f == "alive":
                    w, e = st(p)
                    q = min(p["a"] * (1 - e), 249.0)
                    g.add(Dot(ax.c2p(w % 360, q), radius=0.06, color=P.GREEN))
                else:
                    r = p["removed_at"]
                    g.add(Dot(ax.c2p(p["dvarpi_deg"][r] % 360, max(p["a"] * (1 - p["e"][r]), 2.0)),
                              radius=0.06, color=P.RED))
            return g

        dots_m = always_redraw(dots)

        def alive_n():
            return sum(1 for p in parts if fate(p) == "alive")

        read = always_redraw(lambda: VGroup(
            layout.label(f"t = {clock.get_value():,.0f} Myr", font_size=20, color=P.FG),
            layout.label(f"orbits left: {alive_n()} of {len(parts)}", font_size=20, color=P.GREEN),
        ).arrange(RIGHT, buff=0.5).move_to(np.array([3.9, -1.95, 0])))
        note = layout.label("flat, orbit-averaged model; Neptune removes q < 30 AU", font_size=14,
                            color=P.MUTED)
        note.move_to(np.array([-6.9, -2.8, 0]), aligned_edge=LEFT)

        cap = _say(self, f"A scattered disk: {len(parts)} orbits, a = 250-600 AU, perihelia 32-50 AU, "
                         "pointing every which way", hold=False)
        self.play(FadeIn(sun), Create(nep), Create(p9_orb), FadeIn(p9_lbl), FadeIn(note))
        self.play(FadeIn(swarm_m), Create(ax), FadeIn(ax_l), FadeIn(band), FadeIn(band_lbl),
                  FadeIn(dots_m), FadeIn(read))
        timing.hold_to_read(self, cap, settle=0.4)
        _unsay(self, cap)

        cap = _say(self, "Run the computed evolution: orbits that reach Neptune are scattered (red)",
                   hold=False)
        self.play(clock.animate.set_value(800.0), run_time=10.0, rate_func=rate_functions.linear)
        _unsay(self, cap)
        cap = _say(self, "The survivors rock about Δϖ = 180°: opposite Planet Nine's perihelion",
                   hold=False)
        self.play(clock.animate.set_value(tgrid[-1]), run_time=8.0, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)

        rb = sw["r_bar"][-1]
        mdw = sw["mean_dvarpi_deg"][-1]
        facts = VGroup(
            layout.label(f"clustering R̄ = {rb:.2f}, centred at Δϖ = {mdw:.0f}°", font_size=18,
                         color=P.ORANGE),
            layout.label(f"median perihelion: {sw['median_q_start']:.0f} AU → {sw['median_q_end']:.0f} AU",
                         font_size=18, color=P.ORANGE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12).move_to(np.array([3.9, -2.05, 0]))
        cap = _say(self, "Clustered opposite Planet Nine's perihelion, with perihelia lifted off Neptune",
                   hold=False)
        read.clear_updaters()
        self.play(FadeOut(read), FadeIn(facts))
        timing.hold_to_read(self, cap, facts, settle=0.5)
        _unsay(self, cap)

        # the real distant objects, measured against this Planet Nine
        real = VGroup(*[
            Arrow(sun_p, sun_p + 1.6 * _u(np.deg2rad(w)), buff=0, color=P.GREEN, stroke_width=4,
                  tip_length=0.15)
            for w in d["etno_dvarpi_deg"]])
        cap = _say(self, "Now the ten real distant objects, measured from this Planet Nine's orbit",
                   hold=False)
        swarm_m.clear_updaters()
        dots_m.clear_updaters()
        self.play(swarm_m.animate.set_stroke(opacity=0.2), LaggedStart(*[GrowArrow(r) for r in real],
                                                                  lag_ratio=0.1), run_time=1.8)
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        cap = _say(self, f"They point the same way (R̄ = {d['etno_r_bar']:.2f}): the pattern Planet Nine "
                         "was proposed to explain")
        _unsay(self, cap, note)
        layout.show_takeaway(
            self, "An eccentric Planet Nine traps belt orbits anti-aligned: the clustering clue.")
