"""Iorio (arXiv:1512.05288) -- where along its orbit can the planet be?

Saturn's perihelion is pinned by Cassini ranging to a fraction of a
milliarcsecond per century. A distant planet's tide makes that perihelion
precess at a rate that falls as the cube of the planet's distance, so the bound
on Saturn's supplementary precession says how far along its eccentric orbit the
planet must currently be.

Everything drawn comes from anim.json -> papers -> p9-2016-iorio-precession:
the per-planet rates and bounds, Saturn's rate against true anomaly (the
crate's quadrupole rate for a body parked at the distance r(f)), and the
true anomalies where that rate crosses Saturn's bound.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    VGroup,
    ValueTracker,
    VMobject,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2016-iorio-precession"


def decade_label(k):
    return {-4: "0.0001", -3: "0.001", -2: "0.01", -1: "0.1", 0: "1", 1: "10", 2: "100"}[k]


def orbit_arc(a, e, f0, f1, color, stroke_width=5.0, n=90):
    """The part of a Kepler ellipse between true anomalies f0 and f1 (degrees)."""
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners([orbits.orbit_point(a, e, np.radians(f))
                             for f in np.linspace(f0, f1, n)])
    return m


class Iorio2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        orbit = d["orbit"]
        planets = d["planets"]
        sat = d["saturn"]
        f_lo, f_hi = d["allowed_from_deg"], d["allowed_to_deg"]
        pub = d["published"]

        self.add(paper.scene_header(CRATE))

        # 1. which planet can feel it: predicted precession against the bound
        lo, hi = -4, 2
        ax = widgets.axes([0, len(planets), 1], [lo, hi, 1], x_length=7.6, y_length=4.6,
                          shift_down=-0.1)
        ax.shift(LEFT * 1.9)
        ax.x_axis.set_opacity(0)
        base = Line(ax.c2p(0, lo), ax.c2p(len(planets), lo), color=P.MUTED, stroke_width=2)
        furniture = VGroup(base)
        for k in range(lo, hi + 1):
            furniture.add(layout.label(decade_label(k), font_size=14, color=P.MUTED)
                          .next_to(ax.c2p(0, k), LEFT, buff=0.12))
        ylab = layout.label("perihelion precession  (mas per century)", font_size=15)
        ylab.rotate(np.pi / 2).next_to(ax, LEFT, buff=0.85)
        furniture.add(ylab)

        names, bars, bounds = VGroup(), VGroup(), VGroup()
        for k, pl in enumerate(planets):
            x = k + 0.5
            names.add(layout.label(pl["name"], font_size=15).next_to(ax.c2p(x, lo), DOWN, buff=0.14))
            y0 = np.log10(pl["rate_aphelion_mas_cy"])
            y1 = np.log10(pl["rate_perihelion_mas_cy"])
            bar = Line(ax.c2p(x, y0), ax.c2p(x, y1), color=P.ORANGE, stroke_width=9)
            mean = Dot(ax.c2p(x, np.log10(pl["rate_averaged_mas_cy"])), radius=0.07,
                       color=P.FG).set_z_index(3)
            bars.add(VGroup(bar, mean))
            yb = np.log10(pl["bound_mas_cy"])
            bounds.add(Line(ax.c2p(x - 0.32, yb), ax.c2p(x + 0.32, yb), color=P.RED,
                            stroke_width=4))

        key = VGroup(
            VGroup(Line(LEFT * 0.2, RIGHT * 0.2, color=P.ORANGE, stroke_width=9),
                   layout.label("induced by the planet,\nperihelion to aphelion",
                                font_size=14, line_spacing=0.9)).arrange(RIGHT, buff=0.15),
            VGroup(Dot(radius=0.07, color=P.FG),
                   layout.label("orbit average", font_size=14)).arrange(RIGHT, buff=0.15),
            VGroup(Line(LEFT * 0.2, RIGHT * 0.2, color=P.RED, stroke_width=4),
                   layout.label("largest rate the\nephemeris allows", font_size=14,
                                color=P.RED, line_spacing=0.9)).arrange(RIGHT, buff=0.15),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        key.next_to(ax, RIGHT, buff=0.35).shift(UP * 1.2)

        cap = layout.caption(
            f"A {orbit['mass_earth']:.0f} M⊕ planet on the {orbit['a_au']:.0f} AU orbit "
            "makes every planet's perihelion precess", font_size=22)
        self.play(Create(ax), FadeIn(furniture), FadeIn(names), run_time=1.0)
        self.play(FadeIn(bars, lag_ratio=0.15), FadeIn(key[0]), FadeIn(key[1]), FadeIn(cap),
                  run_time=1.4)
        timing.hold_to_read(self, cap, key[0], settle=0.8)

        ratio = {pl["name"]: pl["rate_averaged_mas_cy"] / pl["bound_mas_cy"] for pl in planets}
        runner_up = max((n for n in ratio if n != "Saturn"), key=ratio.get)
        cap2 = layout.caption(
            f"Saturn, ranged by Cassini, is the bound the planet reaches first; "
            f"{runner_up} is close behind", font_size=22)
        ring = Polygon(ax.c2p(5.05, -1.55), ax.c2p(5.95, -1.55), ax.c2p(5.95, 0.85),
                       ax.c2p(5.05, 0.85), color=P.FG, stroke_width=1.6)
        self.play(FadeIn(bounds, lag_ratio=0.15), FadeIn(key[2]), FadeOut(cap), FadeIn(cap2),
                  run_time=1.2)
        self.play(Create(ring), run_time=0.6)
        timing.hold_to_read(self, cap2, key[2], settle=1.0)
        self.play(FadeOut(VGroup(ax, furniture, names, bars, bounds, key, ring, cap2)))

        # 2. Saturn's rate against where the planet is on its orbit
        f = np.array(sat["f_deg"])
        rate = np.array(sat["rate_mas_cy"])
        bound = sat["bound_mas_cy"]
        top = 4.0
        ax2, labels2 = widgets.labeled_axes(
            [0, 360, 60], [0, top, 1], x_label="true anomaly of the planet  (0° = perihelion)",
            y_label="Saturn's induced precession  (mas per century)", y_rotate=True,
            numbers=True, x_length=10.4, y_length=4.3, shift_down=-0.2, font_size=22)

        def band(a, b, color, opacity):
            p0, p1 = ax2.c2p(a, 0), ax2.c2p(b, top)
            r = Polygon(p0, [p1[0], p0[1], 0], p1, [p0[0], p1[1], 0], stroke_width=0)
            return r.set_fill(color, opacity=opacity)

        excluded = VGroup(band(0, f_lo, P.RED, 0.16), band(f_hi, 360, P.RED, 0.16))
        allowed = band(f_lo, f_hi, P.TEAL, 0.13)
        curve = widgets.curve(ax2, f, rate, color=P.ORANGE, stroke_width=3.5)
        bline = DashedLine(ax2.c2p(0, bound), ax2.c2p(360, bound), color=P.RED, stroke_width=2.5)
        blab = layout.label(f"Saturn's bound  {bound:.2f} mas per century", font_size=14,
                            color=P.RED).next_to(ax2.c2p(180, bound), UP, buff=0.1)

        cap3 = layout.caption(
            "The pull falls as distance cubed, so it depends on where the planet is now",
            font_size=22)
        self.play(Create(ax2), FadeIn(labels2), run_time=0.9)
        self.play(Create(curve), FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, settle=0.6)

        cap4 = layout.caption(
            "Too fast near perihelion: only the far side of the orbit stays under the bound",
            font_size=22)
        span = Line(ax2.c2p(f_lo, 2.3), ax2.c2p(f_hi, 2.3), color=P.TEAL, stroke_width=3)
        span_l = layout.label(f"reproduced: {f_lo:.0f}° to {f_hi:.0f}°", font_size=15,
                              color=P.TEAL).next_to(span, UP, buff=0.1)
        pspan = Line(ax2.c2p(pub["allowed_from_deg"], 1.5), ax2.c2p(pub["allowed_to_deg"], 1.5),
                     color=P.FG, stroke_width=3)
        pspan_l = layout.label(
            f"paper: {pub['allowed_from_deg']:.0f}° to {pub['allowed_to_deg']:.0f}°",
            font_size=15).next_to(pspan, UP, buff=0.1)
        self.play(Create(bline), FadeIn(blab), FadeOut(cap3), FadeIn(cap4), run_time=1.0)
        self.play(FadeIn(excluded), FadeIn(allowed), run_time=0.8)
        self.play(Create(span), FadeIn(span_l), Create(pspan), FadeIn(pspan_l), run_time=1.0)
        timing.hold_to_read(self, cap4, span_l, pspan_l, settle=1.2)
        self.play(FadeOut(VGroup(ax2, labels2, excluded, allowed, curve, bline, blab, span,
                                 span_l, pspan, pspan_l, cap4)))

        # 3. the same verdict on the orbit itself
        scale = 5.6 / (2.0 * orbit["a_au"])
        a, e = orbit["a_au"] * scale, orbit["e"]
        shift = np.array([2.2, 0.15, 0.0])
        red = orbit_arc(a, e, -(360 - f_hi), f_lo, P.RED)
        teal = orbit_arc(a, e, f_lo, f_hi, P.TEAL)
        sun = orbits.sun(radius=0.11)
        view = VGroup(red, teal, sun).shift(shift)
        ticks = VGroup()
        for fp in (pub["allowed_from_deg"], pub["allowed_to_deg"]):
            pt = orbits.orbit_point(a, e, np.radians(fp)) + shift
            out = pt - sun.get_center()
            out = out / np.linalg.norm(out)
            ticks.add(Line(pt - 0.16 * out, pt + 0.16 * out, color=P.FG, stroke_width=3))
        peri = layout.label(f"perihelion  {orbit['perihelion_au']:.0f} AU", font_size=14,
                            color=P.RED)
        peri.next_to(orbits.orbit_point(a, e, 0.0) + shift, RIGHT, buff=0.15)
        aph = layout.label(f"aphelion  {orbit['aphelion_au']:.0f} AU", font_size=14, color=P.TEAL)
        aph.next_to(orbits.orbit_point(a, e, np.pi) + shift, LEFT, buff=0.15)
        sun_l = layout.label("Sun and Saturn", font_size=13, color=P.MUTED)
        sun_l.next_to(sun, DOWN, buff=0.12)
        notes = VGroup(
            layout.label("ruled out by Saturn's perihelion", font_size=17, color=P.RED),
            layout.label(f"allowed: beyond {d['critical_distance_au']:.0f} AU", font_size=17,
                         color=P.TEAL),
            layout.label(f"white ticks: the paper's {pub['allowed_from_deg']:.0f}° and "
                         f"{pub['allowed_to_deg']:.0f}°", font_size=17),
        ).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
        notes.move_to([0, 2.3, 0]).to_edge(LEFT, buff=0.6)
        self.play(FadeIn(sun), FadeIn(sun_l), Create(red), Create(teal), run_time=1.6)
        self.play(FadeIn(peri), FadeIn(aph), FadeIn(ticks), FadeIn(notes), run_time=0.8)
        timing.hold_to_read(self, notes, settle=0.6)

        # 4. fly the planet round its orbit at Kepler speed, reading Saturn's drift
        f_tab, r_tab = f, rate
        mean = ValueTracker(np.pi)

        def true_anom_deg():
            nu = np.degrees(orbits.nu_from_mean_anomaly(e, mean.get_value())) % 360.0
            return nu

        def rate_now():
            return float(np.interp(true_anom_deg(), f_tab, r_tab))

        planet = always_redraw(lambda: Dot(
            orbits.orbit_point(a, e, np.radians(true_anom_deg())) + shift, radius=0.11,
            color=P.BLUE).set_z_index(6))
        tether = always_redraw(lambda: DashedLine(
            sun.get_center(), orbits.orbit_point(a, e, np.radians(true_anom_deg())) + shift,
            color=P.MUTED, stroke_width=1.5))

        def meter():
            v = rate_now()
            col = P.RED if v > bound else P.TEAL
            head = layout.label("Saturn's perihelion drift now", font_size=16, color=P.MUTED)
            val = layout.label(f"{v:.2f} mas per century", font_size=26, color=col,
                               weight="BOLD")
            verdict = layout.label("too fast: excluded" if v > bound else "under the bound: allowed",
                                   font_size=16, color=col)
            g = VGroup(head, val, verdict).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
            return g.move_to([0, -1.0, 0]).to_edge(LEFT, buff=0.6)

        readout = always_redraw(meter)
        bound_note = layout.label(f"Saturn's bound: {bound:.2f}", font_size=15, color=P.RED)
        bound_note.move_to([0, -1.9, 0]).to_edge(LEFT, buff=0.6)
        cap5 = layout.caption(
            "Moving at Kepler speed, the planet dawdles far out where Saturn cannot feel it",
            font_size=22)
        self.play(FadeIn(planet), FadeIn(tether), FadeIn(readout), FadeIn(bound_note),
                  FadeIn(cap5), run_time=0.8)
        self.play(mean.animate.set_value(3 * np.pi), run_time=9.0, rate_func=lambda t: t)
        share = d["allowed_time_fraction"]
        cap6 = layout.caption(
            f"It spends {100 * share:.0f}% of every orbit on the allowed arc: "
            "the planet is most likely there", font_size=22)
        self.play(FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=0.8)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "Saturn says the planet is on the far half of its orbit, not near perihelion.")
