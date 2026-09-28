"""Sefilian & Touma (2018) -- shepherding in a self-gravitating disc.

A massive, eccentric disc of trans-Neptunian objects forces the apsides of the
orbits embedded in it. The scene draws the crate's disc (its own rings), then
the two computed results: how fast the disc turns an orbit compared with the
giant planets, and the secular phase plane that decides which orbits the disc
holds. The crate's linear model confines orbits aligned with the disc and only
at low eccentricity; the paper finds anti-aligned confinement of eccentric
orbits. Both are shown. Data: anim.json -> papers -> p9-2019-selfgrav-disk.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Annulus,
    Axes,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    MoveAlongPath,
    Polygon,
    Scene,
    VGroup,
    VMobject,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2019-selfgrav-disk"


class LogAxes(VGroup):
    """Log-log axes between arbitrary limits, with labelled ticks at the given
    values. ``p(x, y)`` maps data to the scene; ``curve`` draws a polyline."""

    def __init__(self, x_lim, y_lim, x_ticks, y_ticks, x_label, y_label,
                 x_length=9.6, y_length=4.4, centre=(0.3, 0.35, 0), y_fmt=None):
        super().__init__()
        self.x0, self.y0 = np.log10(x_lim[0]), np.log10(y_lim[0])
        self.ax = Axes(x_range=[0, np.log10(x_lim[1]) - self.x0, 10],
                       y_range=[0, np.log10(y_lim[1]) - self.y0, 10],
                       x_length=x_length, y_length=y_length,
                       axis_config={"color": P.MUTED, "include_tip": False,
                                    "include_ticks": False})
        self.ax.move_to(centre)
        self.add(self.ax)
        for v in x_ticks:
            at = self.p(v, y_lim[0])
            self.add(Line(at, at + DOWN * 0.08, color=P.MUTED, stroke_width=1.5))
            self.add(layout.label(f"{v:g}", font_size=15, color=P.MUTED)
                     .next_to(at, DOWN, buff=0.14))
        widest = 0.0
        for v in y_ticks:
            at = self.p(x_lim[0], v)
            self.add(Line(at, at + LEFT * 0.08, color=P.MUTED, stroke_width=1.5))
            text = y_fmt(v) if y_fmt else f"{v:g}"
            lab = layout.label(text, font_size=15, color=P.MUTED).next_to(at, LEFT, buff=0.14)
            widest = max(widest, lab.width)
            self.add(lab)
        self.add(layout.label(x_label, font_size=17).next_to(self.ax, DOWN, buff=0.45))
        self.add(layout.label(y_label, font_size=17).rotate(np.pi / 2)
                 .next_to(self.ax, LEFT, buff=0.3 + widest))

    def p(self, x, y):
        return self.ax.c2p(np.log10(x) - self.x0, np.log10(y) - self.y0)

    def curve(self, xs, ys, color, stroke_width=3.0):
        m = VMobject(color=color, stroke_width=stroke_width)
        m.set_points_as_corners([self.p(x, y) for x, y in zip(xs, ys)])
        return m

    def band(self, x0, x1, y0, y1, color, opacity=0.13):
        a, b = self.p(x0, y0), self.p(x1, y1)
        r = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0, color=color)
        return r.set_fill(color, opacity=opacity)


def period_label(yr):
    """Tick text for a period in years: Myr below a Gyr, Gyr above."""
    return f"{yr / 1e6:g} Myr" if yr < 1e9 else f"{yr / 1e9:g} Gyr"


def polyline(points, color, stroke_width=2.2, opacity=1.0):
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners(points)
    return m.set_stroke(opacity=opacity)


class SelfgravDisk2019(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        disk = d["disk"]
        a_test = d["a_test_au"]

        self.add(paper.scene_header(CRATE))

        # 1. the disc itself: the crate's confocal, apsidally aligned rings
        centre = np.array([-2.6, 0.15, 0.0])
        scale = 2.55 / disk["a_out_au"]
        sig = np.array([r["sigma_earth_per_au2"] for r in disk["rings"]])
        rings = VGroup()
        for r, s in zip(disk["rings"], sig):
            weight = np.log10(s / sig.min()) / np.log10(sig.max() / sig.min())
            rings.add(orbits.ellipse_orbit(
                r["a_au"] * scale, r["e"], color=P.ORANGE, varpi=np.radians(r["varpi_deg"]),
                stroke_width=1.2 + 2.2 * weight, opacity=0.35 + 0.55 * weight).shift(centre))
        sun = orbits.sun(radius=0.07).move_to(centre)
        apse = Line(centre, centre + RIGHT * disk["a_out_au"] * (1 - disk["e"]) * scale,
                    color=P.ORANGE, stroke_width=1.5).set_stroke(opacity=0.7)
        facts = VGroup(
            layout.label("the disc", font_size=20, color=P.ORANGE, weight="BOLD"),
            layout.label(f"mass  {disk['mass_earth']:.0f} M⊕", font_size=18),
            layout.label(f"from {disk['a_in_au']:.0f} to {disk['a_out_au']:.0f} AU", font_size=18),
            layout.label(f"ring eccentricity  {disk['e']:.1f}", font_size=18),
            layout.label("every ring shares one apsidal line", font_size=18),
            layout.label("surface density falls as 1/a²  (line weight)", font_size=16,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT).move_to([3.6, 0.4, 0])
        cap = layout.caption("No planet: a massive, lopsided disc of small bodies beyond Neptune",
                             font_size=22)
        self.play(FadeIn(sun), LaggedStart(*[Create(r) for r in rings], lag_ratio=0.08),
                  FadeIn(cap), run_time=2.0)
        self.play(Create(apse), FadeIn(facts), run_time=0.8)
        timing.hold_to_read(self, cap, facts, settle=0.8)
        self.play(FadeOut(VGroup(rings, sun, apse, facts, cap)), run_time=0.7)

        # 2. computed: who turns the orbit faster, the disc or the planets?
        radii = np.array(d["radii_au"])
        ax = LogAxes((60, 700), (1e6, 1e11), [60, 100, 200, 400, 700],
                     [1e6, 1e7, 1e8, 1e9, 1e10, 1e11],
                     "semi-major axis of the orbit  (AU)", "one turn of the apsidal line",
                     y_length=4.3, centre=(0.1, 0.4, 0), y_fmt=period_label)
        heavy = ax.curve(radii, d["disk_period_yr"], P.ORANGE)
        light = ax.curve(radii, d["light_disk_period_yr"], P.ORANGE, stroke_width=2.0)
        light.set_stroke(opacity=0.55)
        planets = ax.curve(radii, d["planets_period_vs_a_yr"], P.PURPLE)
        age = DashedLine(ax.p(60, d["age_yr"]), ax.p(700, d["age_yr"]), color=P.RED,
                         stroke_width=2)
        lo, hi = d["paper_period_range_yr"]
        band = ax.band(60, 700, lo, hi, P.TEAL, opacity=0.10)
        tags = VGroup(
            layout.label(f"{disk['mass_earth']:.0f} M⊕ disc", font_size=15, color=P.ORANGE)
            .next_to(ax.p(700, d["disk_period_yr"][-1]), RIGHT, buff=0.1),
            layout.label(f"{d['light_disk_earth']:.0f} M⊕ disc", font_size=15, color=P.ORANGE)
            .next_to(ax.p(700, d["light_disk_period_yr"][-1]), RIGHT, buff=0.1),
            layout.label("giant planets", font_size=15, color=P.PURPLE)
            .next_to(ax.p(700, d["planets_period_vs_a_yr"][-1]), RIGHT, buff=0.1),
            layout.label("age of the Solar System", font_size=14, color=P.RED)
            .next_to(ax.p(60, d["age_yr"]), UP + RIGHT, buff=0.08),
            layout.label("paper: 100–1000 Myr", font_size=14, color=P.TEAL)
            .next_to(ax.p(60, lo), UP + RIGHT, buff=0.08),
        )
        cap2 = layout.caption(
            f"Computed for an e = {d['e_test']:.1f} orbit: how long one turn of the apsidal line takes",
            font_size=22)
        self.play(FadeIn(ax), FadeIn(cap2), run_time=0.9)
        self.play(Create(planets), FadeIn(tags[2]), run_time=1.0)
        self.play(Create(heavy), Create(light), FadeIn(tags[0]), FadeIn(tags[1]), run_time=1.2)
        self.play(FadeIn(band), Create(age), FadeIn(tags[3]), FadeIn(tags[4]), run_time=0.7)
        timing.hold_to_read(self, cap2, tags, settle=0.6)

        cross_y = float(np.interp(np.log10(d["crossover_au"]), np.log10(radii),
                                  np.log10(d["disk_period_yr"])))
        cross = Dot(ax.p(d["crossover_au"], 10 ** cross_y), radius=0.08, color=P.FG)
        mark = Dot(ax.p(a_test, d["period_yr"]), radius=0.08, color=P.GREEN)
        mark_lab = layout.label(f"{d['period_myr']:.0f} Myr at {a_test:.0f} AU", font_size=15,
                                color=P.GREEN).next_to(mark, DOWN + RIGHT, buff=0.08)
        cross_lab = layout.label(f"{d['crossover_au']:.0f} AU", font_size=15, color=P.FG)
        cross_lab.next_to(cross, DOWN + RIGHT, buff=0.06)
        cap3 = layout.caption(
            f"Beyond {d['crossover_au']:.0f} AU the {disk['mass_earth']:.0f} M⊕ disc "
            "outpaces the giant planets", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), FadeIn(cross), FadeIn(cross_lab), FadeIn(mark),
                  FadeIn(mark_lab), run_time=0.8)
        timing.hold_to_read(self, cap3, mark_lab, settle=1.0)
        self.play(FadeOut(VGroup(ax, heavy, light, planets, age, band, tags, cross, cross_lab,
                                 mark, mark_lab, cap3)), run_time=0.7)

        # 3. computed: the secular phase plane at 250 AU
        origin = np.array([-2.9, 0.15, 0.0])
        unit = 2.75
        rim = Circle(radius=unit, color=P.MUTED, stroke_width=1.5).move_to(origin)
        e_lo = d["etno_e_min"]
        etno_zone = Annulus(inner_radius=e_lo * unit, outer_radius=unit, color=P.GREEN,
                            stroke_width=0).set_fill(P.GREEN, opacity=0.10).move_to(origin)
        axes_lines = VGroup(
            Line(origin + LEFT * unit, origin + RIGHT * unit, color=P.MUTED, stroke_width=1),
            Line(origin + DOWN * unit, origin + UP * unit, color=P.MUTED, stroke_width=1),
        ).set_stroke(opacity=0.5)
        ends = VGroup(
            layout.label("← disc apse", font_size=14, color=P.MUTED)
            .next_to(origin + RIGHT * unit, RIGHT, buff=0.08),
            layout.label("e = 1", font_size=14, color=P.MUTED)
            .next_to(origin + UP * unit, RIGHT, buff=0.1).shift(DOWN * 0.12),
        )
        trajs = VGroup()
        for t in d["trajectories"]:
            pts = [origin + unit * np.array([k, h, 0.0]) for k, h in zip(t["k"], t["h"])]
            inside = [p for p, k, h in zip(pts, t["k"], t["h"]) if k * k + h * h <= 1.0]
            col = P.ORANGE if t["librates"] else P.TEAL
            trajs.add(polyline(pts if len(inside) == len(pts) else inside, col,
                               stroke_width=2.4 if t["librates"] else 1.8,
                               opacity=0.95 if t["librates"] else 0.7))
        ang = np.radians(d["forced_dvarpi_deg"])
        forced = Dot(origin + unit * d["e_forced"] * np.array([np.cos(ang), np.sin(ang), 0]),
                     radius=0.07, color=P.FG).set_z_index(4)
        key = VGroup(
            layout.label("phase plane at 250 AU", font_size=19, weight="BOLD"),
            layout.label("radius = eccentricity", font_size=15, color=P.MUTED),
            layout.label("angle = apse relative to the disc's", font_size=15, color=P.MUTED),
            layout.label(f"held by the disc:  e < {d['e_libration_max']:.2f}", font_size=17,
                         color=P.ORANGE),
            layout.label("circulating: not held", font_size=17, color=P.TEAL),
            layout.label(f"observed distant orbits:  e ≥ {e_lo:.2f}", font_size=17,
                         color=P.GREEN),
            layout.label(f"held orbits point {abs(d['forced_dvarpi_deg']):.0f}° from the disc apse",
                         font_size=17),
            layout.label(f"paper: {d['paper_forced_dvarpi_deg']:.0f}°, anti-aligned",
                         font_size=17, color=P.MUTED),
        ).arrange(DOWN, buff=0.19, aligned_edge=LEFT).move_to([3.75, 0.2, 0])
        cap4 = layout.caption("Computed: which orbits the disc's gravity holds in place",
                              font_size=22)
        self.play(FadeIn(rim), FadeIn(axes_lines), FadeIn(ends), FadeIn(key[:3]), FadeIn(cap4),
                  run_time=0.9)
        self.play(LaggedStart(*[Create(t) for t in trajs], lag_ratio=0.12), FadeIn(forced),
                  FadeIn(key[3:5]), run_time=2.2)
        timing.hold_to_read(self, cap4, key[:5], settle=0.6)
        riders = VGroup(*[Dot(t.get_start(), radius=0.06, color=t.get_color()).set_z_index(5)
                          for t in trajs])
        cap4b = layout.caption("Follow each orbit: small loops keep their apse near the disc's apse",
                               font_size=22)
        self.play(FadeIn(riders), FadeOut(cap4), FadeIn(cap4b), run_time=0.5)
        self.play(*[MoveAlongPath(r, t) for r, t in zip(riders, trajs)], run_time=4.0,
                  rate_func=linear)
        timing.hold_to_read(self, cap4b, settle=0.3)
        self.play(FadeOut(riders), run_time=0.3)
        cap4 = cap4b
        cap5 = layout.caption(
            f"This linear model holds {d['n_etno_librating']} of "
            f"{len(d['etno_eccentricities'])} observed orbits; the paper's full model holds them",
            font_size=22)
        self.play(FadeIn(etno_zone), FadeIn(key[5:]), FadeOut(cap4), FadeIn(cap5), run_time=0.9)
        timing.hold_to_read(self, cap5, key[5:], settle=1.2)
        self.play(FadeOut(cap5), run_time=0.4)

        layout.show_takeaway(
            self, f"A {disk['mass_earth']:.0f} M⊕ disc can steer distant orbits, if that much mass is there.")
