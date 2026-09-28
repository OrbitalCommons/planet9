"""Batygin et al. (2019) -- the Planet Nine hypothesis (review).

Three years after the proposal, the review revises the planet: about half the
mass, closer in, on a rounder orbit. A smaller semi-major axis also moves the
edge of the clustered region, and a closer aphelion makes the planet brighter
where it spends most of its time. The orbits, the critical semi-major axes
and the brightness ranges are the crate's own (anim.json -> papers ->
p9-2019-review).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Circle,
    Create,
    DashedLine,
    DashedVMobject,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Rectangle,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2019-review"

# Published in the paper; drawn only as a labelled comparison.
PAPER_A_C = 250.0
# Distant-object orbits drawn in the top view must fit the frame.
MAX_APHELION_AU = 1350.0


class TopView:
    """The ecliptic plane seen from the north, turned by ``rotate_deg``."""

    def __init__(self, au_per_unit, sun_at, rotate_deg=0.0):
        self.k = 1.0 / au_per_unit
        self.sun_at = np.array(sun_at, dtype=float)
        t = np.deg2rad(rotate_deg)
        self.rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

    def p(self, x, y):
        v = self.rot @ np.array([x, y])
        return self.sun_at + self.k * np.array([v[0], v[1], 0.0])

    def orbit(self, obj, color, width=2.2, opacity=1.0, dashed=False):
        m = VMobject(color=color, stroke_width=width)
        m.set_points_as_corners([self.p(x, y) for x, y, _ in obj["track"]])
        m.set_stroke(opacity=opacity)
        return DashedVMobject(m, num_dashes=70) if dashed else m


def scale_bar(view, au, at):
    a = np.array(at, dtype=float)
    b = a + np.array([au * view.k, 0, 0])
    bar = VGroup(Line(a, b, color=P.MUTED, stroke_width=2),
                 Line(a + UP * 0.06, a + DOWN * 0.06, color=P.MUTED, stroke_width=2),
                 Line(b + UP * 0.06, b + DOWN * 0.06, color=P.MUTED, stroke_width=2))
    lab = layout.label(f"{au:.0f} AU", font_size=13, color=P.MUTED).next_to(bar, DOWN, buff=0.08)
    return VGroup(bar, lab)


def elements_text(o):
    return (f"{o['mass_earth']:.0f} Earth masses,  a = {o['a']:.0f} AU,  e = {o['e']:.2f},  "
            f"tilt {o['i_deg']:.0f}°")


class Review2019(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        old = d["orbit_2016"]
        new = d["orbit_2019"]
        fits = d["fits"]
        objs = d["objects"]

        self.add(paper.scene_header(CRATE))

        # 1. the revised planet, to scale
        view = TopView(au_per_unit=300.0, sun_at=(-2.2, -0.1, 0.0),
                       rotate_deg=180.0 - old["varpi_deg"])
        sun = Dot(view.sun_at, radius=0.06, color=P.SUN).set_z_index(5)
        neptune = Circle(radius=30.0 * view.k, color=P.MUTED, stroke_width=1.2)
        neptune.move_to(view.sun_at)
        # The most elongated orbits reach past the frame edge at this scale.
        fitting = [o for o in objs if o["a"] * (1 + o["e"]) <= MAX_APHELION_AU]
        kbos = VGroup(*[view.orbit(o, P.GREEN, width=1.2, opacity=0.3)
                        for o in fitting])
        orbit_old = view.orbit(old, P.MUTED, width=2.6, dashed=True)
        family = VGroup(*[view.orbit(f, P.BLUE, width=1.4, opacity=0.45) for f in fits])
        orbit_new = view.orbit(new, P.BLUE, width=3.6)
        bar = scale_bar(view, 500.0, (3.0, -2.55, 0))
        rows_old = VGroup(
            layout.label("2016", font_size=17, color=P.MUTED, weight="BOLD"),
            layout.label(elements_text(old), font_size=15, color=P.MUTED),
            layout.label(f"perihelion {old['q']:.0f} AU,  aphelion "
                         f"{old['a'] * (1 + old['e']):,.0f} AU", font_size=15, color=P.MUTED),
        ).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        rows_new = VGroup(
            layout.label("2019", font_size=17, color=P.BLUE, weight="BOLD"),
            layout.label(elements_text(new), font_size=15, color=P.BLUE),
            layout.label(f"perihelion {new['q']:.0f} AU,  aphelion "
                         f"{new['a'] * (1 + new['e']):,.0f} AU", font_size=15, color=P.BLUE),
            layout.label(f"thin lines: the other best fits, a = "
                         f"{min(f['a'] for f in fits):.0f}-{max(f['a'] for f in fits):.0f} AU",
                         font_size=13, color=P.BLUE).set_opacity(0.8),
        ).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        rows = VGroup(rows_old, rows_new).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
        rows.move_to([2.15, 1.3, 0], aligned_edge=LEFT)
        kbo_lab = layout.label(
            f"the distant objects ({len(objs) - len(fitting)} too elongated to fit)",
            font_size=14, color=P.GREEN)
        kbo_lab.set_opacity(0.8).move_to([-4.6, 2.6, 0])
        cap = layout.caption(f"The 2016 planet: {old['mass_earth']:.0f} Earth masses on a long, "
                             f"eccentric orbit",
                             font_size=22)
        self.play(FadeIn(sun), Create(neptune), FadeIn(bar), FadeIn(kbos), FadeIn(kbo_lab))
        self.play(Create(orbit_old), FadeIn(rows_old), FadeIn(cap), run_time=1.8)
        timing.hold_to_read(self, cap, rows_old, settle=0.6)
        cap2 = layout.caption(
            "The 2019 revision: half the mass, closer in, on a rounder orbit", font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(Create(family), Create(orbit_new), FadeIn(rows_new), FadeIn(cap2),
                  run_time=2.2)
        timing.hold_to_read(self, cap2, rows_new, settle=1.2)
        self.play(FadeOut(VGroup(sun, neptune, bar, kbos, kbo_lab, orbit_old, family,
                                 orbit_new, rows, cap2)))

        # 2. where the planet's grip starts
        x0, x1, a_max = -6.0, 6.2, 800.0

        def ax_x(a):
            return x0 + (x1 - x0) * a / a_max

        y = 0.2
        axis = Line([x0, y, 0], [x1, y, 0], color=P.MUTED, stroke_width=2)
        ticks = VGroup()
        for a in range(0, int(a_max) + 1, 100):
            ticks.add(Line([ax_x(a), y - 0.08, 0], [ax_x(a), y + 0.08, 0], color=P.MUTED,
                           stroke_width=1.5))
            ticks.add(layout.label(f"{a}", font_size=13, color=P.MUTED)
                      .move_to([ax_x(a), y - 0.32, 0]))
        axis_lab = layout.label("semi-major axis of a distant object (AU)", font_size=16,
                                color=P.FG).move_to([0.1, y - 0.75, 0])
        shown = [o for o in objs if o["a"] <= a_max]
        off = [o for o in objs if o["a"] > a_max]
        marks = VGroup(*[Line([ax_x(o["a"]), y + 0.12, 0], [ax_x(o["a"]), y + 0.62, 0],
                              color=P.GREEN, stroke_width=3.5) for o in shown])
        marks_lab = layout.label(
            f"the {len(objs)} distant objects" + (f"  ({len(off)} beyond {a_max:.0f} AU)"
                                                  if off else ""),
            font_size=15, color=P.GREEN)
        marks_lab.move_to([ax_x(500), y + 0.95, 0])

        def edge(a_c, color, text, level, side):
            ln = DashedLine([ax_x(a_c), y - 0.05, 0], [ax_x(a_c), y - 1.35 - 0.45 * level, 0],
                            color=color, stroke_width=2.5)
            lab = layout.label(text, font_size=15, color=color)
            lab.next_to(ln.get_end(), side, buff=0.1)
            return VGroup(ln, lab)

        # Two edges 25 AU apart: the earlier one labelled to its left, the new
        # one to its right and lower, so neither line crosses the other's label.
        edge_old = edge(d["a_c_2016"], P.MUTED, f"2016 planet:  a_c = {d['a_c_2016']:.0f} AU",
                        0, LEFT)
        edge_new = edge(d["a_c"], P.BLUE,
                        f"2019 planet:  a_c = {d['a_c']:.0f} AU   (paper: about "
                        f"{PAPER_A_C:.0f} AU)", 1, RIGHT)
        zone = Rectangle(width=ax_x(a_max) - ax_x(d["a_c"]), height=1.2, stroke_width=0)
        zone.set_fill(P.BLUE, opacity=0.10).move_to([(ax_x(d["a_c"]) + ax_x(a_max)) / 2,
                                                    y + 0.55, 0])
        zone_lab = layout.label("held anti-aligned by the planet", font_size=15, color=P.BLUE)
        zone_lab.move_to([ax_x(560), y + 1.45, 0])
        free_lab = layout.label("free to precess", font_size=15, color=P.MUTED)
        free_lab.move_to([ax_x(110), y + 1.45, 0])
        cap3 = layout.caption(
            "Beyond a critical distance a_c the planet locks orbits in place; inside it, not",
            font_size=22)
        self.play(Create(axis), FadeIn(ticks), FadeIn(axis_lab))
        self.play(LaggedStart(*[Create(m) for m in marks], lag_ratio=0.08), FadeIn(marks_lab),
                  run_time=1.4)
        self.play(FadeIn(zone), FadeIn(zone_lab), FadeIn(free_lab), Create(edge_new),
                  FadeIn(cap3), run_time=1.4)
        self.play(Create(edge_old))
        timing.hold_to_read(self, cap3, edge_new, settle=1.0)
        n_in = sum(1 for o in objs if o["a"] >= d["a_c"])
        cap4 = layout.caption(
            f"{n_in} of the {len(objs)} clustered objects lie beyond the revised edge",
            font_size=22)
        self.play(FadeOut(cap3), run_time=0.4)
        self.play(FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=1.0)
        self.play(FadeOut(VGroup(axis, ticks, axis_lab, marks, marks_lab, edge_old, edge_new,
                                 zone, zone_lab, free_lab, cap4)))

        # 3. a closer planet is a brighter one
        v0, v1 = 13.0, 27.0
        bx0, bx1 = -3.4, 6.4

        def bx(v):
            return bx0 + (bx1 - bx0) * (v - v0) / (v1 - v0)

        base = -1.75
        vaxis = Line([bx0, base, 0], [bx1, base, 0], color=P.MUTED, stroke_width=2)
        vticks = VGroup()
        for v in range(int(v0), int(v1) + 1, 2):
            vticks.add(Line([bx(v), base - 0.07, 0], [bx(v), base + 0.07, 0], color=P.MUTED,
                            stroke_width=1.5))
            vticks.add(layout.label(f"{v}", font_size=13, color=P.MUTED)
                       .move_to([bx(v), base - 0.28, 0]))
        vlab = layout.label("apparent magnitude V   (fainter →)", font_size=16, color=P.FG)
        vlab.move_to([bx(20.0), base - 0.62, 0])
        bars, names = VGroup(), VGroup()
        for k, b in enumerate(d["brightness"]):
            yy = 1.85 - 0.62 * k
            lo, hi = b["v_perihelion_bright"], b["v_aphelion_faint"]
            r = Rectangle(width=bx(hi) - bx(lo), height=0.34, stroke_width=0)
            r.set_fill(P.BLUE, opacity=0.75).move_to([(bx(lo) + bx(hi)) / 2, yy, 0])
            bars.add(r)
            text = b["label"].replace(" ME,", " Earth masses,")
            names.add(layout.label(text, font_size=15, color=P.BLUE)
                      .move_to([bx0 - 0.15, yy, 0], aligned_edge=RIGHT))
        surveys = VGroup()
        for s_k, s in enumerate(d["surveys"]):
            ln = DashedLine([bx(s["v_limit"]), base, 0], [bx(s["v_limit"]), 2.3, 0],
                            color=P.PURPLE, stroke_width=2)
            lab = layout.label(s["name"], font_size=13, color=P.PURPLE)
            lab.next_to(ln.get_end(), UP, buff=0.06).shift(UP * 0.28 * (s_k % 2))
            surveys.add(VGroup(ln, lab))
        pluto = next(m for m in dataio.solar_system()["magnitudes"] if m["name"] == "Pluto")
        pl = DashedLine([bx(pluto["app_mag"]), base, 0], [bx(pluto["app_mag"]), 2.3, 0],
                        color=P.MUTED, stroke_width=2)
        pl_lab = layout.label(f"Pluto today, V = {pluto['app_mag']:.1f}", font_size=13,
                              color=P.MUTED).next_to(pl.get_end(), UP, buff=0.06)
        key = layout.label("each bar: from perihelion (bright end) to aphelion (faint end)",
                           font_size=14, color=P.FG)
        key.move_to([bx(20.0), -1.1, 0])
        cap5 = layout.caption(
            f"The revised planets span V ≈ "
            f"{min(b['v_perihelion_bright'] for b in d['brightness']):.0f}-"
            f"{max(b['v_aphelion_faint'] for b in d['brightness']):.0f}: faint, but inside "
            f"the reach of deep surveys",
            font_size=22)
        self.play(Create(vaxis), FadeIn(vticks), FadeIn(vlab), Create(pl), FadeIn(pl_lab))
        self.play(LaggedStart(*[FadeIn(r, shift=RIGHT * 0.2) for r in bars], lag_ratio=0.2),
                  FadeIn(names), FadeIn(key), FadeIn(cap5), run_time=1.6)
        self.play(LaggedStart(*[Create(s) for s in surveys], lag_ratio=0.2), run_time=1.4)
        timing.hold_to_read(self, cap5, key, settle=1.4)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "A smaller, closer planet: brighter, and within reach of the next surveys.")
