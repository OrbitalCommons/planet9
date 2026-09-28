"""Sheppard & Trujillo (2016) -- new extreme trans-Neptunian objects.

The first new objects found after the Planet Nine proposal, by a survey built
to look for them. Two are detached and distant enough to join the clustered
sample: 2014 SR349 lines up with it, 2013 FT28 points the opposite way. Both
share the argument-of-perihelion clustering. A third, 2014 FE72, reaches the
outer Oort cloud. The orbits and the circular statistics are the crate's own
(anim.json -> papers -> p9-2016-sheppard-etnos).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Rectangle,
    ReplacementTransform,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2016-sheppard-etnos"


class TopView:
    """The ecliptic plane seen from the north: AU -> scene units, turned by
    ``rotate_deg``."""

    def __init__(self, au_per_unit, sun_at, rotate_deg=0.0):
        self.k = 1.0 / au_per_unit
        self.sun_at = np.array(sun_at, dtype=float)
        t = np.deg2rad(rotate_deg)
        self.rot = np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

    def p(self, x, y):
        v = self.rot @ np.array([x, y])
        return self.sun_at + self.k * np.array([v[0], v[1], 0.0])

    def orbit(self, obj, color, width=2.2, opacity=1.0):
        m = VMobject(color=color, stroke_width=width)
        m.set_points_as_corners([self.p(x, y) for x, y, _ in obj["track"]])
        return m.set_stroke(opacity=opacity)

    def direction(self, lon_deg, length_au):
        t = np.deg2rad(lon_deg)
        return self.p(length_au * np.cos(t), length_au * np.sin(t))

    def aphelion(self, obj):
        far = max(obj["track"], key=lambda q: q[0] ** 2 + q[1] ** 2)
        return self.p(far[0], far[1])


def plate(label, buff=0.04):
    box = Rectangle(width=label.width + 2 * buff, height=label.height + 2 * buff,
                    stroke_width=0).move_to(label)
    return box.set_fill(P.BG, opacity=0.85)


def tagged(text, at, font_size=15, color=P.GREEN, weight="NORMAL"):
    lab = layout.label(text, font_size=font_size, color=color, weight=weight).move_to(at)
    return VGroup(plate(lab), lab)


def scale_bar(view, au, at):
    a = np.array(at, dtype=float)
    b = a + np.array([au * view.k, 0, 0])
    bar = VGroup(Line(a, b, color=P.MUTED, stroke_width=2),
                 Line(a + UP * 0.06, a + DOWN * 0.06, color=P.MUTED, stroke_width=2),
                 Line(b + UP * 0.06, b + DOWN * 0.06, color=P.MUTED, stroke_width=2))
    lab = layout.label(f"{au:,.0f} AU", font_size=13, color=P.MUTED).next_to(bar, DOWN, buff=0.08)
    return VGroup(bar, lab)


def ring(centre, radius, title, subtitle):
    """A compass rose for an orbital angle: 0 deg to the right, increasing
    anticlockwise, with a plain-words note on what the angle is measured from."""
    c = np.array(centre, dtype=float)
    g = VGroup(Circle(radius=radius, color=P.MUTED, stroke_width=1.5).move_to(c))
    for deg in (0, 90, 180, 270):
        u = np.array([np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg)), 0.0])
        g.add(Line(c + u * (radius - 0.08), c + u * (radius + 0.08), color=P.MUTED,
                   stroke_width=1.5))
        g.add(layout.label(f"{deg}°", font_size=12, color=P.MUTED)
              .move_to(c + u * (radius + 0.32)))
    g.add(layout.label(title, font_size=17, color=P.FG).move_to(c + UP * (radius + 1.0)))
    g.add(layout.label(subtitle, font_size=13, color=P.MUTED).move_to(c + UP * (radius + 0.68)))
    return g


def spoke(centre, radius, deg, color, width=2.5, opacity=1.0):
    c = np.array(centre, dtype=float)
    u = np.array([np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg)), 0.0])
    line = Line(c, c + u * radius, color=color, stroke_width=width).set_stroke(opacity=opacity)
    dot = Dot(c + u * radius, radius=0.06, color=color).set_opacity(opacity)
    return VGroup(line, dot)


class SheppardEtnos2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        known = d["known"]
        new = d["new"]
        members = [o for o in new if o["in_sample"]]
        sr349 = next(o for o in new if o["name"] == "2014 SR349")
        ft28 = next(o for o in new if o["anti_aligned"])
        fe72 = next(o for o in new if o["outer_oort"])
        turn = -d["varpi_before"]["mean_deg"]

        self.add(paper.scene_header(CRATE))

        # 1. the known cluster, and the two new members
        near = TopView(au_per_unit=200.0, sun_at=(-0.6, 0.3, 0.0), rotate_deg=turn)
        sun = Dot(near.sun_at, radius=0.06, color=P.SUN).set_z_index(5)
        old = VGroup(*[near.orbit(o, P.GREEN, width=1.6, opacity=0.45) for o in known])
        fresh = VGroup(*[near.orbit(o, P.GREEN, width=3.2) for o in members])
        bar = scale_bar(near, 500.0, (-6.3, -2.35, 0))
        old_lab = layout.label(f"the {d['n_before']} orbits of early 2016", font_size=15,
                               color=P.GREEN).set_opacity(0.7)
        old_lab.move_to([-5.0, 2.8, 0])
        cap = layout.caption("Early 2016: six distant orbits, perihelia all on one side",
                             font_size=22)
        self.play(FadeIn(sun), FadeIn(bar))
        self.play(LaggedStart(*[Create(o) for o in old], lag_ratio=0.15), FadeIn(old_lab),
                  FadeIn(cap), run_time=2.0)
        timing.hold_to_read(self, cap, settle=0.5)

        header = ["new object", "a (AU)", "q (AU)", "tilt"]
        table = VGroup()
        for r, row in enumerate([header] + [
                [o["name"], f"{o['a']:,.0f}", f"{o['q']:.0f}", f"{o['i_deg']:.0f}°"]
                for o in new]):
            joins = r > 0 and new[r - 1]["in_sample"]
            colour = P.MUTED if r == 0 else (P.GREEN if joins else P.FG)
            for c, text in enumerate(row):
                cell = layout.label(text, font_size=15, color=colour,
                                    weight="BOLD" if joins else "NORMAL")
                x = 3.0 + [0.0, 2.2, 3.05, 3.75][c]
                cell.move_to([x, 2.6 - 0.36 * r, 0], aligned_edge=LEFT if c == 0 else RIGHT)
                table.add(cell)
        note = layout.label(
            f"bold: a > {d['a_threshold']:.0f} AU and q > {d['q_threshold']:.0f} AU",
            font_size=13, color=P.MUTED)
        note.next_to(table, DOWN, buff=0.18, aligned_edge=LEFT)
        names = VGroup(
            tagged("2014 SR349", near.aphelion(sr349) + UP * 0.25 + LEFT * 0.3, weight="BOLD"),
            tagged("2013 FT28", near.aphelion(ft28) + DOWN * 0.3 + RIGHT * 0.2, weight="BOLD"),
        )
        cap2 = layout.caption(
            f"2014 SR349 joins the cluster; 2013 FT28 points "
            f"{d['ft28_offset_deg']:.0f}° away from it", font_size=22)
        self.play(FadeIn(table, lag_ratio=0.02), FadeIn(note), FadeOut(cap))
        self.play(LaggedStart(*[Create(o) for o in fresh], lag_ratio=0.5), FadeIn(cap2),
                  run_time=2.4)
        self.play(FadeIn(names))
        timing.hold_to_read(self, cap2, settle=1.4)

        # 2. zoom out: the outer Oort cloud object
        far = TopView(au_per_unit=620.0, sun_at=(3.2, 1.9, 0.0), rotate_deg=turn)
        old_far = VGroup(*[far.orbit(o, P.GREEN, width=1.2, opacity=0.45) for o in known])
        fresh_far = VGroup(*[far.orbit(o, P.GREEN, width=2.0) for o in members])
        sun_far = Dot(far.sun_at, radius=0.05, color=P.SUN).set_z_index(5)
        bar_far = scale_bar(far, 2000.0, (2.6, -2.35, 0))
        big = far.orbit(fe72, P.GREEN, width=3.2)
        aph = max(np.hypot(x, y) for x, y, _ in fe72["track"])
        big_lab = VGroup(
            layout.label("2014 FE72", font_size=17, color=P.GREEN, weight="BOLD"),
            layout.label(f"a = {fe72['a']:,.0f} AU,  perihelion {fe72['q']:.0f} AU",
                         font_size=15, color=P.GREEN),
            layout.label(f"reaches {aph:,.0f} AU from the Sun", font_size=15, color=P.GREEN),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        big_lab.move_to([-4.3, 1.6, 0])
        cap3 = layout.caption(
            "2014 FE72: the first outer Oort cloud object with a perihelion beyond Neptune",
            font_size=22)
        self.play(FadeOut(VGroup(table, note, names, old_lab, cap2)), run_time=0.6)
        self.play(ReplacementTransform(old, old_far), ReplacementTransform(fresh, fresh_far),
                  ReplacementTransform(sun, sun_far), ReplacementTransform(bar, bar_far),
                  run_time=1.6)
        self.play(Create(big), FadeIn(big_lab), FadeIn(cap3), run_time=2.2)
        timing.hold_to_read(self, cap3, big_lab, settle=1.0)
        self.play(FadeOut(VGroup(old_far, fresh_far, sun_far, bar_far, big, big_lab, cap3)))

        # 3. the two angles: argument and longitude of perihelion
        left, right, rad = (-3.4, -0.35, 0.0), (3.4, -0.35, 0.0), 1.6
        rose_w = ring(left, rad, "argument of perihelion ω",
                      "counted from where the orbit crosses the ecliptic")
        rose_v = ring(right, rad, "longitude of perihelion ϖ",
                      "counted from a fixed direction in space")
        old_w = VGroup(*[spoke(left, rad, o["omega_deg"], P.GREEN, 1.8, 0.45) for o in known])
        old_v = VGroup(*[spoke(right, rad, o["varpi_deg"], P.GREEN, 1.8, 0.45) for o in known])
        new_w = VGroup(*[spoke(left, rad, o["omega_deg"], P.GREEN, 3.5) for o in members])
        new_v = VGroup(*[spoke(right, rad, o["varpi_deg"], P.GREEN, 3.5) for o in members])
        tags = VGroup()
        for o in members:
            for centre, key in ((left, "omega_deg"), (right, "varpi_deg")):
                u = np.array([np.cos(np.deg2rad(o[key])), np.sin(np.deg2rad(o[key])), 0.0])
                side = RIGHT if u[0] >= 0 else LEFT
                lab = layout.label(o["name"], font_size=13, color=P.GREEN, weight="BOLD")
                lab.next_to(np.array(centre) + u * (rad + 0.12), side, buff=0.08)
                if abs(u[0]) < 0.4:
                    lab.shift(side * 0.4)
                else:
                    lab.shift(UP * 0.22 * np.sign(u[1]))
                tags.add(VGroup(plate(lab), lab))
        w = d["omega_after"]
        stat_w = layout.label(
            f"all {w['n']} clustered around {w['mean_deg']:.0f}°  (spread {w['std_deg']:.0f}°)",
            font_size=15, color=P.FG).move_to([left[0], -2.8, 0])
        stat_v = layout.label(
            f"2013 FT28 sits {d['ft28_offset_deg']:.0f}° from the rest  (paper: about 180°)",
            font_size=15, color=P.FG).move_to([right[0], -2.8, 0])
        cap4 = layout.caption(
            "Both new orbits share the clustering in ω, but one is flipped in ϖ",
            font_size=22)
        self.play(FadeIn(rose_w), FadeIn(rose_v))
        self.play(FadeIn(old_w, lag_ratio=0.1), FadeIn(old_v, lag_ratio=0.1), run_time=1.0)
        self.play(FadeIn(new_w), FadeIn(new_v), FadeIn(tags), FadeIn(cap4), run_time=1.2)
        self.play(FadeIn(stat_w), FadeIn(stat_v))
        timing.hold_to_read(self, cap4, stat_w, stat_v, settle=1.2)
        self.play(FadeOut(cap4), FadeOut(stat_w), FadeOut(stat_v))

        layout.show_takeaway(
            self, f"The sample grows from {d['n_before']} to {d['n_after']}, "
                  "and gains its first anti-aligned orbit.")
