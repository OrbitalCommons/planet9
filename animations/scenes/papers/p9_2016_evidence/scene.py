"""Batygin & Brown (2016) -- "Evidence for a Distant Giant Planet in the Solar System".

The founding paper. Six dynamically stable distant Kuiper belt objects have
perihelia that point the same way and orbital planes that tilt the same way;
the chance of both at once is a few in 100,000; a ten-Earth-mass planet on an
eccentric orbit on the *opposite* side of the Sun would hold them there.
Everything drawn here is the crate's own (anim.json -> papers ->
p9-2016-evidence): the six orbits to scale, the Monte Carlo null and the
nominal perturber.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arrow,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Polygon,
    Rectangle,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2016-evidence"


class TopView:
    """The ecliptic plane seen from the north: AU -> scene units, with the
    view turned by ``rotate_deg`` so the long axis of the picture runs across
    the frame."""

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

    def circle(self, radius_au, color=P.MUTED, width=1.2):
        return Circle(radius=radius_au * self.k, color=color,
                      stroke_width=width).move_to(self.sun_at)

    def direction(self, lon_deg, length_au):
        t = np.deg2rad(lon_deg)
        return self.p(length_au * np.cos(t), length_au * np.sin(t))

    def aphelion(self, obj):
        far = max(obj["track"], key=lambda q: q[0] ** 2 + q[1] ** 2)
        return self.p(far[0], far[1])


def scale_bar(view, au, at):
    a = np.array(at, dtype=float)
    b = a + np.array([au * view.k, 0, 0])
    bar = VGroup(Line(a, b, color=P.MUTED, stroke_width=2),
                 Line(a + UP * 0.06, a + DOWN * 0.06, color=P.MUTED, stroke_width=2),
                 Line(b + UP * 0.06, b + DOWN * 0.06, color=P.MUTED, stroke_width=2))
    lab = layout.label(f"{au:.0f} AU", font_size=13, color=P.MUTED).next_to(bar, DOWN, buff=0.08)
    return VGroup(bar, lab)


def plate(label, buff=0.04):
    """A background plate so a label stays readable where orbits cross it."""
    box = Rectangle(width=label.width + 2 * buff, height=label.height + 2 * buff,
                    stroke_width=0).move_to(label)
    return box.set_fill(P.BG, opacity=0.85)


def north_arrow(view, at):
    """Which way ecliptic longitude 0 points in the turned view."""
    a = np.array(at, dtype=float)
    tip = a + (view.direction(0.0, 1.0) - view.sun_at) / view.k * 0.55
    arrow = Arrow(a, tip, buff=0, color=P.MUTED, stroke_width=2,
                  max_tip_length_to_length_ratio=0.25)
    lab = layout.label("ecliptic longitude 0°", font_size=13, color=P.MUTED)
    lab.next_to(arrow, RIGHT, buff=0.12)
    return VGroup(arrow, lab)


class Evidence2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        p9 = d["p9"]

        self.add(paper.scene_header(CRATE))

        # One picture for the whole scene, turned so the shared perihelion
        # direction points right: the six orbits reach out to the left, and
        # the planet's orbit (beat 3) to the right.
        view = TopView(au_per_unit=210.0, sun_at=(-0.7, 0.3, 0.0),
                       rotate_deg=-d["mean_varpi_deg"])
        sun = Dot(view.sun_at, radius=0.06, color=P.SUN).set_z_index(5)
        neptune = view.circle(30.0)
        orbits_g = VGroup(*[view.orbit(o, P.GREEN, width=2.0, opacity=0.9) for o in objs])
        world = VGroup(neptune, orbits_g, sun)
        compass = north_arrow(view, (-6.5, 2.95, 0))

        # 1. six real orbits, to scale
        names = VGroup()
        for o in objs:
            far = view.aphelion(o)
            out = far - view.sun_at
            out /= np.linalg.norm(out)
            lab = layout.label(o["name"], font_size=15, color=P.GREEN)
            lab.move_to(far + out * 0.22 + LEFT * 0.5 * lab.width * abs(out[0]))
            names.add(VGroup(plate(lab), lab))
        header = ["", "a (AU)", "q (AU)", "tilt"]
        cells = [[o["name"], f"{o['a']:.0f}", f"{o['q']:.0f}", f"{o['i_deg']:.0f}°"] for o in objs]
        table = VGroup()
        for r, row in enumerate([header] + cells):
            for c, text in enumerate(row):
                cell = layout.label(text, font_size=15, color=P.MUTED if r == 0 else P.FG)
                x = 2.7 + [0.0, 2.05, 2.95, 3.65][c]
                cell.move_to([x, 1.9 - 0.36 * r, 0], aligned_edge=LEFT if c == 0 else RIGHT)
                table.add(cell)
        bar = scale_bar(view, 500.0, (3.4, -2.3, 0))
        cap = layout.caption(
            f"The {d['n_sample']} most distant stable orbits known in 2016, seen from above",
            font_size=22)
        self.play(FadeIn(sun), Create(neptune), FadeIn(bar), FadeIn(compass))
        self.play(LaggedStart(*[Create(o) for o in orbits_g], lag_ratio=0.25),
                  FadeIn(cap), run_time=2.6)
        self.play(FadeIn(names, lag_ratio=0.1), FadeIn(table, lag_ratio=0.02), run_time=1.0)
        timing.hold_to_read(self, cap, settle=1.6)

        peri = VGroup(*[
            Arrow(view.sun_at, view.direction(o["varpi_deg"], 330.0), buff=0, color=P.ORANGE,
                  stroke_width=2.5, max_tip_length_to_length_ratio=0.1)
            for o in objs])
        peri_note = layout.label("perihelion directions", font_size=15, color=P.ORANGE)
        peri_note.next_to(peri, DOWN, buff=0.15).shift(RIGHT * 0.6)
        spread = max(o["varpi_deg"] for o in objs) - min(o["varpi_deg"] for o in objs)
        cap2 = layout.caption(
            f"All six perihelia fall within {spread:.0f}° of each other; "
            f"all six planes tilt the same way", font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(LaggedStart(*[Create(a) for a in peri], lag_ratio=0.12), FadeIn(peri_note),
                  FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, settle=1.0)
        self.play(FadeOut(VGroup(world, names, table, peri, peri_note, bar, compass, cap2)))

        # 2. how often chance does this: the crate's Monte Carlo null
        null = np.array(d["null_trials"])
        rv, rp = d["r_bar_varpi"], d["r_pole"]
        y_lo = 0.93
        ax, labels = widgets.labeled_axes(
            [0, 1, 0.2], [y_lo, 1.0, 0.01], x_label="alignment of the perihelia",
            y_label="alignment of the orbital planes", y_rotate=True, numbers=True,
            x_length=6.4, y_length=4.3, shift_down=-0.25)
        plot = VGroup(ax, labels).shift(LEFT * 2.9)
        dots = VGroup(*[Dot(ax.c2p(x, max(y, y_lo)), radius=0.022, color=P.MUTED).set_opacity(0.75)
                        for x, y in null])
        a, b = ax.c2p(rv, rp), ax.c2p(1.0, 1.0)
        corner = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=P.GREEN, stroke_width=1.5)
        corner.set_fill(P.GREEN, opacity=0.15)
        obs = Dot(a, radius=0.08, color=P.GREEN).set_z_index(4)
        obs_lab = layout.label("the six real orbits", font_size=14, color=P.GREEN)
        obs_lab.next_to(obs, DOWN, buff=0.12)
        cap3 = layout.caption(
            f"{len(null)} random sets of six orbits: none is this aligned in both",
            font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(dots, lag_ratio=0.002), FadeIn(cap3), run_time=1.6)
        self.play(FadeIn(corner), FadeIn(obs), FadeIn(obs_lab))

        rows = VGroup(
            layout.label(f"perihelia this aligned:  {100 * d['p_varpi']:.2f}%", font_size=19),
            layout.label(f"planes this aligned:  {100 * d['p_pole']:.2f}%", font_size=19),
            layout.label(f"both at once:  {100 * d['p_joint']:.4f}%", font_size=22,
                         color=P.GREEN, weight="BOLD"),
            layout.label("paper:  0.007%", font_size=17, color=P.MUTED),
            layout.label(f"({d['n_trials'] / 1e6:.0f} million random trials)", font_size=13,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        rows.move_to([4.2, 0.45, 0])
        self.play(FadeIn(rows, lag_ratio=0.2), run_time=1.4)
        timing.hold_to_read(self, cap3, rows, settle=1.2)
        self.play(FadeOut(VGroup(plot, dots, corner, obs, obs_lab, rows, cap3)))

        # 3. the perturber that would do it
        p9_orbit = view.orbit(p9, P.BLUE, width=3.0)
        peri_dot = Dot(view.p(*p9["peri_xyz"][:2]), radius=0.05, color=P.BLUE)
        p9_lab = layout.label("Planet Nine", font_size=18, color=P.BLUE, weight="BOLD")
        elems = VGroup(
            layout.label(f"{d['mass']:.0f} Earth masses", font_size=15, color=P.BLUE),
            layout.label(f"a = {d['a']:.0f} AU    e = {d['e']:.1f}", font_size=15, color=P.BLUE),
            layout.label(f"perihelion {p9['q']:.0f} AU", font_size=15, color=P.BLUE),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        tag = VGroup(p9_lab, elems).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        tag.move_to([2.5, 0.3, 0])
        kbo_lab = layout.label("the six orbits", font_size=16, color=P.GREEN)
        kbo_lab.move_to([-5.2, 1.0, 0])
        cap4 = layout.caption(
            f"The proposed planet: its perihelion {d['p9_offset_deg']:.0f}° away "
            f"from theirs, on the far side of the Sun", font_size=22)
        self.play(FadeIn(world), FadeIn(bar), FadeIn(compass), FadeIn(kbo_lab))
        self.play(Create(p9_orbit), FadeIn(cap4), run_time=2.2)
        self.play(FadeIn(peri_dot), FadeIn(tag))
        timing.hold_to_read(self, cap4, elems, settle=1.2)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, "Six aligned orbits, odds of a few in 100,000, and a planet that would explain them.")
