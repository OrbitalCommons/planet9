"""Becker et al. (2018) -- discovery and dynamical analysis of 2015 BP519.

The Dark Energy Survey finds a distant object whose orbit is tilted 54 degrees
to the plane of the planets, about twice the tilt of any of the clustered
orbits. The known giant planets cannot tilt an orbit that far; a distant
inclined planet can. The orbits, the secular inclination histories and the
perihelion offsets are the crate's own (anim.json -> papers -> p9-2018-bp519).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arc,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Rectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2018-bp519"


def tilt_line(sun_at, i_deg, reach, color, width, opacity=1.0):
    """An orbital plane seen edge-on along its node: a line through the Sun at
    the orbit's inclination to the plane of the planets."""
    t = np.deg2rad(i_deg)
    end = np.array(sun_at) + reach * np.array([np.cos(t), np.sin(t), 0.0])
    return Line(sun_at, end, color=color, stroke_width=width).set_stroke(opacity=opacity)


def plate(label, buff=0.04):
    box = Rectangle(width=label.width + 2 * buff, height=label.height + 2 * buff,
                    stroke_width=0).move_to(label)
    return box.set_fill(P.BG, opacity=0.85)


def ring(centre, radius, title):
    c = np.array(centre, dtype=float)
    g = VGroup(Circle(radius=radius, color=P.MUTED, stroke_width=1.5).move_to(c))
    for deg in (0, 90, 180, 270):
        u = np.array([np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg)), 0.0])
        g.add(Line(c + u * (radius - 0.08), c + u * (radius + 0.08), color=P.MUTED,
                   stroke_width=1.5))
        g.add(layout.label(f"{deg}°", font_size=12, color=P.MUTED)
              .move_to(c + u * (radius + 0.32)))
    g.add(layout.label(title, font_size=17, color=P.FG).move_to(c + UP * (radius + 0.75)))
    return g


def spoke(centre, radius, deg, color, width=2.5, opacity=1.0):
    c = np.array(centre, dtype=float)
    u = np.array([np.cos(np.deg2rad(deg)), np.sin(np.deg2rad(deg)), 0.0])
    line = Line(c, c + u * radius, color=color, stroke_width=width).set_stroke(opacity=opacity)
    dot = Dot(c + u * radius, radius=0.06, color=color).set_opacity(opacity)
    return VGroup(line, dot)


class Bp519Discovery2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        bp = d["bp519"]
        sample = d["sample"]

        self.add(paper.scene_header(CRATE))

        # 1. tilt: every orbital plane seen edge-on, the plane of the planets flat
        sun_at = np.array([-4.3, -1.95, 0.0])
        reach = 6.0
        plane = DashedLine(sun_at + LEFT * 2.2, sun_at + RIGHT * 7.4, color=P.ORANGE,
                           stroke_width=1.8).set_stroke(opacity=0.85)
        plane_lab = layout.label("plane of the planets", font_size=15, color=P.ORANGE)
        plane_lab.next_to(plane.get_end(), DOWN, buff=0.12, aligned_edge=RIGHT)
        sun = Dot(sun_at, radius=0.07, color=P.SUN).set_z_index(5)
        others = VGroup(*[tilt_line(sun_at, o["i_deg"], reach, P.GREEN, 2.0, 0.55)
                          for o in sample])
        tilts = [o["i_deg"] for o in sample]
        pluto = dataio.body("Pluto")
        ladder = VGroup(
            layout.label("tilt to the plane of the planets", font_size=15, color=P.MUTED),
            layout.label(f"Pluto:  {pluto['i_deg']:.0f}°", font_size=17, color=P.FG),
            layout.label(f"the {len(sample)} clustered orbits:  {min(tilts):.0f}°-"
                         f"{max(tilts):.0f}°", font_size=17, color=P.GREEN),
        ).arrange(DOWN, buff=0.18, aligned_edge=LEFT)
        ladder.move_to([2.4, 2.3, 0], aligned_edge=UP + LEFT)
        cap = layout.caption(
            "Each orbit's plane seen edge-on: the clustered orbits tilt only modestly",
            font_size=22)
        self.play(FadeIn(sun), Create(plane), FadeIn(plane_lab))
        self.play(LaggedStart(*[Create(o) for o in others], lag_ratio=0.1),
                  FadeIn(ladder[:3]), FadeIn(cap), run_time=2.2)
        timing.hold_to_read(self, cap, ladder, settle=0.8)

        new = tilt_line(sun_at, bp["i_deg"], reach, P.GREEN, 4.5)
        arc = Arc(radius=1.3, start_angle=0, angle=np.deg2rad(bp["i_deg"]), color=P.GREEN,
                  stroke_width=2.5, arc_center=sun_at)
        mid = np.deg2rad(bp["i_deg"] / 2)
        arc_lab = layout.label(f"{bp['i_deg']:.0f}°", font_size=18, color=P.GREEN,
                               weight="BOLD").move_to(sun_at + 1.65 * np.array(
                                   [np.cos(mid), np.sin(mid), 0.0]))
        bp_rows = VGroup(
            layout.label(f"2015 BP519:  {bp['i_deg']:.0f}°", font_size=19, color=P.GREEN,
                         weight="BOLD"),
            layout.label(f"a = {bp['a']:.0f} AU,  e = {bp['e']:.2f},  "
                         f"perihelion {bp['q']:.0f} AU", font_size=14, color=P.GREEN),
            layout.label(f"found by the Dark Energy Survey at H = {d['h_mag']:.1f}",
                         font_size=14, color=P.MUTED),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        bp_rows.next_to(ladder, DOWN, buff=0.35, aligned_edge=LEFT)
        cap2 = layout.caption(
            f"2015 BP519: tilted {bp['i_deg']:.0f}°, about twice as steep as any of them",
            font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(Create(new), Create(arc), FadeIn(arc_lab), FadeIn(cap2), run_time=1.8)
        self.play(FadeIn(bp_rows, lag_ratio=0.2))
        timing.hold_to_read(self, cap2, bp_rows, settle=1.2)
        self.play(FadeOut(VGroup(plane, plane_lab, sun, others, new, arc, arc_lab, ladder,
                                 bp_rows, cap2)))

        # 2. what tilts an orbit that far
        runs = d["pumped"]
        ctrl = d["control"]
        t_end = ctrl["t_myr"][-1]
        ax, labels = widgets.labeled_axes(
            [0, t_end, 10], [0, 80, 20], x_label="time (million years)",
            y_label="tilt from Planet Nine's orbit plane (°)", y_rotate=True, numbers=True,
            x_length=10.0, y_length=4.2, shift_down=-0.3)
        target = DashedLine(ax.c2p(0, d["i_deg"]), ax.c2p(t_end, d["i_deg"]), color=P.GREEN,
                            stroke_width=2)
        target_lab = layout.label(f"{d['i_deg']:.0f}°: the tilt of 2015 BP519, for scale",
                                  font_size=15, color=P.GREEN)
        target_lab.next_to(ax.c2p(0, d["i_deg"]), UP, buff=0.08, aligned_edge=LEFT).shift(
            RIGHT * 0.15)
        flat = widgets.curve(ax, ctrl["t_myr"], ctrl["i_deg"], color=P.MUTED, stroke_width=3)
        flat_lab = layout.label("no Planet Nine: the tilt never changes", font_size=15,
                                color=P.MUTED)
        flat_lab.next_to(ax.c2p(t_end, ctrl["i_deg"][-1]), DOWN, buff=0.12, aligned_edge=RIGHT)
        pumped = VGroup(*[widgets.curve(ax, r["t_myr"], r["i_deg"], color=P.ORANGE,
                                        stroke_width=3.5) for r in runs])
        pumped_lab = layout.label("with Planet Nine", font_size=15, color=P.ORANGE)
        pumped_lab.next_to(ax.c2p(t_end, max(r["i_deg"][-1] for r in runs)), UP, buff=0.1,
                           aligned_edge=RIGHT)
        p9 = d["p9"]
        cap3 = layout.caption(
            f"A scattered orbit (a = {d['particle']['a']:.0f} AU) followed for {t_end:.0f} "
            f"million years, with and without the planet", font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(Create(target), FadeIn(target_lab), FadeIn(cap3))
        self.play(Create(flat), FadeIn(flat_lab), run_time=1.2)
        self.play(*[Create(c) for c in pumped], FadeIn(pumped_lab), run_time=2.4)
        timing.hold_to_read(self, cap3, settle=0.8)
        starts = " and ".join(f"{r['start_deg']:.0f}°" for r in runs)
        cap4 = layout.caption(
            f"A {p9['mass']:.0f} Earth-mass planet at {p9['a']:.0f} AU lifts orbits starting "
            f"at {starts} past {d['pumped_i_deg']:.0f}°", font_size=22)
        self.play(FadeOut(cap3), run_time=0.4)
        self.play(FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=1.2)
        self.play(FadeOut(VGroup(ax, labels, target, target_lab, flat, flat_lab, pumped,
                                 pumped_lab, cap4)))

        # 3. and its perihelion points the way the cluster's do
        centre, rad = (-2.6, 0.05, 0.0), 1.8
        rose = ring(centre, rad, "longitude of perihelion")
        old = VGroup(*[spoke(centre, rad, o["varpi_deg"], P.GREEN, 1.8, 0.5) for o in sample])
        mine = spoke(centre, rad, bp["varpi_deg"], P.GREEN, 4.0)
        u = np.array([np.cos(np.deg2rad(bp["varpi_deg"])),
                      np.sin(np.deg2rad(bp["varpi_deg"])), 0.0])
        mine_lab = layout.label("2015 BP519", font_size=14, color=P.GREEN, weight="BOLD")
        mine_lab.next_to(np.array(centre) + u * (rad + 0.1), LEFT if u[0] < 0 else RIGHT,
                         buff=0.1).shift(UP * 0.2)
        mine_tag = VGroup(plate(mine_lab), mine_lab)
        rows = VGroup(
            layout.label(f"centre of the cluster:  {d['cluster_varpi_deg']:.0f}°", font_size=18),
            layout.label(f"2015 BP519:  {bp['varpi_deg']:.0f}°,  "
                         f"{abs(d['varpi_offset_deg']):.0f}° away", font_size=18, color=P.GREEN),
            layout.label(f"alignment of the sample:  {d['r_bar_before']:.2f}  →  "
                         f"{d['r_bar_after']:.2f}", font_size=18),
            layout.label("(0 = every direction, 1 = one direction)", font_size=14,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        rows.move_to([3.2, 0.1, 0])
        cap5 = layout.caption(
            "Its perihelion falls on the cluster's side of the sky, toward the edge of the group",
            font_size=22)
        self.play(FadeIn(rose), FadeIn(old, lag_ratio=0.1), run_time=1.2)
        self.play(FadeIn(mine), FadeIn(mine_tag), FadeIn(rows, lag_ratio=0.2), FadeIn(cap5),
                  run_time=1.4)
        timing.hold_to_read(self, cap5, rows, settle=1.2)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "The first distant orbit steep enough to need something tilting it.")
