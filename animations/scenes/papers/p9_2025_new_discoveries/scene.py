"""Cheng, Li & Yang (2025) -- 2017 OF201, a dwarf-planet candidate on an
extremely wide orbit (with Ammonite / 2023 KQ14, Chen et al. 2025).

Two new distant objects whose perihelia point away from the cluster that
motivates Planet Nine. Reproduced in p9-2025-new-discoveries from the JPL
orbits: the orbit paths, perihelion longitudes and the alignment statistics
before and after each addition are the crate's own
(anim.json -> papers -> p9-2025-new-discoveries).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Arrow,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Scene,
    Transform,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2025-new-discoveries"


class Dial(VGroup):
    """A compass of orbital longitude: 0° to the right, increasing anticlockwise."""

    def __init__(self, radius=1.5, centre=(0.0, 0.0, 0.0), title=None):
        super().__init__()
        self.radius = radius
        self.centre = np.array(centre, dtype=float)
        self.add(Circle(radius=radius, color=P.MUTED, stroke_width=1.5).move_to(self.centre))
        for deg in range(0, 360, 30):
            a, b = self.p(deg, radius), self.p(deg, radius + (0.12 if deg % 90 == 0 else 0.06))
            self.add(Line(a, b, color=P.MUTED, stroke_width=1.2))
        for deg in (0, 90, 180, 270):
            lab = layout.label(f"{deg}°", font_size=13, color=P.MUTED)
            self.add(lab.move_to(self.p(deg, radius + 0.36)))
        if title:
            t = layout.label(title, font_size=15, color=P.FG)
            self.add(t.move_to(self.centre + UP * (radius + 0.8)))

    def p(self, deg, r=None):
        r = self.radius if r is None else r
        t = np.deg2rad(deg)
        return self.centre + r * np.array([np.cos(t), np.sin(t), 0.0])

    def dot(self, deg, color, radius=0.075):
        return Dot(self.p(deg), radius=radius, color=color).set_z_index(3)

    def resultant(self, stat, color=P.FG):
        """The mean of the unit vectors: direction = mean longitude, length = R."""
        return Arrow(self.centre, self.p(stat["mean_varpi_deg"], self.radius * stat["r_bar"]),
                     buff=0, color=color, stroke_width=4.5,
                     max_tip_length_to_length_ratio=0.2).set_z_index(4)


def fit_paths(paths, box):
    pts = np.array([p for path in paths for p in path])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    s = min((box[1] - box[0]) / (hi[0] - lo[0]), (box[3] - box[2]) / (hi[1] - lo[1]))
    mid = 0.5 * (lo + hi)
    origin = np.array([0.5 * (box[0] + box[1]) - s * mid[0],
                       0.5 * (box[2] + box[3]) - s * mid[1], 0.0])
    return s, origin


def path_curve(path, s, origin, color, stroke_width=2.0, opacity=1.0):
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners([origin + s * np.array([x, y, 0.0]) for x, y in path])
    return m.set_stroke(opacity=opacity)


def far_point(path, s, origin):
    x, y = max(path, key=lambda p: p[0] ** 2 + p[1] ** 2)
    return origin + s * np.array([x, y, 0.0])


def stat_row(text, stat, color, weight="NORMAL"):
    return layout.label(
        f"{text}:  R = {stat['r_bar']:.2f},  p = {100 * stat['rayleigh_p']:.1f}%",
        font_size=17, color=color, weight=weight)


class NewDiscoveries2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        base, of, am = d["baseline"], d["of201"], d["ammonite"]
        before, with_of, with_both = d["stress"][0], d["stress"][1], d["stress"][3]

        self.add(paper.scene_header(CRATE))

        paths = [o["path"] for o in base] + [of["path"], am["path"]]
        s, origin = fit_paths(paths, (-6.6, 0.0, -2.5, 3.05))
        sun = orbits.sun(radius=0.06).move_to(origin)
        swarm = VGroup(*[path_curve(o["path"], s, origin, P.GREEN, 1.6, 0.75) for o in base])
        bar_au = 500.0
        bar = Line([0, 0, 0], [s * bar_au, 0, 0], color=P.MUTED, stroke_width=2)
        bar.move_to([-5.6, -2.55, 0])
        bar_lab = layout.label(f"{bar_au:.0f} AU", font_size=14, color=P.MUTED)
        bar_lab.next_to(bar, UP, buff=0.08)

        dial = Dial(radius=1.45, centre=(4.2, 0.95, 0.0), title="longitude of perihelion ϖ")
        dots = VGroup(*[dial.dot(o["varpi_deg"], P.GREEN) for o in base])
        mean = dial.resultant(before)
        rows = VGroup(
            stat_row(f"{before['n']} objects", before, P.GREEN),
            stat_row("with 2017 OF201", with_of, P.RED),
            stat_row("with Ammonite too", with_both, P.RED, weight="BOLD"),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        rows.move_to([4.2, -1.95, 0])

        # 1. the sample as it stood
        cap = layout.caption(
            f"The {before['n']} distant orbits seen from above: perihelia bunched on one side",
            font_size=22)
        self.play(FadeIn(sun), FadeIn(bar), FadeIn(bar_lab), FadeIn(dial), FadeIn(cap))
        self.play(LaggedStart(*[Create(c) for c in swarm], lag_ratio=0.1),
                  LaggedStart(*[FadeIn(x, scale=1.8) for x in dots], lag_ratio=0.1),
                  run_time=2.2)
        self.play(Create(mean), FadeIn(rows[0]))
        timing.hold_to_read(self, cap, rows[0], settle=0.8)

        # 2. 2017 OF201
        o_of = path_curve(of["path"], s, origin, P.RED, 3.0)
        of_tag = layout.label(
            f"2017 OF201\na = {of['a']:.0f} AU,  q = {of['q']:.0f} AU", font_size=16,
            color=P.RED, weight="BOLD", line_spacing=0.9)
        of_tag.move_to([-2.4, 2.7, 0])
        cap2 = layout.caption(
            f"2017 OF201 reaches {of['aphelion']:,.0f} AU and points "
            f"{of['offset_deg']:.0f}° away from the cluster", font_size=22)
        self.play(Create(o_of), FadeOut(cap), FadeIn(cap2), run_time=2.0)
        self.play(FadeIn(of_tag), FadeIn(dial.dot(of["varpi_deg"], P.RED, radius=0.1)),
                  Transform(mean, dial.resultant(with_of)), FadeIn(rows[1]))
        timing.hold_to_read(self, cap2, of_tag, rows[1], settle=1.0)

        # 3. Ammonite
        o_am = path_curve(am["path"], s, origin, P.RED, 3.0)
        am_tag = layout.label(
            f"Ammonite\na = {am['a']:.0f} AU,  q = {am['q']:.0f} AU", font_size=16,
            color=P.RED, weight="BOLD", line_spacing=0.9)
        am_tip = far_point(am["path"], s, origin)
        am_tag.move_to([-0.3, am_tip[1] + 0.3, 0])
        am_ptr = Line(am_tag.get_left() + LEFT * 0.08, am_tip, color=P.RED, stroke_width=1.5)
        am_tag = VGroup(am_tag, am_ptr)
        cap3 = layout.caption(
            f"Ammonite (Chen et al. 2025), a fourth Sedna-like orbit, points "
            f"{am['offset_deg']:.0f}° away", font_size=22)
        self.play(Create(o_am), FadeOut(cap2), FadeIn(cap3), run_time=1.6)
        self.play(FadeIn(am_tag), FadeIn(dial.dot(am["varpi_deg"], P.RED, radius=0.1)),
                  Transform(mean, dial.resultant(with_both)), FadeIn(rows[2]))
        timing.hold_to_read(self, cap3, am_tag, rows[2], settle=1.2)
        self.play(FadeOut(cap3), FadeOut(bar), FadeOut(bar_lab))

        layout.show_takeaway(
            self, f"Both newcomers point off the cluster: alignment falls from "
                  f"{before['r_bar']:.2f} to {with_both['r_bar']:.2f}.")
