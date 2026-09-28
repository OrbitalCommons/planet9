"""Linder & Mordasini (2016) -- evolution and magnitudes of candidate Planet Nine.

Cooling models of a small ice giant give its present temperature and radius,
and from those its brightness in every band: faint in reflected sunlight, which
falls as the fourth power of distance, but bright near 20 µm from its own heat.
Reproduced in p9-2016-linder-evolution with the workspace reflected-light
photometry and a single-temperature thermal model; the curves are the crate's
own and the open markers are the paper's tabulated magnitudes
(anim.json -> papers -> p9-2016-linder-evolution).
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
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2016-linder-evolution"


class Plot(VGroup):
    """Axes with the origin at the lower-left corner of the data range, optional
    log10 scales, and tick labels and axis titles in the house style."""

    def __init__(self, x_range, y_range, x_ticks, y_ticks, x_label, y_label, width=10.0,
                 height=4.2, centre=(0.35, 0.55), x_log=False, y_log=False, x_fmt="{:g}",
                 y_fmt="{:g}", **kwargs):
        super().__init__(**kwargs)
        self.x_log, self.y_log = x_log, y_log
        self.x0, self.x1 = (self._s(v, x_log) for v in x_range)
        self.y0, self.y1 = (self._s(v, y_log) for v in y_range)
        self.left, self.bottom = centre[0] - width / 2, centre[1] - height / 2
        self.w, self.h = width, height
        frame = VGroup(
            Line([self.left, self.bottom, 0], [self.left + width, self.bottom, 0]),
            Line([self.left, self.bottom, 0], [self.left, self.bottom + height, 0]),
        ).set_stroke(P.MUTED, width=2)
        marks = VGroup()
        for v in x_ticks:
            at = self.p(v, y_range[0])
            marks.add(Line(at, at + DOWN * 0.08, color=P.MUTED, stroke_width=1.5))
            marks.add(layout.label(self._text(x_fmt, v), font_size=15, color=P.MUTED)
                      .next_to(at, DOWN, buff=0.14))
        widest = 0.0
        for v in y_ticks:
            at = self.p(x_range[0], v)
            lab = layout.label(self._text(y_fmt, v), font_size=15, color=P.MUTED)
            lab.next_to(at, LEFT, buff=0.14)
            widest = max(widest, lab.width)
            marks.add(Line(at, at + LEFT * 0.08, color=P.MUTED, stroke_width=1.5), lab)
        xl = layout.label(x_label, font_size=19, color=P.FG)
        xl.move_to([centre[0], self.bottom - 0.66, 0])
        yl = layout.label(y_label, font_size=17, color=P.FG).rotate(np.pi / 2)
        yl.move_to([self.left - widest - 0.45, centre[1], 0])
        self.add(frame, marks, xl, yl)

    @staticmethod
    def _s(v, log):
        return float(np.log10(v)) if log else float(v)

    @staticmethod
    def _text(fmt, v):
        return fmt(v) if callable(fmt) else fmt.format(v)

    def p(self, x, y):
        """Scene point of the data point (x, y), clamped to the plot area."""
        u = (self._s(x, self.x_log) - self.x0) / (self.x1 - self.x0)
        v = (self._s(y, self.y_log) - self.y0) / (self.y1 - self.y0)
        return np.array([self.left + self.w * min(max(u, 0.0), 1.0),
                         self.bottom + self.h * min(max(v, 0.0), 1.0), 0.0])

    def inside(self, x, y):
        sx, sy = self._s(x, self.x_log), self._s(y, self.y_log)
        return (min(self.x0, self.x1) <= sx <= max(self.x0, self.x1)
                and min(self.y0, self.y1) <= sy <= max(self.y0, self.y1))

    def curve(self, xs, ys, color, stroke_width=3.0):
        """Polyline through the data points that fall inside the plot area."""
        runs, cur = [], []
        for x, y in zip(xs, ys):
            if self.inside(x, y):
                cur.append(self.p(x, y))
            elif cur:
                runs.append(cur)
                cur = []
        if cur:
            runs.append(cur)
        g = VGroup()
        for pts in runs:
            if len(pts) > 1:
                g.add(VMobject(color=color, stroke_width=stroke_width)
                      .set_points_as_corners(pts))
        return g

    def hline(self, y, color, x_from=None, x_to=None, stroke_width=1.6):
        a = self.p(self._x(x_from, self.x0), y)
        b = self.p(self._x(x_to, self.x1), y)
        return DashedLine(a, b, color=color, stroke_width=stroke_width)

    def vline(self, x, color, y_from=None, y_to=None, stroke_width=1.6):
        a = self.p(x, self._y(y_from, self.y0))
        b = self.p(x, self._y(y_to, self.y1))
        return DashedLine(a, b, color=color, stroke_width=stroke_width)

    def band(self, x_from, x_to, y_from, y_to, color, opacity=0.16):
        a, b = self.p(x_from, y_from), self.p(x_to, y_to)
        r = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=color, stroke_width=0)
        return r.set_fill(color, opacity=opacity)

    def _x(self, v, scaled):
        return (10 ** scaled if self.x_log else scaled) if v is None else v

    def _y(self, v, scaled):
        return (10 ** scaled if self.y_log else scaled) if v is None else v

def spread_vertically(labels, gap):
    """Push labels apart so neighbours are at least ``gap`` from centre to centre."""
    ordered = sorted(labels, key=lambda m: m.get_center()[1])
    for below, above in zip(ordered, ordered[1:]):
        short = gap - (above.get_center()[1] - below.get_center()[1])
        if short > 0:
            above.shift(UP * short)


def published_marker(point):
    return Circle(radius=0.08, color=P.FG, stroke_width=2.2).move_to(point)


class LinderEvolution2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        dist = np.array(d["distance_au"])
        curves = d["curves"]
        aph = d["aphelion"]

        self.add(paper.scene_header(CRATE))

        def panel(band, y_range, y_ticks, y_label, colour):
            plot = Plot([200, 1200], y_range, [200, 400, 600, 800, 1000, 1200], y_ticks,
                        "distance from the Sun (AU)", y_label)
            lines, names = VGroup(), VGroup()
            for c in curves:
                lines.add(plot.curve(dist, c[band], colour, stroke_width=2.6))
                name = layout.label(f"{c['mass_earth']:.0f} M⊕", font_size=15, color=colour)
                name.next_to(plot.p(dist[-1], c[band][-1]), RIGHT, buff=0.1)
                names.add(name)
            spread_vertically(names, 0.26)
            orbit = VGroup()
            for au, text in ((d["perihelion_au"], "perihelion"), (d["nominal_au"], "a"),
                             (d["aphelion_au"], "aphelion")):
                line = plot.vline(au, P.MUTED, stroke_width=1.2)
                lab = layout.label(f"{text} {au:.0f} AU", font_size=14, color=P.MUTED)
                lab.next_to(line.get_end(), UP, buff=0.06)
                orbit.add(VGroup(line, lab))
            return plot, lines, names, orbit

        # 1. reflected sunlight
        plot, lines, names, orbit = panel("v_mag", [25.5, 16.5], [24, 22, 20, 18],
                                          "V magnitude (brighter upward)", P.SUN)
        self.play(FadeIn(plot), FadeIn(orbit))
        cap = layout.caption("In reflected sunlight Planet Nine fades as distance to the fourth power",
                             font_size=22)
        self.play(Create(lines), FadeIn(names), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.6)

        marks = VGroup(published_marker(plot.p(d["nominal_au"], d["v_10me_700au_published"])))
        for row in aph:
            marks.add(published_marker(plot.p(d["aphelion_au"], row["v_published"])))
        worst = max(abs(r["v_mag"] - r["v_published"]) for r in aph)
        key = VGroup(published_marker([0, 0, 0]),
                     layout.label("paper's cooling models", font_size=15, color=P.FG))
        key.arrange(RIGHT, buff=0.12).move_to(plot.p(450, 24.6))
        cap2 = layout.caption(
            f"10 Earth masses at {d['nominal_au']:.0f} AU: V = {d['v_10me_700au']:.1f} here, "
            f"{d['v_10me_700au_published']:.1f} in the paper; all within {worst:.1f} mag",
            font_size=22)
        self.play(FadeIn(marks), FadeIn(key), FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.8)
        self.play(FadeOut(VGroup(plot, lines, names, orbit, marks, key, cap2)))

        # 2. its own heat, near 20 µm
        plot, lines, names, orbit = panel("q_mag", [17.0, 3.0], [16, 14, 12, 10, 8, 6, 4],
                                          "Q magnitude at 20 µm (brighter upward)", P.BLUE)
        self.play(FadeIn(plot), FadeIn(orbit))
        cap3 = layout.caption(
            f"At 20 µm it shines by its own heat: {d['t_eff_10me_k']:.0f} K for 10 Earth masses",
            font_size=22)
        self.play(Create(lines), FadeIn(names), FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, settle=0.6)

        marks = VGroup(published_marker(plot.p(d["nominal_au"], d["q_10me_700au_published"])))
        for row in aph:
            marks.add(published_marker(plot.p(d["aphelion_au"], row["q_published"])))
        gaps = [r["q_mag"] - r["q_published"] for r in aph]
        key = VGroup(published_marker([0, 0, 0]),
                     layout.label("paper's cooling models", font_size=15, color=P.FG))
        key.arrange(RIGHT, buff=0.12).move_to(plot.p(800, 4.0))
        cap4 = layout.caption(
            f"The single-temperature model here is {min(gaps):.1f}-{max(gaps):.1f} mag fainter "
            "than the paper's", font_size=22)
        self.play(FadeIn(marks), FadeIn(key), FadeOut(cap3), FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=0.8)

        heavy = aph[-1]
        cap5 = layout.caption(
            f"At {heavy['mass_earth']:.0f} Earth masses it reaches Q = {heavy['q_published']:.1f}: "
            "past infrared surveys would have seen it", font_size=22)
        ring = Circle(radius=0.2, color=P.RED, stroke_width=3).move_to(
            plot.p(d["aphelion_au"], heavy["q_published"]))
        self.play(Create(ring), FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=0.8)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Faint in sunlight, bright in its own heat: look in the far-infrared.")
