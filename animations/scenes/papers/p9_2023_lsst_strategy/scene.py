"""Schwamb et al. (2023) -- tuning the LSST observing strategy for solar system
science.

Rubin's scheduler can trade depth, revisits and sky coverage against each
other. Scored against the predicted Planet Nines, only one of those knobs
matters: every predicted planet inside the footprint is bright enough and
revisited often enough to be linked, so what LSST can find is set by how far
north the survey reaches. Reproduced in p9-2023-lsst-strategy: the footprint,
the scored population and the strategy scans are the crate's own
(anim.json -> papers -> p9-2023-lsst-strategy).
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
    Rectangle,
    Scene,
    SurroundingRectangle,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2023-lsst-strategy"


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


def _readout(title, value, note, colour):
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


class LsstStrategy2023(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        fp = d["footprint"]
        pop = d["population"]
        self.add(paper.scene_header(CRATE))

        # 1. the baseline footprint against the predicted planets
        m = sky.SkyMap(width=9.6, dec_range=(-80, 80), centre=(-1.55, 0.55, 0.0))
        ecl, gal = m.reference_curves()
        found = [s for s in pop if s["p_discover"] >= 0.5]
        missed = [s for s in pop if s["p_discover"] < 0.5]
        df = m.dots(found, color=P.TEAL, radius=0.028)
        dm = m.dots(missed, color=P.TEAL, radius=0.028)
        self.play(FadeIn(m), Create(ecl), Create(gal), run_time=1.0)
        cap = layout.caption(f"{len(pop)} predicted Planet Nines", font_size=22)
        self.play(FadeIn(VGroup(df, dm), lag_ratio=0.02), FadeIn(cap), run_time=1.2)
        timing.hold_to_read(self, cap, settle=0.2)
        band = m.dec_band(fp["dec_lo"], fp["dec_hi"], opacity=0.2)
        foot = _readout("LSST baseline", f"dec {fp['dec_lo']:.0f}° to +{fp['dec_hi']:.0f}°",
                        f"galactic plane (|b| < {fp['galactic_lat_min_deg']:.0f}°) thin",
                        P.PURPLE)
        foot.next_to(m.frame, RIGHT, buff=0.25).align_to(m.frame, UP)
        cap2 = layout.caption(
            f"Rubin sees the southern sky to r = {d['single_visit_depth']:.1f} per visit, "
            f"some {d['visits_per_field']} usable visits in ten years", font_size=22)
        self.play(FadeIn(band), FadeIn(foot), FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.2)
        frac = _readout("predicted orbits LSST would find", f"{100 * d['fraction']:.0f}%",
                        None, P.TEAL)
        frac.next_to(foot, DOWN, buff=0.3)
        cap3 = layout.caption(
            "Inside the footprint every one is linked; outside, none", font_size=22)
        self.play(dm.animate.set_color(P.MUTED).set_opacity(0.4), FadeIn(frac), FadeOut(cap2),
                  FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.8)
        self.play(FadeOut(VGroup(m, ecl, gal, band, df, dm, foot, frac, cap3)))

        # 2. turn each knob: which one moves the answer?
        knobs = [
            ("by_depth", "depth_r", "depth per visit (r)", [21, 22, 23, 24, 25],
             d["single_visit_depth"], "{:g}"),
            ("by_visits", "visits", "visits needed to link", [1, 10, 20, 30],
             d["visits_for_linking"], "{:g}"),
            ("by_dec_limit", "dec_max_deg", "northern limit (dec)", [-20, 0, 20, 40],
             fp["dec_hi"], "{:g}°"),
        ]
        plots, curves, bases = VGroup(), VGroup(), VGroup()
        for k, (key, xk, xl, ticks, base, fmt) in enumerate(knobs):
            xs = np.array(d[key][xk], dtype=float)
            ys = np.array(d[key]["fraction"])
            pl = Plot([ticks[0], ticks[-1]], [0, 0.8], ticks, [0, 0.2, 0.4, 0.6, 0.8], xl,
                      "found" if k == 0 else "", centre=(-4.35 + 4.45 * k, 0.75), width=3.4,
                      height=3.0, x_fmt=fmt, y_fmt=lambda v: f"{100 * v:.0f}%")
            plots.add(pl)
            curves.add(pl.curve(xs, ys, P.TEAL, stroke_width=3.5))
            bases.add(Dot(pl.p(base, float(np.interp(base, xs, ys))), radius=0.07,
                          color=P.FG))
        self.play(FadeIn(plots))
        cap4 = layout.caption("Change one choice at a time and re-score the planets",
                              font_size=22)
        self.play(FadeIn(bases), FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=0.2)

        caps = [
            f"Deeper or shallower images: {100 * min(d['by_depth']['fraction']):.0f}-"
            f"{100 * max(d['by_depth']['fraction']):.0f}%",
            "More or fewer revisits: no change at all",
            f"How far north it looks: {100 * min(d['by_dec_limit']['fraction']):.0f}-"
            f"{100 * max(d['by_dec_limit']['fraction']):.0f}%",
        ]
        prev = cap4
        for k, (key, xk, _, _, _, _) in enumerate(knobs):
            xs = np.array(d[key][xk], dtype=float)
            ys = np.array(d[key]["fraction"])
            t = ValueTracker(float(xs[0]))
            pl = plots[k]
            rider = always_redraw(lambda pl=pl, xs=xs, ys=ys, t=t: Dot(
                pl.p(t.get_value(), float(np.interp(t.get_value(), xs, ys))), radius=0.09,
                color=P.TEAL))
            c = layout.caption(caps[k], font_size=22)
            self.play(FadeIn(rider), FadeOut(prev), FadeIn(c), run_time=0.5)
            self.play(Create(curves[k]), t.animate.set_value(float(xs[-1])), run_time=2.2,
                      rate_func=lambda a: a)
            self.remove(rider)
            timing.hold_to_read(self, c, settle=0.1)
            prev = c
        k0, k1 = knobs[2][3][0], knobs[2][3][-1]
        corner_a, corner_b = plots[2].p(k0, 0.0), plots[2].p(k1, 0.8)
        box = Rectangle(width=corner_b[0] - corner_a[0] + 0.9,
                        height=corner_b[1] - corner_a[1] + 1.25, color=P.TEAL,
                        stroke_width=2.5).move_to(
            0.5 * (corner_a + corner_b) + np.array([-0.2, -0.3, 0]))
        cap5 = layout.caption(
            f"Paper: most cadences shift the metrics by under "
            f"{100 * d['published_typical_spread']:.0f}%; for Planet Nine, sky coverage is what counts",
            font_size=22)
        self.play(Create(box), FadeOut(prev), FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Rubin will reach about half the predicted orbits; the footprint sets the rest.")
