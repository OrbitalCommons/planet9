"""Naess et al. (2021) -- the Atacama Cosmology Telescope: a search for Planet 9.

At millimetre wavelengths Planet Nine's own heat is in the Rayleigh-Jeans
regime, so its flux falls only as the inverse square of distance. ACT's maps
are sensitive to 4-12 mJy at 150 GHz depending on position; wherever the
planet's predicted flux beats that limit, the null search rules it out.
Reproduced in p9-2021-act-mm: the flux curves, the limits and the reach in
mass and distance are the crate's own (anim.json -> papers -> p9-2021-act-mm).
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
    SurroundingRectangle,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2021-act-mm"


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


def _crossing(dist, flux, limit):
    """Distance where a falling flux curve drops to ``limit`` (interpolated)."""
    return float(np.interp(-limit, -np.asarray(flux), np.asarray(dist)))


class ActMm2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        dist = np.array(d["distance_au"])
        c5, c10 = d["cases"]
        lo, nom, hi = d["limit_lo_mjy"], d["limit_nominal_mjy"], d["limit_hi_mjy"]
        self.add(paper.scene_header(CRATE))

        # 1. the planet's millimetre glow against ACT's noise floor
        plot = Plot([200, 1100], [0.5, 60.0], [200, 400, 600, 800, 1000], [1, 3, 10, 30],
                    "distance from the Sun (AU)", "150 GHz flux density (mJy)", y_log=True,
                    centre=(0.35, 0.5), height=4.3)
        f10 = plot.curve(dist, c10["flux_150_mjy"], P.BLUE)
        f5 = plot.curve(dist, c5["flux_150_mjy"], P.BLUE, stroke_width=2.2)
        l10 = layout.label(f"{c10['mass_earth']:.0f} Earth masses", font_size=16, color=P.BLUE)
        l10.next_to(plot.p(960, np.interp(960, dist, c10["flux_150_mjy"])), UP, buff=0.22)
        l5 = layout.label(f"{c5['mass_earth']:.0f} Earth masses", font_size=16, color=P.BLUE)
        l5.next_to(plot.p(960, np.interp(960, dist, c5["flux_150_mjy"])), DOWN, buff=0.22)
        self.play(FadeIn(plot))
        cap = layout.caption(
            f"Its own heat at {c5['temp_k']:.0f} K: the flux falls only as 1/distance²",
            font_size=22)
        self.play(Create(f10), Create(f5), FadeIn(l10), FadeIn(l5), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.3)

        band = plot.band(200, 1100, lo, hi, P.PURPLE, opacity=0.14)
        band_lab = layout.label(f"ACT limit: {lo:.0f}-{hi:.0f} mJy across the map",
                                font_size=16, color=P.PURPLE)
        band_lab.next_to(plot.p(1100, hi), UP + LEFT, buff=0.08)
        cap2 = layout.caption("The maps rule out anything brighter than their noise allows",
                              font_size=22)
        self.play(FadeIn(band), FadeIn(band_lab), FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.2)

        t = ValueTracker(lo)
        limit = always_redraw(lambda: plot.hline(t.get_value(), P.PURPLE, stroke_width=2.4))

        def hits():
            g = VGroup()
            for c in (c5, c10):
                x = _crossing(dist, c["flux_150_mjy"], t.get_value())
                g.add(Dot(plot.p(x, t.get_value()), radius=0.08, color=P.RED))
                g.add(Line(plot.p(200, t.get_value()), plot.p(x, t.get_value()), color=P.RED,
                           stroke_width=6).set_opacity(0.55))
            return g

        reach = always_redraw(hits)
        tag = always_redraw(lambda: layout.label(
            f"limit {t.get_value():.0f} mJy: a {c5['mass_earth']:.0f} M⊕ planet is excluded to "
            f"{_crossing(dist, c5['flux_150_mjy'], t.get_value()):.0f} AU",
            font_size=18, color=P.FG).move_to(plot.p(650, 40.0)))
        cap3 = layout.caption("Deep patches reach far; shallow ones only nearby (red: excluded)",
                              font_size=22)
        self.play(FadeIn(limit), FadeIn(reach), FadeIn(tag), FadeOut(cap2), FadeIn(cap3))
        self.play(t.animate.set_value(hi), run_time=3.0)
        self.play(t.animate.set_value(lo), run_time=2.0)
        timing.hold_to_read(self, cap3, settle=0.2)
        self.remove(limit, reach, tag)
        self.play(FadeOut(VGroup(plot, f10, f5, l10, l5, band, band_lab, cap3)))

        # 2. what that rules out, in mass and distance
        rv = d["reach_vs_mass"]
        m = np.array(rv["mass_earth"])
        deep, shallow = np.array(rv["deep_au"]), np.array(rv["shallow_au"])
        plot2 = Plot([3, 15], [100, 900], [3, 5, 7, 9, 11, 13, 15],
                     [100, 300, 500, 700, 900], "Planet Nine mass (Earth masses)",
                     "distance from the Sun (AU)", centre=(-0.4, 0.5), width=8.6, height=4.3)
        sure = Polygon(*([plot2.p(x, 100) for x in (m[0], m[-1])]
                         + [plot2.p(x, y) for x, y in zip(m[::-1], shallow[::-1])]),
                       stroke_width=0).set_fill(P.RED, opacity=0.5)
        maybe = Polygon(*([plot2.p(x, y) for x, y in zip(m, shallow)]
                          + [plot2.p(x, y) for x, y in zip(m[::-1], deep[::-1])]),
                        stroke_width=0).set_fill(P.RED, opacity=0.22)
        self.play(FadeIn(plot2))
        cap4 = layout.caption("Excluded everywhere ACT looked (dark) or where its map is deep "
                              "(light)", font_size=22)
        self.play(FadeIn(sure), FadeIn(maybe), FadeIn(cap4), run_time=1.2)
        timing.hold_to_read(self, cap4, settle=0.3)

        pubs = VGroup()
        for c in (c5, c10):
            pb = c["published"]
            x = c["mass_earth"]
            bar = Line(plot2.p(x, pb["reach_lo_au"]), plot2.p(x, pb["reach_hi_au"]),
                       color=P.FG, stroke_width=3)
            caps = VGroup(*[Line(plot2.p(x, y) + LEFT * 0.1, plot2.p(x, y) + RIGHT * 0.1,
                                 color=P.FG, stroke_width=3)
                            for y in (pb["reach_lo_au"], pb["reach_hi_au"])])
            pubs.add(bar, caps)
        pub_lab = layout.label("paper's range", font_size=15, color=P.FG)
        pub_lab.next_to(plot2.p(c10["mass_earth"], c10["published"]["reach_hi_au"]), UP,
                        buff=0.1)
        wise = Dot(plot2.p(c5["mass_earth"], c5["wise_w1_reach_au"]), radius=0.08,
                   color=P.PURPLE)
        wise_lab = layout.label(f"WISE 3.4 µm, sunlight only: {c5['wise_w1_reach_au']:.0f} AU",
                                font_size=15, color=P.PURPLE)
        wise_lab.next_to(wise, RIGHT, buff=0.12)
        box = _readout(f"{c5['mass_earth']:.0f} M⊕ excluded to",
                       f"{c5['reach_shallow_au']:.0f}-{c5['reach_deep_au']:.0f} AU",
                       f"paper: {c5['published']['reach_lo_au']:.0f}-"
                       f"{c5['published']['reach_hi_au']:.0f} AU", P.RED)
        box.to_edge(RIGHT, buff=0.35).shift(UP * 1.3)
        cap5 = layout.caption(
            f"A {c5['mass_earth']:.0f} M⊕ planet gives {c5['flux_150_mjy_at_500au']:.1f} mJy at "
            f"500 AU here, {c5['published']['flux_150_mjy_at_500au']:.1f} in the paper's hotter model",
            font_size=22)
        self.play(FadeIn(pubs), FadeIn(pub_lab), FadeIn(box), FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, box, settle=0.5)
        cap6 = layout.caption("Millimetre heat reaches two to three times farther than "
                              "sunlight in WISE", font_size=22)
        self.play(FadeIn(wise), FadeIn(wise_lab), FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=0.3)
        area = _readout("searched", f"{d['area_deg2']:,.0f} deg²",
                        f"{100 * d['sky_fraction']:.0f}% of the sky", P.PURPLE)
        area.next_to(box, DOWN, buff=0.3)
        cap7 = layout.caption(
            f"Paper: about {100 * c5['published']['eliminated']:.0f}% of the parameter space "
            f"gone for {c5['mass_earth']:.0f} M⊕, {100 * c10['published']['eliminated']:.0f}% "
            f"for {c10['mass_earth']:.0f} M⊕", font_size=22)
        self.play(FadeIn(area), FadeOut(cap6), FadeIn(cap7))
        timing.hold_to_read(self, cap7, settle=1.0)
        self.play(FadeOut(cap7))

        layout.show_takeaway(
            self, "No glow above 4-12 mJy: a nearby, massive Planet Nine is ruled out.")
