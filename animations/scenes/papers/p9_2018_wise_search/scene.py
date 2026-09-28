"""Meisner et al. (2017) -- a 3pi search for Planet Nine at 3.4 microns with
WISE and NEOWISE.

The 2016 coadd search covered one patch; this one covers every part of the sky
away from the crowded galactic plane and finds nothing brighter than W1 = 16.7.
What that excludes depends on how Planet Nine shines at 3.4 µm: its own heat
fades as the square of distance, reflected sunlight as the fourth power.
Reproduced in p9-2018-wise-search: the footprint mask, the brightness curves
and the scored population of predicted planets are the crate's own
(anim.json -> papers -> p9-2018-wise-search).
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
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2018-wise-search"


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


class WiseSearch2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pop = d["population"]
        mask = d["mask"]
        depth = d["published_depth_w1"]

        self.add(paper.scene_header(CRATE))

        # 1. the footprint: everything but the galactic plane
        m = sky.SkyMap(width=10.0, dec_range=(-80, 80), centre=(-1.3, 0.55, 0.0))
        ecl, gal = m.reference_curves()
        whole = m.dec_band(-90, 90, opacity=0.16)
        plane = m.cells(mask["ra_centres"], mask["dec_centres"], mask["masked"], color=P.RED,
                        vmax=1.0, max_opacity=0.35, floor=0.5)
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), run_time=1.0)
        dots = m.dots(pop, color=P.TEAL, radius=0.028)
        cap = layout.caption(
            f"{len(pop)} Planet Nines drawn from the predicted orbits", font_size=22)
        self.play(FadeIn(dots, lag_ratio=0.02), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.5)

        cap2 = layout.caption(
            f"Four years of WISE coadds cover the sky; the crowded galactic plane "
            f"(|b| < {d['galactic_mask_deg']:.0f}°) is masked", font_size=22)
        self.play(FadeIn(whole), FadeOut(cap), FadeIn(cap2), run_time=0.9)
        self.play(FadeIn(plane), run_time=1.0)
        timing.hold_to_read(self, cap2, settle=0.4)

        tally = _readout("sky searched", f"{100 * d['sky_fraction']:.0f}%",
                         f"paper: {100 * d['published_sky_fraction']:.0f}%", P.PURPLE)
        tally.next_to(m.frame, RIGHT, buff=0.35).align_to(m.frame, UP)
        cap3 = layout.caption(
            f"Result: nothing moving like Planet Nine brighter than W1 = {depth:.1f}",
            font_size=22)
        self.play(FadeIn(tally), FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.8)
        self.play(FadeOut(VGroup(m, ecl, gal, whole, plane, dots, tally, cap3)))

        # 2. but how bright is Planet Nine at 3.4 µm?  Slide it outward.
        c = d["curve"]
        dist = np.array(c["distance_au"])
        refl = np.array(c["w1_reflected_10me"])
        lum = np.array(c["w1_luminous"])
        plot = Plot([100, 1100], [23.0, 12.0], [200, 400, 600, 800, 1000],
                    [22, 20, 18, 16, 14, 12], "distance from the Sun (AU)",
                    "W1 magnitude (brighter upward)", centre=(0.35, 0.45), height=4.4)
        sun_c = plot.curve(dist, refl, P.SUN)
        lum_c = plot.curve(dist, lum, P.BLUE)
        limit = plot.hline(depth, P.PURPLE)
        limit_lab = layout.label(f"WISE limit W1 = {depth:.1f}", font_size=15, color=P.PURPLE)
        limit_lab.next_to(limit.get_end(), UP, buff=0.06).align_to(limit.get_end(), RIGHT)
        lum_lab = layout.label("its own heat: fades as 1/d²", font_size=16, color=P.BLUE)
        lum_lab.move_to(plot.p(330, 13.1))
        sun_lab = layout.label("reflected sunlight: fades as 1/d⁴", font_size=16, color=P.SUN)
        sun_lab.next_to(plot.p(150, 22.2), RIGHT, buff=0.1)
        self.play(FadeIn(plot), Create(limit), FadeIn(limit_lab))
        cap4 = layout.caption("Two ways a 10 Earth-mass planet can shine at 3.4 µm",
                              font_size=22)
        self.play(Create(lum_c), Create(sun_c), FadeIn(lum_lab), FadeIn(sun_lab), FadeIn(cap4),
                  run_time=1.6)
        timing.hold_to_read(self, cap4, settle=0.3)

        t = ValueTracker(120.0)

        def rider(curve, colour):
            def build():
                x = t.get_value()
                y = float(np.interp(x, dist, curve))
                col = colour if y <= depth else P.RED
                return Dot(plot.p(x, y), radius=0.09, color=col)
            return always_redraw(build)

        cursor = always_redraw(lambda: plot.vline(t.get_value(), P.MUTED, stroke_width=1.2))
        readout = always_redraw(lambda: layout.label(
            f"planet at d = {t.get_value():.0f} AU", font_size=18, color=P.FG).move_to(
                plot.p(900, 12.9)))
        r_sun, r_lum = rider(refl, P.SUN), rider(lum, P.BLUE)
        cap5 = layout.caption("Move it outward: below the dashed line WISE cannot see it",
                              font_size=22)
        self.play(FadeIn(cursor), FadeIn(r_sun), FadeIn(r_lum), FadeIn(readout),
                  FadeOut(cap4), FadeIn(cap5))
        self.play(t.animate.set_value(1050.0), run_time=5.0, rate_func=lambda a: a)
        self.remove(cursor, r_sun, r_lum, readout)

        reach = VGroup()
        for x, side in ((d["reach_reflected_10me_au"], UP + RIGHT), (d["reach_luminous_au"], DOWN)):
            dot = Dot(plot.p(x, depth), radius=0.08, color=P.GREEN)
            lab = layout.label(f"{x:.0f} AU", font_size=16, color=P.GREEN)
            lab.next_to(dot, side, buff=0.1)
            reach.add(dot, lab)
        lo, hi = d["published_reach_au"]
        cap6 = layout.caption(
            f"Reach: {d['reach_reflected_10me_au']:.0f} AU in sunlight, "
            f"{d['reach_luminous_au']:.0f} AU if self-luminous (paper: {lo:.0f}-{hi:.0f} AU)",
            font_size=22)
        self.play(FadeIn(reach), FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=0.8)
        self.play(FadeOut(VGroup(plot, sun_c, lum_c, limit, limit_lab, lum_lab, sun_lab, reach,
                                 cap6)))

        # 3. what that excludes among the predicted planets
        r_lum_au, r_sun_au = d["reach_luminous_au"], d["reach_reflected_10me_au"]
        dd = np.array([s["dist_au"] for s in pop])
        hit = np.array([bool(s["in_luminous_reach"]) for s in pop])
        edges = np.arange(150.0, 1000.01, 50.0)
        n_hit, _ = np.histogram(dd[hit], bins=edges)
        n_all, _ = np.histogram(dd, bins=edges)
        top = int(np.ceil(n_all.max() / 20.0) * 20)
        hplot = Plot([150, 1000], [0, top], [200, 400, 600, 800, 1000],
                     list(range(0, top + 1, top // 4)), "distance of the predicted planet today (AU)",
                     "predicted planets", centre=(-0.6, 0.5), width=9.0, height=4.2)
        bars_hit = _bars(hplot, edges, n_hit, np.zeros_like(n_hit), P.TEAL)
        bars_rest = _bars(hplot, edges, n_all - n_hit, n_hit, P.TEAL)
        self.play(FadeIn(hplot))
        cap7 = layout.caption("Where the predicted planets are today", font_size=22)
        self.play(FadeIn(VGroup(*bars_hit, *bars_rest), lag_ratio=0.05), FadeIn(cap7),
                  run_time=1.4)
        timing.hold_to_read(self, cap7, settle=0.3)

        m_sun = plot_marker(hplot, r_sun_au, top, f"sunlight reach {r_sun_au:.0f} AU", P.SUN)
        m_lum = plot_marker(hplot, r_lum_au, top, f"self-luminous reach {r_lum_au:.0f} AU",
                            P.BLUE)
        f_lum, f_sun = d["fraction_luminous"], d["fraction_reflected"]
        box_lum = _readout("if self-luminous", f"{100 * f_lum:.0f}%", None, P.RED)
        box_sun = _readout("if only sunlit", f"{100 * f_sun:.1f}%", None,
                           P.SUN)
        head = layout.label("ruled out", font_size=18, color=P.FG)
        boxes = VGroup(head, box_lum, box_sun).arrange(DOWN, buff=0.22)
        boxes.to_edge(RIGHT, buff=0.45).align_to(hplot.p(0, 0.95 * top), UP)
        cap8 = layout.caption(
            "Inside the footprint and the reach (red), the planet is ruled out", font_size=22)
        self.play(Create(m_sun), Create(m_lum), FadeOut(cap7), FadeIn(cap8))
        self.play(*[b.animate.set_fill(P.RED, opacity=0.8) for b in bars_hit], FadeIn(boxes))
        timing.hold_to_read(self, cap8, box_lum, box_sun, settle=1.0)
        self.play(FadeOut(cap8))

        layout.show_takeaway(
            self, "Most predicted orbits are excluded only if Planet Nine glows at 3.4 µm.")


def _readout(title, value, note, colour):
    """A boxed number with a readable title (and an optional comparison)."""
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


def _bars(plot, edges, counts, base, colour):
    bars = []
    for k, n in enumerate(counts):
        if n <= 0:
            continue
        a = plot.p(edges[k], base[k])
        b = plot.p(edges[k + 1], base[k] + n)
        bar = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=colour, stroke_width=0.6)
        bars.append(bar.set_fill(colour, opacity=0.55))
    return bars


def plot_marker(plot, x, top, text, colour):
    line = plot.vline(x, colour, y_from=0, y_to=top)
    lab = layout.label(text, font_size=15, color=colour)
    lab.next_to(line.get_end(), UP, buff=0.08)
    return VGroup(line, lab)
