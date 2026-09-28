"""Belyakov, Bernardinelli & Brown (2022) -- limits on the detection of Planet
Nine in the Dark Energy Survey.

DES imaged 5,000 deg² of the southern sky ten times over six years to r = 23.8,
three magnitudes deeper than ZTF. Deep enough for nearly every predicted
Planet Nine -- but only the few whose paths cross its footprint can be tested.
Reproduced in p9-2022-des: the footprint, the population scored by both the
DES and ZTF survey models, the recovery rates and the exclusion bookkeeping
are the crate's own (anim.json -> papers -> p9-2022-des).
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
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2022-des"


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


def _bars(plot, edges, counts, base, colour, opacity):
    g = VGroup()
    for k, n in enumerate(counts):
        if n <= 0:
            continue
        a = plot.p(edges[k], base[k])
        b = plot.p(edges[k + 1], base[k] + n)
        bar = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=colour, stroke_width=0.6)
        g.add(bar.set_fill(colour, opacity=opacity))
    return g


class Des2022(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pop = d["population"]
        ztf = [s for s in pop if s["p_ztf"] >= 0.5]
        new = [s for s in pop if s["p_ztf"] < 0.5 and s["p_des"] >= 0.5]
        rest = [s for s in pop if s["p_ztf"] < 0.5 and s["p_des"] < 0.5]
        self.add(paper.scene_header(CRATE))

        # 1. the predicted planets, and what ZTF already removed
        m = sky.SkyMap(width=9.6, dec_range=(-80, 80), centre=(-1.55, 0.55, 0.0), ra_centre=0.0)
        ecl, gal = m.reference_curves()
        dz = m.dots(ztf, color=P.TEAL, radius=0.028)
        dn = m.dots(new, color=P.TEAL, radius=0.028)
        dr = m.dots(rest, color=P.TEAL, radius=0.028)
        self.play(FadeIn(m), Create(ecl), Create(gal), run_time=1.0)
        cap = layout.caption(f"{len(pop)} predicted Planet Nines on tonight's sky", font_size=22)
        self.play(FadeIn(VGroup(dz, dn, dr), lag_ratio=0.02), FadeIn(cap), run_time=1.3)
        timing.hold_to_read(self, cap, settle=0.2)
        cap2 = layout.caption(
            f"ZTF had already ruled out the bright northern ones: {100 * d['ztf']:.0f}% of them",
            font_size=22)
        self.play(dz.animate.set_color(P.RED).set_opacity(0.3), FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.3)

        # 2. the DES footprint
        foot = VGroup(*[m.box(b["ra_lo"], b["ra_hi"], b["dec_lo"], b["dec_hi"], opacity=0.3,
                              stroke_width=0) for b in d["footprint"]])
        foot_lab = _readout("DES wide survey", f"{d['footprint_area_deg2']:,.0f} deg²",
                            f"{100 * d['footprint_sky_fraction']:.0f}% of the sky, r ≈ {d['depth_r']:.1f}",
                            P.PURPLE)
        foot_lab.next_to(m.frame, RIGHT, buff=0.25).align_to(m.frame, UP)
        cap3 = layout.caption(
            f"DES goes to r = {d['depth_r']:.1f}, three magnitudes deeper, but only in the south",
            font_size=22)
        self.play(FadeIn(foot), FadeIn(foot_lab), FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.3)
        cross = _readout("orbits crossing it", f"{100 * d['crossing_fraction']:.0f}%",
                         f"paper: {100 * d['published_crossing_fraction']:.0f}%", P.TEAL)
        cross.next_to(foot_lab, DOWN, buff=0.25)
        cap4 = layout.caption(
            "Those inside that DES would have linked are newly ruled out", font_size=22)
        self.play(*[x.animate.set_color(P.RED).scale(1.8) for x in dn], FadeIn(cross), FadeOut(cap3),
                  FadeIn(cap4), run_time=1.4)
        timing.hold_to_read(self, cap4, settle=0.8)
        self.play(FadeOut(VGroup(m, ecl, gal, dz, dn, dr, foot, foot_lab, cross, cap4)))

        # 3. deep enough for nearly all of them: the footprint is the limit
        r = np.array([s["r_mag"] for s in pop])
        is_z = np.array([s["p_ztf"] >= 0.5 for s in pop])
        is_n = np.array([s["p_ztf"] < 0.5 and s["p_des"] >= 0.5 for s in pop])
        edges = np.arange(16.0, 24.51, 0.5)
        nz, _ = np.histogram(r[is_z], bins=edges)
        nn, _ = np.histogram(r[is_n], bins=edges)
        na, _ = np.histogram(r, bins=edges)
        top = int(np.ceil(na.max() / 25.0) * 25)
        plot = Plot([16, 24.5], [0, top], [16, 17, 18, 19, 20, 21, 22, 23, 24],
                    list(range(0, top + 1, top // 5)), "brightness r  (fainter →)",
                    "predicted planets", centre=(-0.5, 0.45), width=9.2, height=4.2)
        b_z = _bars(plot, edges, nz, np.zeros_like(nz), P.RED, 0.35)
        b_n = _bars(plot, edges, nn, nz, P.RED, 0.9)
        b_r = _bars(plot, edges, na - nz - nn, nz + nn, P.TEAL, 0.5)
        comp = d["completeness"]
        c_line = plot.curve(comp["r_mag"], np.array(comp["fraction"]) * top, P.PURPLE,
                            stroke_width=2.5)
        c_lab = layout.label("DES detection efficiency", font_size=15, color=P.PURPLE)
        c_lab.next_to(plot.p(21.2, 0.97 * top), UP, buff=0.05)
        self.play(FadeIn(plot))
        cap5 = layout.caption("Every predicted planet by brightness", font_size=22)
        self.play(FadeIn(VGroup(b_z, b_n, b_r), lag_ratio=0.05), FadeIn(cap5), run_time=1.3)
        timing.hold_to_read(self, cap5, settle=0.2)
        cap6 = layout.caption(
            "DES could see nearly all of them; most just never cross its patch of sky",
            font_size=22)
        self.play(Create(c_line), FadeIn(c_lab), FadeOut(cap5), FadeIn(cap6))
        key = VGroup(
            layout.label("ruled out by ZTF", font_size=16, color=P.RED).set_opacity(0.6),
            layout.label("newly ruled out by DES", font_size=16, color=P.RED),
            layout.label("still possible", font_size=16, color=P.TEAL),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        rec = _readout("recovered if it crosses", f"{100 * d['recovery']:.0f}%",
                       f"paper: {100 * d['published_recovery']:.0f}%", P.PURPLE)
        uniq = _readout("new exclusion from DES", f"+{100 * d['des_unique']:.1f}%",
                        f"paper: +{100 * d['published_des_unique']:.0f}%", P.RED)
        side = VGroup(key, rec, uniq).arrange(DOWN, buff=0.3)
        side.to_edge(RIGHT, buff=0.3).align_to(plot.p(16, top), UP)
        self.play(FadeIn(key), FadeIn(rec))
        timing.hold_to_read(self, cap6, key, settle=0.4)
        cap7 = layout.caption(
            f"Together with ZTF, {100 * d['cumulative']:.0f}% of the predicted orbits are now gone",
            font_size=22)
        self.play(FadeIn(uniq), FadeOut(cap6), FadeIn(cap7))
        timing.hold_to_read(self, cap7, uniq, settle=1.0)
        self.play(FadeOut(cap7))

        layout.show_takeaway(
            self, "Deep but narrow: DES adds only a few percent to what ZTF excluded.")
