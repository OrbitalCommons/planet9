"""Meisner et al. (2016) -- searching for Planet Nine with coadded WISE and
NEOWISE-Reactivation images.

Single WISE exposures reach W1 = 15.3; coadding the dozen exposures of a day
reaches 16.66, which for the most luminous model atmosphere of Fortney et al.
(2016) moves the reach from 430 to 800 AU. The search covered the patch of sky
that the Cassini ranging analysis favoured and found nothing. Reproduced in
p9-2016-wise-coadd: depths, brightness curves and reaches are the crate's own
(anim.json -> papers -> p9-2016-wise-coadd).
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
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2016-wise-coadd"


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

class WiseCoadd2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        region = d["region"]

        self.add(paper.scene_header(CRATE))

        # 1. where they looked
        m = sky.SkyMap(width=12.0, dec_range=(-75, 75), centre=(0.0, 0.1, 0.0))
        ecl, gal = m.reference_curves()
        patch = m.box(region["ra_lo"], region["ra_hi"], region["dec_lo"], region["dec_hi"],
                      opacity=0.3)
        patch_lab = layout.label("searched", font_size=16, color=P.PURPLE)
        patch_lab.next_to(patch, UP, buff=0.1)
        key = m.legend([("ecliptic", P.ORANGE), ("galactic plane ±10°", P.PURPLE)])
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), FadeIn(key), run_time=1.0)
        cap = layout.caption(
            f"The patch Cassini ranging favoured: {region['area_deg2']:,.0f} deg² here, "
            f"{region['published_area_deg2']:,.0f} deg² in the paper", font_size=22)
        self.play(FadeIn(patch), FadeIn(patch_lab), FadeOut(key), FadeIn(cap))
        timing.hold_to_read(self, cap, settle=0.8)
        self.play(FadeOut(VGroup(m, ecl, gal, patch, patch_lab, cap)))

        # 2. how stacking deepens the search
        dv = d["depth_vs_frames"]
        plot = Plot([0, 40], [15.0, 17.5], [1, 10, 20, 30, 40], [15.0, 15.5, 16.0, 16.5, 17.0, 17.5],
                    "exposures coadded", "limiting W1 magnitude (deeper upward)",
                    y_fmt="{:.1f}")
        curve = plot.curve(dv["frames"], dv["w1_depth"], P.PURPLE)
        single = Dot(plot.p(1, d["single_depth"]), radius=0.07, color=P.PURPLE)
        single_lab = layout.label(f"one exposure: W1 = {d['single_depth']:.1f}", font_size=16,
                                  color=P.PURPLE)
        single_lab.next_to(single, RIGHT, buff=0.2).shift(DOWN * 0.1)
        self.play(FadeIn(plot), FadeIn(single), FadeIn(single_lab))
        cap2 = layout.caption("Every tenfold increase in exposures gains 1.25 magnitudes",
                              font_size=22)
        self.play(Create(curve), FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, settle=0.4)

        target = plot.hline(d["published_coadd_depth"], P.FG)
        target_lab = layout.label(
            f"paper: 90% complete to W1 = {d['published_coadd_depth']:.2f}", font_size=16,
            color=P.FG)
        target_lab.next_to(target.get_end(), UP, buff=0.08).align_to(target.get_end(), RIGHT)
        day = Dot(plot.p(d["frames_per_coadd"], d["coadd_depth"]), radius=0.08, color=P.GREEN)
        day_lab = layout.label(f"{d['frames_per_coadd']:.0f} exposures", font_size=16,
                               color=P.GREEN)
        day_lab.next_to(day, DOWN, buff=0.15).shift(RIGHT * 0.6)
        cap3 = layout.caption(
            f"A day of WISE passes, about {d['frames_per_coadd']:.0f} exposures, "
            f"reaches the published depth", font_size=22)
        self.play(Create(target), FadeIn(target_lab), FadeIn(day), FadeIn(day_lab),
                  FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.8)
        self.play(FadeOut(VGroup(plot, curve, single, single_lab, target, target_lab, day,
                                 day_lab, cap3)))

        # 3. how far that reaches
        dist = np.array(d["distance_au"])
        refl = d["w1_reflected"]
        plot2 = Plot([100, 1000], [21.0, 12.0], [200, 400, 600, 800, 1000],
                     [20, 18, 16, 14, 12], "distance from the Sun (AU)",
                     "W1 magnitude (brighter upward)")
        lum = plot2.curve(dist, d["w1_luminous"], P.BLUE)
        sun = plot2.curve(dist, refl["default"], P.SUN)
        lum_lab = layout.label("brightest model atmosphere", font_size=16, color=P.BLUE)
        lum_lab.next_to(plot2.p(300, np.interp(300, dist, d["w1_luminous"])), UP, buff=0.3)
        lum_lab.shift(RIGHT * 1.2)
        sun_lab = layout.label(
            f"reflected sunlight only (albedo {refl['albedo_default']:.2f})", font_size=16,
            color=P.SUN)
        sun_lab.next_to(plot2.p(600, 19.6), RIGHT, buff=0.0)
        self.play(FadeIn(plot2))
        cap4 = layout.caption(
            f"How bright a {d['mass_earth']:.0f} Earth-mass planet is at 3.4 µm depends on "
            "its atmosphere", font_size=22)
        self.play(Create(lum), Create(sun), FadeIn(lum_lab), FadeIn(sun_lab), FadeIn(cap4),
                  run_time=1.5)
        timing.hold_to_read(self, cap4, settle=0.5)

        limits = VGroup()
        for depth, text in ((d["single_depth"], "single exposures"),
                            (d["coadd_depth"], "coadds")):
            line = plot2.hline(depth, P.PURPLE)
            lab = layout.label(text, font_size=15, color=P.PURPLE)
            lab.next_to(line.get_end(), UP, buff=0.06).align_to(line.get_end(), RIGHT)
            limits.add(line, lab)
        reach = VGroup()
        for key, depth, corner in (("reach_luminous_single_au", d["single_depth"], DOWN + LEFT),
                                   ("reach_luminous_coadd_au", d["coadd_depth"], DOWN + LEFT),
                                   ("reach_reflected_single_au", d["single_depth"], UP + RIGHT),
                                   ("reach_reflected_coadd_au", d["coadd_depth"], UP + RIGHT)):
            dot = Dot(plot2.p(d[key], depth), radius=0.07, color=P.GREEN)
            lab = layout.label(f"{d[key]:.0f} AU", font_size=15, color=P.GREEN)
            lab.next_to(dot, corner, buff=0.06)
            reach.add(dot, lab)
        cap5 = layout.caption(
            f"Coadds push the reach from {d['reach_luminous_single_au']:.0f} to "
            f"{d['reach_luminous_coadd_au']:.0f} AU (paper: "
            f"{d['published_single_reach_au']:.0f} to {d['published_coadd_reach_au']:.0f})",
            font_size=22)
        self.play(FadeIn(limits), FadeIn(reach), FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Nothing to W1 = 16.66: a self-luminous Planet Nine is not in that patch.")
