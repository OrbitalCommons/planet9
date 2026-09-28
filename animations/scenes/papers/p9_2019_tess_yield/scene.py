"""Payne, Holman & Pál (2019) -- a TESS search for distant solar system objects:
yield estimates.

TESS stares at each sector for 27 days. A single 30-minute full-frame image
reaches only I ~ 18, but shifting ~1,300 of them along a trial orbit and adding
them up makes a slow mover pile up in one place while the noise averages down,
reaching I ~ 22. The yield note asks how much of the plausible Planet Nine box
that depth covers. Reproduced in p9-2019-tess-yield: the stacking law, the
brightness curves, the detectable fractions and the drift limit are the
crate's own (anim.json -> papers -> p9-2019-tess-yield).
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
    Square,
    SurroundingRectangle,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2019-tess-yield"


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


# Schematic frame contents for the shift-and-stack explainer: fixed star
# positions (frame units) and the mover's step per frame. Not data.
STARS = [(-0.45, 0.38), (0.28, 0.52), (0.5, -0.12), (-0.2, -0.45), (0.05, 0.05), (-0.55, -0.1)]
MOVER0, STEP = (-0.42, 0.2), (0.22, -0.1)


def _frame(k, size=1.9):
    box = Square(side_length=size, color=P.MUTED, stroke_width=1.5).set_fill("#16171f", 1.0)
    stars = VGroup(*[Dot([x * size / 1.4, y * size / 1.4, 0], radius=0.05, color=P.FG)
                     for x, y in STARS])
    mx, my = MOVER0[0] + k * STEP[0], MOVER0[1] + k * STEP[1]
    mover = Dot([mx * size / 1.4, my * size / 1.4, 0], radius=0.05, color=P.BLUE)
    mover.set_opacity(0.45)
    return VGroup(box, stars, mover)


def _readout(title, value, note, colour):
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


class TessYield2019(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        bx = d["box"]
        self.add(paper.scene_header(CRATE))

        # 1. the trick: shift the frames along the motion, then add them up
        n = 4
        frames = VGroup(*[_frame(k) for k in range(n)]).arrange(RIGHT, buff=0.45)
        frames.move_to(UP * 1.35)
        stamps = VGroup(*[layout.label(f"day {7 * k}", font_size=15, color=P.MUTED)
                          .next_to(frames[k], DOWN, buff=0.12) for k in range(n)])
        cap = layout.caption("In one TESS image a distant planet is lost in the noise",
                             font_size=22)
        self.play(FadeIn(frames, lag_ratio=0.15), FadeIn(stamps), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.3)

        size = frames[0][0].width
        target = DOWN * 1.35
        shifted = VGroup(*[f.copy() for f in frames])
        for f in shifted:
            f[0].set_fill(opacity=0.0).set_stroke(opacity=0.3)
            f[1].set_opacity(0.4)
        cap2 = layout.caption(
            "Slide each frame back along a trial orbit and add them: the planet piles up",
            font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2))
        self.play(*[f.animate.move_to(target - k * np.array([STEP[0], STEP[1], 0]) * size / 1.4)
                    for k, f in enumerate(shifted)], run_time=1.8)
        stacked = Dot(shifted[0][2].get_center(), radius=0.09, color=P.BLUE)
        star_note = layout.label("stars smear out", font_size=16, color=P.FG)
        star_note.next_to(shifted, RIGHT, buff=0.5).shift(UP * 0.3)
        p9_note = layout.label("the mover adds up", font_size=16, color=P.BLUE)
        p9_note.next_to(star_note, DOWN, buff=0.25, aligned_edge=LEFT)
        self.play(FadeIn(stacked, scale=2.0), FadeIn(star_note), FadeIn(p9_note))
        timing.hold_to_read(self, cap2, settle=0.6)
        self.play(FadeOut(VGroup(frames, stamps, shifted, stacked, star_note, p9_note, cap2)))

        # 2. how deep: noise falls as 1/sqrt(N)
        dv = d["depth_vs_frames"]
        lf, dep = np.array(dv["log_frames"]), np.array(dv["depth"])
        plot = Plot([1, 2000], [17.5, 23.0], [1, 10, 100, 1000], [18, 19, 20, 21, 22, 23],
                    "images stacked", "limiting magnitude I (deeper upward)", x_log=True,
                    centre=(0.35, 0.5), height=4.2)
        curve = plot.curve(10 ** lf, dep, P.PURPLE)
        self.play(FadeIn(plot))
        t = ValueTracker(0.0)
        rider = always_redraw(lambda: Dot(
            plot.p(10 ** t.get_value(), float(np.interp(t.get_value(), lf, dep))), radius=0.08,
            color=P.PURPLE))
        tag = always_redraw(lambda: layout.label(
            f"{10 ** t.get_value():,.0f} images: I = {np.interp(t.get_value(), lf, dep):.1f}",
            font_size=18, color=P.FG).move_to(plot.p(8, 22.4)))
        cap3 = layout.caption("Each tenfold more images reaches 1.25 magnitudes fainter",
                              font_size=22)
        self.play(FadeIn(rider), FadeIn(tag), FadeIn(cap3))
        self.play(Create(curve), t.animate.set_value(lf[-1]), run_time=4.0,
                  rate_func=lambda a: a)
        timing.hold_to_read(self, cap3, settle=0.2)
        pub = plot.band(1, 2000, d["published_depth"] - d["published_depth_sigma"],
                        d["published_depth"] + d["published_depth_sigma"], P.FG, opacity=0.08)
        pub_lab = layout.label(
            f"paper: I = {d['published_depth']:.1f} ± {d['published_depth_sigma']:.1f}",
            font_size=16, color=P.FG)
        pub_lab.next_to(plot.p(2000, d["published_depth"] - d["published_depth_sigma"]),
                        DOWN + LEFT, buff=0.1)
        gain = d["stacked_depth"] - d["single_frame_depth"]
        cap4 = layout.caption(
            f"A full sector of {d['frames_per_sector']:,.0f} images: I = {d['stacked_depth']:.1f}, "
            f"{10 ** (0.4 * gain):.0f} times fainter than one image reaches", font_size=22)
        self.play(FadeIn(pub), FadeIn(pub_lab), FadeOut(cap3), FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=0.6)
        self.remove(rider, tag)
        self.play(FadeOut(VGroup(plot, curve, pub, pub_lab, cap4)))

        # 3. is that deep enough for Planet Nine?
        dist = np.array(d["distance_au"])
        vlo, vhi = np.array(d["v_mass_lo"]), np.array(d["v_mass_hi"])
        plot2 = Plot([200, 1000], [24.0, 16.0], [200, 400, 600, 800, 1000],
                     [24, 22, 20, 18, 16], "distance from the Sun (AU)",
                     "reflected-light magnitude (brighter upward)", centre=(0.35, 0.5),
                     height=4.2)
        pts = [plot2.p(x, y) for x, y in zip(dist, vhi)] + \
              [plot2.p(x, y) for x, y in zip(dist[::-1], vlo[::-1])]
        band = Polygon(*pts, stroke_width=0).set_fill(P.TEAL, opacity=0.45)
        band_lab = layout.label(f"Planet Nine, {bx['mass_lo']:.0f}-{bx['mass_hi']:.0f} Earth masses",
                                font_size=16, color=P.TEAL)
        band_lab.next_to(plot2.p(212, 21.0), RIGHT, buff=0.1)
        depth = plot2.band(200, 1000, d["published_depth"] - d["published_depth_sigma"],
                           d["published_depth"] + d["published_depth_sigma"], P.PURPLE,
                           opacity=0.18)
        depth_line = plot2.hline(d["stacked_depth"], P.PURPLE)
        depth_lab = layout.label(f"one sector: I = {d['stacked_depth']:.1f}", font_size=16,
                                 color=P.PURPLE)
        depth_lab.next_to(plot2.p(1000, d["stacked_depth"] - 0.5), UP + LEFT, buff=0.05)
        self.play(FadeIn(plot2))
        cap5 = layout.caption("Sunlight reflected from the planet fades as distance to the 4th power",
                              font_size=22)
        self.play(FadeIn(band), FadeIn(band_lab), FadeIn(cap5), run_time=1.2)
        timing.hold_to_read(self, cap5, settle=0.2)
        cap6 = layout.caption("It stays brighter than the stack's limit out to about 700 AU",
                              font_size=22)
        self.play(FadeIn(depth), Create(depth_line), FadeIn(depth_lab), FadeOut(cap5),
                  FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=0.3)

        lim = d["motion_limit_au"]
        wall = plot2.vline(lim, P.RED)
        wall_lab = layout.label(f"drift < {d['min_displacement_px']:.0f} pixels per sector",
                                font_size=15, color=P.RED)
        wall_lab.next_to(plot2.p(lim, 16.0), DOWN + LEFT, buff=0.1)
        cap7 = layout.caption(
            f"Beyond {lim:.0f} AU (paper: {d['published_motion_limit_au']:.0f}) it barely moves "
            f"in 27 days: no longer told apart from a star", font_size=22)
        self.play(Create(wall), FadeIn(wall_lab), FadeOut(cap6), FadeIn(cap7))
        timing.hold_to_read(self, cap7, settle=0.3)

        pbox = plot2.band(bx["distance_lo"], bx["distance_hi"], 24.0, 16.0, P.TEAL, opacity=0.07)
        frac = _readout(
            f"of the {bx['distance_lo']:.0f}-{bx['distance_hi']:.0f} AU box detectable",
            f"{100 * d['detectable_fraction']:.0f}%",
            f"{100 * d['detectable_fraction_shallow']:.0f}-"
            f"{100 * d['detectable_fraction_deep']:.0f}% for ±0.5 mag", P.TEAL)
        frac.move_to(plot2.p(505, 17.3))
        cap8 = layout.caption("A forecast, not a search: TESS could reach most of the box",
                              font_size=22)
        self.play(FadeIn(pbox), FadeIn(frac), FadeOut(cap7), FadeIn(cap8))
        timing.hold_to_read(self, cap8, frac, settle=1.0)
        self.play(FadeOut(cap8))

        layout.show_takeaway(
            self, "Stacked TESS sectors reach I = 22: deep enough for Planet Nine to 700 AU.")
