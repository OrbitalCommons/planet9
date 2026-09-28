"""Fienga et al. (2016) -- where Cassini ranging lets Planet Nine be.

The planet is placed at each true anomaly of the Batygin & Brown orbit, the
ephemeris is refitted, and the Earth-Saturn range residuals measured by Cassini
either grow (the planet cannot be there) or do not. The paper does this with
the full INPOP fit; the reproduction crate uses an analytic stand-in for the
fit (the tidal range signal along Saturn's Cassini-era arc with the part an
orbit fit absorbs removed).

Everything drawn comes from anim.json -> papers -> p9-2016-cassini-ranging:
the reproduced post-fit curve, the residual floor, the reproduced favoured
position, and the paper's intervals as labelled comparisons.
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
    ReplacementTransform,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2016-cassini-ranging"


class Plot(VGroup):
    """Axes anchored at their lower-left corner, with optional log scales and
    hand-placed tick labels. ``c2p`` takes data values."""

    def __init__(self, x, y, x_ticks, y_ticks, x_label, y_label, size=(10.0, 4.2),
                 centre=(0.0, 0.2), xlog=False, ylog=False, tick_size=15):
        super().__init__()
        self.xlog, self.ylog = xlog, ylog
        self.x, self.y = x, y
        self.w, self.h = size
        self.corner = np.array([centre[0] - self.w / 2, centre[1] - self.h / 2, 0.0])
        self.add(Line(self.c2p(x[0], y[0]), self.c2p(x[1], y[0]), color=P.MUTED, stroke_width=2),
                 Line(self.c2p(x[0], y[0]), self.c2p(x[0], y[1]), color=P.MUTED, stroke_width=2))
        xt, yt = VGroup(), VGroup()
        for v, text in x_ticks.items():
            p = self.c2p(v, y[0])
            self.add(Line(p, p + DOWN * 0.08, color=P.MUTED, stroke_width=2))
            xt.add(layout.label(text, font_size=tick_size, color=P.MUTED).next_to(p, DOWN, buff=0.14))
        for v, text in y_ticks.items():
            p = self.c2p(x[0], v)
            self.add(Line(p, p + LEFT * 0.08, color=P.MUTED, stroke_width=2))
            yt.add(layout.label(text, font_size=tick_size, color=P.MUTED).next_to(p, LEFT, buff=0.14))
        self.add(xt, yt)
        self.add(layout.label(x_label, font_size=18).next_to(xt, DOWN, buff=0.14)
                 .set_x(self.corner[0] + self.w / 2))
        self.add(layout.label(y_label, font_size=15).rotate(np.pi / 2).next_to(yt, LEFT, buff=0.14)
                 .set_y(self.corner[1] + self.h / 2))

    @staticmethod
    def _f(v, log):
        return np.log10(v) if log else v

    def c2p(self, x, y):
        fx = (self._f(x, self.xlog) - self._f(self.x[0], self.xlog)) / (
            self._f(self.x[1], self.xlog) - self._f(self.x[0], self.xlog))
        fy = (self._f(y, self.ylog) - self._f(self.y[0], self.ylog)) / (
            self._f(self.y[1], self.ylog) - self._f(self.y[0], self.ylog))
        return self.corner + np.array([fx * self.w, fy * self.h, 0.0])

    def band(self, a, b, color, opacity):
        p0, p1 = self.c2p(a, self.y[0]), self.c2p(b, self.y[1])
        r = Polygon(p0, [p1[0], p0[1], 0], p1, [p0[0], p1[1], 0], stroke_width=0)
        return r.set_fill(color, opacity=opacity)


def signed(nu):
    """True anomaly on the paper's (-180, 180] convention."""
    return nu - 360.0 if nu > 180.0 else nu


def signed_intervals(intervals):
    """[0, 360) intervals re-expressed on (-180, 180], merged across the seam."""
    out = sorted((signed(a + 1e-9), signed(b)) for a, b in intervals)
    merged = []
    for a, b in out:
        if merged and abs(a - merged[-1][1]) < 1e-6:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    return merged


def orbit_arc(a, e, f0, f1, color, stroke_width=5.0, n=80):
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners([orbits.orbit_point(a, e, np.radians(f))
                             for f in np.linspace(f0, f1, n)])
    return m


class CassiniRanging2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        orbit, pub = d["orbit"], d["published"]
        floor = d["floor_m"]

        nu = np.array([signed(v) for v in d["curve"]["nu_deg"]])
        post = np.array(d["curve"]["postfit_m"])
        pre = np.array(d["curve"]["prefit_m"])
        order = np.argsort(nu)
        nu, post, pre = nu[order], post[order], pre[order]
        nu = np.concatenate([[-180.0], nu])
        post = np.concatenate([[post[-1]], post])
        pre = np.concatenate([[pre[-1]], pre])

        forbidden = signed_intervals(pub["excluded_deg"])
        noticed = signed_intervals(d["excluded_deg"])
        fav_lo, fav_hi = pub["favored_interval_deg"]

        self.add(paper.scene_header(CRATE))

        # 1. the raw tug on Saturn's range, then what survives the orbit refit
        ax = Plot(
            (-180, 180), (3, 3000),
            {v: f"{v}°" for v in range(-180, 181, 60)},
            {10: "10", 30: "30", 100: "100", 300: "300", 1000: "1000"},
            "true anomaly of Planet Nine  (0° = perihelion)",
            "Earth-Saturn range signal  (m)",
            size=(10.6, 4.0), centre=(0.3, 0.35), ylog=True)
        curve = widgets.curve(ax, nu, post, color=P.ORANGE, stroke_width=3.5)
        raw = widgets.curve(ax, nu, pre, color=P.ORANGE, stroke_width=3.5)
        raw_lab = layout.label("raw tug over the Cassini years", font_size=16, color=P.ORANGE)
        raw_lab.next_to(ax.c2p(0, pre.max()), UP, buff=0.08)
        cap = layout.caption(
            f"Reproduced: a {orbit['mass_earth']:.0f} M⊕ planet at each point of its "
            f"{orbit['a_au']:.0f} AU orbit tugs Saturn's range", font_size=22)
        self.play(FadeIn(ax), run_time=0.9)
        self.play(Create(raw), FadeIn(raw_lab), FadeIn(cap), run_time=1.8)
        timing.hold_to_read(self, cap, raw_lab, settle=0.6)

        ghost = raw.copy().set_stroke(P.MUTED, width=1.6, opacity=0.6)
        post_lab = layout.label("left after refitting Saturn's orbit", font_size=16,
                                color=P.ORANGE)
        post_lab.next_to(ax.c2p(0, post.max()), UP, buff=0.08)
        capb = layout.caption(
            "Refitting Saturn's orbit soaks up most of it; only the leftover can be seen",
            font_size=22)
        self.add(ghost)
        self.play(ReplacementTransform(raw, curve), ReplacementTransform(raw_lab, post_lab),
                  FadeOut(cap), FadeIn(capb), run_time=1.8)
        timing.hold_to_read(self, capb, settle=0.8)
        self.play(FadeOut(post_lab), run_time=0.4)
        cap = capb

        # 2. the residual floor decides what would have been noticed
        fl = DashedLine(ax.c2p(-180, floor), ax.c2p(180, floor), color=P.PURPLE,
                        stroke_width=2.5)
        fl_lab = layout.label(f"residual floor  {floor:.0f} m", font_size=14, color=P.PURPLE)
        fl_lab.next_to(ax.c2p(-180, floor), UP, buff=0.08, aligned_edge=LEFT).shift(RIGHT * 0.1)
        span = "  and  ".join(f"{a:.0f}° to {b:.0f}°" for a, b in noticed)
        cap2 = layout.caption(
            f"Reproduced: from {span} the signal clears Cassini's residual floor",
            font_size=22)
        self.play(Create(fl), FadeIn(fl_lab), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        timing.hold_to_read(self, cap2, fl_lab, settle=0.8)

        # 3. the paper's verdict from the full ephemeris fit
        bands = VGroup(*[ax.band(a, b, P.RED, 0.2) for a, b in forbidden])
        green = ax.band(fav_lo, fav_hi, P.GREEN, 0.3)
        wide = max(forbidden, key=lambda iv: iv[1] - iv[0])
        forb_lab = layout.label("paper: forbidden", font_size=15, color=P.RED)
        forb_lab.next_to(ax.c2p(0.5 * (wide[0] + wide[1]), ax.y[1]), UP, buff=0.1)
        fav_lab = layout.label(f"paper: favoured  {pub['favored_deg']:.1f}°", font_size=15,
                               color=P.GREEN)
        fav_lab.next_to(ax.c2p(0.5 * (fav_lo + fav_hi), ax.y[1]), UP, buff=0.1)
        mine = Dot(ax.c2p(d["favored_deg"], d["favored_postfit_m"]), radius=0.09,
                   color=P.ORANGE).set_z_index(4)
        mine_lab = layout.label(f"reproduced:\n{d['favored_deg']:.1f}°", font_size=15,
                                color=P.ORANGE, line_spacing=0.9)
        mine_lab.next_to(mine, UP, buff=0.12).shift(RIGHT * 0.45)
        cap3 = layout.caption(
            "The paper's full refit forbids the near side and improves the fit at one spot",
            font_size=22)
        self.play(FadeIn(bands), FadeIn(forb_lab), FadeOut(cap2), FadeIn(cap3), run_time=1.0)
        self.play(FadeIn(green), FadeIn(fav_lab), run_time=0.8)
        self.play(FadeIn(mine), FadeIn(mine_lab), run_time=0.6)
        timing.hold_to_read(self, cap3, forb_lab, fav_lab, mine_lab, settle=1.2)
        self.play(FadeOut(VGroup(ax, ghost, curve, fl, fl_lab, bands, green, forb_lab, fav_lab,
                                 mine, mine_lab, cap3)))

        # 4. the same verdict on the orbit
        scale = 5.6 / (2.0 * orbit["a_au"])
        a, e = orbit["a_au"] * scale, orbit["e"]
        shift = np.array([2.3, 0.2, 0.0])
        ghost = orbits.ellipse_orbit(a, e, color=P.MUTED, stroke_width=1.6, opacity=0.7)
        arcs = VGroup(*[orbit_arc(a, e, f0, f1, P.RED) for f0, f1 in forbidden])
        fav = orbit_arc(a, e, fav_lo, fav_hi, P.GREEN, stroke_width=7.0)
        sun = orbits.sun(radius=0.11)
        here = Dot(orbits.orbit_point(a, e, np.radians(d["favored_deg"])), radius=0.08,
                   color=P.ORANGE).set_z_index(4)
        VGroup(ghost, arcs, fav, sun, here).shift(shift)
        sun_l = layout.label("Sun and Saturn", font_size=13, color=P.MUTED)
        sun_l.next_to(sun, DOWN, buff=0.12)
        notes = VGroup(
            layout.label("paper: forbidden", font_size=18, color=P.RED),
            layout.label(
                f"paper: favoured, {fav_lo:.0f}° to {fav_hi:.0f}°,\n"
                f"about {pub['favored_distance_au']:.0f} AU from the Sun",
                font_size=18, color=P.GREEN, line_spacing=0.9),
            layout.label(
                f"reproduced: {d['favored_deg']:.1f}°,\n"
                f"{d['favored_distance_au']:.0f} AU from the Sun",
                font_size=18, color=P.ORANGE, line_spacing=0.9),
            layout.label("grey: not excluded", font_size=18, color=P.MUTED),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        notes.move_to([0, 1.0, 0]).to_edge(LEFT, buff=0.6)
        self.play(FadeIn(sun), FadeIn(sun_l), Create(ghost), run_time=1.0)
        self.play(Create(arcs), FadeIn(notes[0]), run_time=1.0)
        self.play(Create(fav), FadeIn(notes[1]), run_time=0.8)
        self.play(FadeIn(here), FadeIn(notes[2]), FadeIn(notes[3]), run_time=0.8)
        timing.hold_to_read(self, notes, settle=0.8)

        layout.show_takeaway(
            self, "Cassini forbids the near side of the orbit and favours a spot "
                  f"~{pub['favored_distance_au']:.0f} AU out.")
