"""Finale -- how much of the prediction the searches have erased, where the
rest hides, and how much of it Rubin will settle.

CarvingTheSky puts one Brown & Batygin (2021) reference population on tonight's
sky and erases it survey by survey (ZTF, DES, Pan-STARRS1). WhereItHides shows
the survivors in distance and brightness and why they are far and faint.
WhatRubinTakes runs the Rubin/LSST baseline over the survivors and shows what it
leaves. Every number is read from ``anim.json -> preface -> finale``
(``crates/p9-anim-data/src/preface/finale.rs``); the space-telescope scenes that
follow are in space_strategy.py.
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
    Dot,
    FadeIn,
    FadeOut,
    GrowFromEdge,
    Line,
    Rectangle,
    Scene,
    Transform,
    UpdateFromAlphaFunc,
    VGroup,
    VMobject,
    ManimColor,
    interpolate_color,
)

import p9_manim as P
from p9_manim import dataio, layout, ledger, orbits, sky, timing

DOT_OPACITY = 0.85
# draw every SHOW-th orbit (all of them enter every quoted fraction)
SHOW = 2
SHORT = {"Pan-STARRS1": "PS1"}


def _finale():
    return dataio.section("preface")["finale"]


def _header(text):
    badge = layout.concept_badge("FINALE")
    title = layout.label(text, font_size=24, color=P.FG, weight="BOLD")
    title.to_edge(UP, buff=0.32)
    return VGroup(badge, title)


def _survival(orbits_, keys):
    s = np.ones(len(orbits_["ra"]))
    for k in keys:
        s *= 1.0 - np.asarray(orbits_[k])
    return s


def _carve(dots, before, after, flash=P.RED, run_time=2.6):
    """Flash each dot toward ``flash`` in proportion to the share of it this
    step removes, then dim it to its new surviving weight."""
    lost = np.where(before > 0, (before - after) / np.maximum(before, 1e-12), 0.0)
    teal, hot = ManimColor(P.TEAL), ManimColor(flash)

    def update(group, alpha):
        red = min(alpha / 0.45, 1.0)
        fade = max((alpha - 0.45) / 0.55, 0.0)
        for d, f, b, a in zip(group, lost, before, after):
            peak = interpolate_color(teal, hot, f * red)
            d.set_fill(interpolate_color(peak, teal, fade),
                       opacity=DOT_OPACITY * (b + (a - b) * fade))

    return UpdateFromAlphaFunc(dots, update, run_time=run_time)


class _Gauge(VGroup):
    """A horizontal 'share of the prediction' bar filled left to right."""

    def __init__(self, title, width=7.6, height=0.26, centre=(0.9, 2.62, 0)):
        super().__init__()
        self.w, self.h = width, height
        self.track = Rectangle(width=width, height=height, color=P.MUTED, stroke_width=1.2)
        self.track.set_fill(P.TEAL, opacity=0.12).move_to(centre)
        self.title = layout.label(title, font_size=15, color=P.FG)
        self.title.next_to(self.track, LEFT, buff=0.2)
        self.add(self.track, self.title)
        self.filled = 0.0
        self.readout = None

    def x(self, frac):
        return self.track.get_left()[0] + self.w * frac

    def segment(self, frac, color, opacity=0.75):
        seg = Rectangle(width=max(self.w * frac, 1e-3), height=self.h, stroke_width=0)
        seg.set_fill(color, opacity=opacity)
        seg.move_to([self.x(self.filled) + self.w * frac / 2, self.track.get_center()[1], 0])
        self.filled += frac
        return seg

    def seg_label(self, seg, text, color):
        lab = layout.label(text, font_size=12, color=color)
        if lab.width + 0.1 < seg.width:
            return lab.move_to(seg)
        return lab.next_to(seg, DOWN, buff=0.06)

    def value(self, text, color=P.FG):
        lab = layout.label(text, font_size=17, color=color, weight="BOLD")
        return lab.next_to(self.track, RIGHT, buff=0.2)


def _survivor_map(d, width=12.0):
    m = sky.SkyMap(width=width, dec_range=(-60, 60), centre=(0.0, 0.05, 0.0))
    ecl, gal = m.reference_curves()
    s = dataio.section("sky")
    e_ra, e_dec = min(s["ecliptic"], key=lambda p: abs(p[0] - 20.0))
    g_ra, g_dec = min(s["galactic_plane"], key=lambda p: abs(p[0] - 330.0) + abs(p[1] - 40.0))
    tags = VGroup(
        layout.label("ecliptic", font_size=12, color=P.ORANGE).next_to(m.p(e_ra, e_dec), DOWN,
                                                                       buff=0.12),
        layout.label("Milky Way", font_size=12, color=P.PURPLE).next_to(m.p(g_ra, g_dec), RIGHT,
                                                                        buff=0.12),
    )
    o = d["orbits"]
    dots = VGroup(*[Dot(m.p(ra, dec), radius=0.026, color=P.TEAL).set_opacity(DOT_OPACITY)
                    for ra, dec in zip(o["ra"][::SHOW], o["dec"][::SHOW])])
    return m, VGroup(ecl, gal, tags), dots


class CarvingTheSky(Scene):
    """The prediction on tonight's sky, erased survey by survey."""

    def construct(self):
        d = _finale()
        o = d["orbits"]
        self.add(_header("What the searches have ruled out"))

        m, refs, dots = _survivor_map(d)
        self.play(FadeIn(m), Create(refs), run_time=1.0)
        cap = layout.caption(
            "Each dot: one orbit the 2021 fit allows, "
            "placed where the planet would be tonight.", font_size=20)
        self.play(FadeIn(dots, lag_ratio=0.0005), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=1.0)

        gauge = _Gauge("ruled out")
        value = gauge.value("0%")
        self.play(FadeIn(gauge), FadeIn(value), FadeOut(cap), run_time=0.6)

        keys, before = [], np.ones(len(dots))
        for sv in d["surveys"]:
            keys.append(sv["key"])
            after = _survival(o, keys)[::SHOW]
            if sv["key"] == "p_des":
                foot = VGroup(*[m.box(ra - 1, ra + 1, dec - 1, dec + 1,
                                      color=P.PURPLE,
                                      opacity=0.16, stroke_width=0)
                                for ra, dec in d["des_cells"] if dec > -61])
                text = (f"DES ({sv['year']}): {sv['area_deg2']:.0f} deg² of southern sky "
                        f"to r {sv['depth_r']:.1f}. Deep, but narrow.")
            else:
                foot = m.dec_band(sv["dec_min"], 90, color=P.PURPLE, opacity=0.08,
                                  stroke_width=1.2)
                text = (f"{sv['name']} ({sv['year']}): everything north of "
                        f"{sv['dec_min']:.0f}°, down to V {sv['depth_v']:.1f}.")
            cap = layout.caption(text, font_size=20)
            share = sv["cumulative"] - gauge.filled
            seg = gauge.segment(share, P.RED, opacity=0.35 + 0.2 * len(keys))
            seg_lab = gauge.seg_label(seg, SHORT.get(sv["name"], sv["name"]), P.FG)
            new_value = gauge.value(f"{100 * sv['cumulative']:.0f}%", P.RED)
            self.play(FadeIn(foot), FadeIn(cap), run_time=0.8)
            self.play(_carve(dots, before, after), GrowFromEdge(seg, LEFT),
                      Transform(value, new_value), run_time=2.6)
            self.play(FadeIn(seg_lab), run_time=0.3)
            timing.hold_to_read(self, cap, settle=0.8)
            self.play(FadeOut(foot), FadeOut(cap), run_time=0.5)
            before = after

        published = ledger.entry("p9-2024-panstarrs")["result"]["published"]
        tick = Line(UP * 0.24, DOWN * 0.24, color=P.FG, stroke_width=2.5)
        tick.move_to([gauge.x(published), gauge.track.get_center()[1], 0])
        tick_lab = layout.label(f"paper: {100 * published:.0f}%", font_size=12, color=P.FG)
        tick_lab.next_to(tick, UP, buff=0.05)
        cap = layout.caption(
            f"Our reproduction: {100 * d['surveys'][-1]['cumulative']:.0f}% ruled out; "
            f"the Pan-STARRS1 paper reports {100 * published:.0f}%.", font_size=20)
        self.play(Create(tick), FadeIn(tick_lab), FadeIn(cap))
        timing.hold_to_read(self, cap, settle=0.8)

        cap2 = layout.caption(
            f"The {100 * d['hiding']['remaining']:.0f}% that survives is still spread "
            "around the sky: position alone no longer narrows it.", font_size=20)
        self.play(FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.8)
        self.play(FadeOut(cap2))

        layout.show_takeaway(
            self, f"ZTF, DES and Pan-STARRS1 have erased "
                  f"{100 * d['surveys'][-1]['cumulative']:.0f}% of the predicted orbits.")


class _DistV(VGroup):
    """Distance (AU) across, magnitude down (bright at the top), drawn by hand
    so the axes sit on the plot edges rather than crossing at zero."""

    D0, D1, V0, V1 = 100.0, 1100.0, 16.0, 25.0

    def __init__(self, left=-5.9, right=2.3, bottom=-2.05, top=2.45):
        super().__init__()
        self.box = (left, right, bottom, top)
        frame = VGroup(Line([left, bottom, 0], [right, bottom, 0]),
                       Line([left, bottom, 0], [left, top, 0])).set_stroke(P.MUTED, 1.6)
        ticks = VGroup()
        for dist in range(100, 1101, 200):
            p = self.c2p(dist, self.V1)
            ticks.add(Line(p, p + DOWN * 0.08).set_stroke(P.MUTED, 1.6),
                      layout.label(f"{dist}", font_size=14, color=P.FG)
                      .next_to(p, DOWN, buff=0.14))
        for v in range(16, 26, 2):
            p = self.c2p(self.D0, v)
            ticks.add(Line(p, p + LEFT * 0.08).set_stroke(P.MUTED, 1.6),
                      layout.label(f"{v}", font_size=14, color=P.FG)
                      .next_to(p, LEFT, buff=0.14))
        xlab = layout.label("distance from the Sun today (AU)", font_size=16, color=P.FG)
        xlab.next_to(ticks, DOWN, buff=0.12).set_x((left + right) / 2)
        ylab = layout.label("brightness V  (fainter ↓)", font_size=16, color=P.FG)
        ylab.rotate(np.pi / 2).next_to(ticks, LEFT, buff=0.15).set_y((bottom + top) / 2)
        self.add(frame, ticks, xlab, ylab)

    def c2p(self, dist, v):
        left, right, bottom, top = self.box
        x = left + (dist - self.D0) / (self.D1 - self.D0) * (right - left)
        y = top - (v - self.V0) / (self.V1 - self.V0) * (top - bottom)
        return np.array([x, y, 0.0])


class WhereItHides(Scene):
    """The survivors in distance and brightness, and why they are far and faint."""

    def construct(self):
        d = _finale()
        o, h = d["orbits"], d["hiding"]
        self.add(_header("Where it could still hide: far and faint"))

        # 1. the prediction in distance x brightness, carved by the searches
        ax = _DistV()
        self.play(FadeIn(ax), run_time=1.0)
        dist, v = np.asarray(o["dist"])[::SHOW], np.asarray(o["v"])[::SHOW]
        keep = (dist >= 100) & (dist <= 1100) & (v >= 16) & (v <= 25)
        dots = VGroup(*[Dot(ax.c2p(x, y), radius=0.022, color=P.TEAL).set_opacity(DOT_OPACITY)
                        for x, y in zip(dist[keep], v[keep])])
        cap = layout.caption("The same orbits by distance and brightness: farther means fainter.",
                             font_size=20)
        self.play(FadeIn(dots, lag_ratio=0.0005), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.6)

        limits = VGroup()
        for sv in d["surveys"]:
            depth = sv.get("depth_v")
            if depth is None:
                continue
            ln = DashedLine(ax.c2p(100, depth), ax.c2p(1100, depth), color=P.RED,
                            stroke_width=2, dash_length=0.1)
            lab = layout.label(f"{sv['name']} limit V {depth:.1f}", font_size=12, color=P.RED)
            lab.next_to(ax.c2p(1100, depth), LEFT, buff=0.08).shift(UP * 0.14)
            limits.add(VGroup(ln, lab))
        survive = _survival(o, [s["key"] for s in d["surveys"]])[::SHOW][keep]
        cap2 = layout.caption(
            "Each search erased what was brighter than its limit, inside its own sky.",
            font_size=20)
        self.play(FadeOut(cap), FadeIn(cap2), Create(limits), run_time=0.9)
        self.play(_carve(dots, np.ones(len(dots)), survive), run_time=2.6)
        timing.hold_to_read(self, cap2, settle=0.8)

        # 2. familiar anchors for the typical survivor
        vm, dm = h["v_median_survivors"], h["dist_median_survivors"]
        vline = DashedLine(ax.c2p(dm, 25), ax.c2p(dm, 16), color=P.TEAL, stroke_width=2,
                           dash_length=0.08)
        hline = DashedLine(ax.c2p(100, vm), ax.c2p(1100, vm), color=P.TEAL, stroke_width=2,
                           dash_length=0.08)
        panel = VGroup(
            layout.label("typical survivor", font_size=16, color=P.TEAL, weight="BOLD"),
            layout.label(f"V {vm:.1f}", font_size=26, color=P.FG, weight="BOLD"),
            layout.label(f"{h['fainter_than_pluto']:,.0f}× fainter than Pluto",
                         font_size=14, color=P.FG),
            layout.label(f"(Pluto: V {h['pluto_v']:.1f})", font_size=12, color=P.MUTED),
            layout.label(f"{dm:.0f} AU", font_size=26, color=P.FG, weight="BOLD"),
            layout.label(f"{h['times_neptune']:.0f}× Neptune's distance",
                         font_size=14, color=P.FG),
            layout.label(f"(Neptune: {h['neptune_au']:.0f} AU)", font_size=12, color=P.MUTED),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        panel[4].shift(DOWN * 0.25)
        panel[5:].shift(DOWN * 0.25)
        panel.move_to([4.95, 0.3, 0])
        cap3 = layout.caption("Half of the surviving probability is fainter, and farther, "
                              "than these lines.", font_size=20)
        self.play(Create(vline), Create(hline), FadeIn(panel, lag_ratio=0.1), FadeOut(cap2),
                  FadeIn(cap3), run_time=1.4)
        timing.hold_to_read(self, cap3, panel, settle=1.0)
        self.play(*[FadeOut(x) for x in [ax, dots, limits, vline, hline, panel, cap3]],
                  run_time=0.7)

        # 3. why far: Kepler's clock on the 2021 best-fit orbit
        c = d["orbit_clock"]
        scale = 5.2 / (c["q"] + c["big_q"])
        sun_x = -0.6
        pts = np.array(c["points"]) * scale + np.array([sun_x, 0.2])
        ell = VMobject(color=P.BLUE, stroke_width=2.2)
        ell.set_points_smoothly([[x, y, 0] for x, y in np.vstack([pts, pts[:1]])])
        sun = orbits.sun(radius=0.12).move_to([sun_x, 0.2, 0])
        ring = Circle(radius=c["a"] * scale, color=P.MUTED, stroke_width=1.4)
        ring.move_to(sun).set_stroke(opacity=0.8)
        ring_lab = layout.label(f"mean distance a = {c['a']:.0f} AU", font_size=13, color=P.MUTED)
        ring_lab.next_to(ring, UP, buff=0.08)
        nep = Circle(radius=h["neptune_au"] * scale, color=P.FG, stroke_width=1.4).move_to(sun)
        nep_lab = layout.label("Neptune", font_size=12, color=P.FG)
        nep_lab.next_to(nep, DOWN, buff=0.14)
        r = np.hypot(*(np.array(c["points"]).T))
        ticks = VGroup(*[Dot([x, y, 0], radius=0.07,
                             color=P.TEAL if rr > c["a"] else P.MUTED)
                         for (x, y), rr in zip(pts, r)])
        peri = layout.label(f"perihelion {c['q']:.0f} AU\nV {c['v_peri']:.1f}", font_size=13,
                            color=P.FG).next_to(ring, RIGHT, buff=0.15)
        aph = layout.label(f"aphelion {c['big_q']:.0f} AU\nV {c['v_aph']:.1f}", font_size=13,
                           color=P.FG).next_to([*pts[len(pts) // 2], 0], LEFT, buff=0.2)
        cap4 = layout.caption(
            f"The 2021 best-fit orbit, one dot every {c['period_yr'] / len(pts):.0f} years: "
            "it crawls where it is far.", font_size=20)
        self.play(Create(ell), FadeIn(sun), Create(ring), FadeIn(ring_lab), Create(nep),
                  FadeIn(nep_lab), FadeIn(cap4), run_time=1.2)
        self.play(FadeIn(ticks, lag_ratio=0.15), FadeIn(peri), FadeIn(aph), run_time=2.0)
        timing.hold_to_read(self, cap4, peri, aph, settle=0.8)

        bars = VGroup()
        rows = [("every predicted orbit", h["beyond_a_share_prior"], P.MUTED),
                ("the survivors", h["beyond_a_share_survivors"], P.TEAL)]
        head = layout.label("share beyond its mean distance", font_size=14, color=P.FG)
        for name, frac, col in rows:
            track = Rectangle(width=2.6, height=0.24, color=P.MUTED, stroke_width=1)
            fill = Rectangle(width=2.6 * frac, height=0.24, stroke_width=0)
            fill.set_fill(col, opacity=0.8).align_to(track, LEFT)
            val = layout.label(f"{100 * frac:.0f}%", font_size=14, color=col, weight="BOLD")
            val.next_to(track, RIGHT, buff=0.1)
            lab = layout.label(name, font_size=13, color=P.FG).next_to(track, UP, buff=0.06)
            lab.align_to(track, LEFT)
            bars.add(VGroup(lab, track, fill, val))
        bars.arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        box = VGroup(head, bars).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        box.move_to([4.9, -1.6, 0])
        cap5 = layout.caption(
            f"It spends {100 * c['time_beyond_a']:.0f}% of each orbit beyond a, "
            "and the searches caught the near side first.", font_size=20)
        self.play(FadeOut(cap4), FadeIn(cap5), FadeIn(box, lag_ratio=0.1), run_time=1.0)
        timing.hold_to_read(self, cap5, box, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, f"What survives is V ≈ {vm:.0f} and ~{round(dm, -1):.0f} AU away, near aphelion.")


class WhatRubinTakes(Scene):
    """The Rubin/LSST baseline run over the survivors, and what it leaves."""

    def construct(self):
        d = _finale()
        o, rb = d["orbits"], d["rubin"]
        self.add(_header("What Rubin will settle, and what it leaves"))

        m, refs, dots = _survivor_map(d)
        survive = _survival(o, [s["key"] for s in d["surveys"]])[::SHOW]
        for dot, s in zip(dots, survive):
            dot.set_opacity(DOT_OPACITY * s)
        gauge = _Gauge("the prediction")
        excluded = gauge.segment(d["surveys"][-1]["cumulative"], P.RED, opacity=0.8)
        ex_lab = gauge.seg_label(excluded, "ruled out", P.FG)
        cap = layout.caption(
            f"Start from the survivors: {100 * d['hiding']['remaining']:.0f}% "
            "of the prediction.", font_size=20)
        self.play(FadeIn(m), Create(refs), FadeIn(dots), FadeIn(gauge), FadeIn(excluded),
                  FadeIn(ex_lab), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.6)

        foot = m.dec_band(rb["dec_min"], rb["dec_max"], color=P.PURPLE, opacity=0.12,
                          stroke_width=1.6)
        edge = layout.label(f"Rubin's northern limit, {rb['dec_max']:+.0f}°", font_size=13,
                            color=P.PURPLE)
        edge.add_background_rectangle(color=P.BG, opacity=0.85, buff=0.04)
        edge.next_to(m.p(235, rb["dec_max"]), UP, buff=0.06)
        cap2 = layout.caption(
            f"Rubin/LSST: ten years south of {rb['dec_max']:+.0f}°, to r {rb['depth_r']:.1f}, "
            f"skipping the Milky Way (|b| < {rb['gal_b_min']:.0f}°).", font_size=20)
        self.play(FadeIn(foot), FadeIn(edge), FadeOut(cap), FadeIn(cap2), run_time=0.9)
        timing.hold_to_read(self, cap2, settle=0.6)

        after = survive * (1.0 - np.asarray(o["p_rubin"])[::SHOW])
        seg = gauge.segment(rb["of_prediction"], P.PURPLE, opacity=0.8)
        seg_lab = gauge.seg_label(seg, "Rubin", P.FG)
        left = gauge.value(f"{100 * rb['left_after']:.0f}% left", P.TEAL)
        cap3 = layout.caption(
            f"It should find the planet on {100 * rb['of_survivors']:.0f}% of the surviving "
            "orbits.", font_size=20)
        self.play(FadeOut(cap2), FadeIn(cap3), run_time=0.4)
        self.play(_carve(dots, survive, after, flash=P.PURPLE), GrowFromEdge(seg, LEFT),
                  FadeIn(seg_lab), FadeIn(left), run_time=2.8)
        timing.hold_to_read(self, cap3, settle=0.8)

        cap4 = layout.caption(
            f"Of what it misses, {100 * rb['left_north_share']:.0f}% is north of "
            f"{rb['dec_max']:+.0f}° and {100 * rb['left_plane_share']:.0f}% is in the "
            "Milky Way's plane.", font_size=20)
        self.play(FadeOut(foot), FadeOut(cap3), FadeIn(cap4), run_time=0.8)
        timing.hold_to_read(self, cap4, settle=1.2)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, "Rubin can settle half of what is left; the rest hides north and in the plane.")
