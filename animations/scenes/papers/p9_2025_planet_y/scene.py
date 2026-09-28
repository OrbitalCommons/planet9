"""Siraj, Chyba & Tremaine (2025) -- the mean plane of the distant Kuiper belt.

Each distant orbit's plane wobbles (its pole precesses) around a local
"forced" plane, so the average plane of many orbits is that forced plane. With
only the known planets it is the invariable plane at every distance. The paper
measures the mean plane bin by bin and finds it tilted by about 15 degrees at
80-200 AU but flat at 50-80 AU, and shows a Mercury-to-Earth-mass planet at
100-200 AU could hold such a warp. The scene draws the crate's forced-plane
profiles and its (mass, distance) maps of the tilt each bin would carry.
Data: anim.json -> papers -> p9-2025-planet-y.
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
    Line,
    Polygon,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-planet-y"


def rect(a, b, color, opacity):
    return Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0).set_fill(
        color, opacity=opacity)


class TiltMap(VGroup):
    """A (log mass) x (semi-major axis) grid of cells shaded by tilt."""

    def __init__(self, masses, axes_au, tilts, vmax, width=4.6, height=3.3, title=""):
        super().__init__()
        self.lm = np.log10(masses)
        self.a = np.asarray(axes_au)
        self.w, self.h = width, height
        dl = self.lm[1] - self.lm[0]
        da = self.a[1] - self.a[0]
        self.l0, self.l1 = self.lm[0] - dl / 2, self.lm[-1] + dl / 2
        self.a0, self.a1 = self.a[0] - da / 2, self.a[-1] + da / 2
        vals = np.asarray(tilts).reshape(len(self.a), len(self.lm))
        frame = Rectangle(width=width, height=height, color=P.MUTED, stroke_width=1.2)
        frame.move_to([width / 2, height / 2, 0])
        self.frame = frame
        cells = VGroup()
        for iy, a in enumerate(self.a):
            for ix, lm in enumerate(self.lm):
                f = min(vals[iy, ix] / vmax, 1.0)
                cells.add(rect(self.p(10 ** (lm - dl / 2), a - da / 2),
                               self.p(10 ** (lm + dl / 2), a + da / 2), P.ORANGE, 0.05 + 0.85 * f))
        self.add(cells, frame)
        self.add(layout.label(title, font_size=16, weight="BOLD").next_to(frame, UP, buff=0.12))

    def p(self, m, a):
        x = (np.log10(m) - self.l0) / (self.l1 - self.l0) * self.w
        y = (a - self.a0) / (self.a1 - self.a0) * self.h
        return np.array([x, y, 0.0]) + self.frame.get_corner(DOWN + LEFT)


class PlanetY2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        earth = d["profiles"][-1]

        self.add(paper.scene_header(CRATE))

        # 1. mechanism: orbit poles circle a forced pole; their average is the mean plane
        c = np.array([-3.3, 0.1, 0.0])
        k = 2.5 / 20.0  # scene units per degree of tilt
        ring = VGroup(*[Circle(radius=r * k, color=P.MUTED, stroke_width=1).move_to(c)
                        .set_stroke(opacity=0.5) for r in (10, 20)])
        ring_lab = VGroup(*[layout.label(f"{r}°", font_size=13, color=P.MUTED)
                            .next_to(c + UP * r * k, UP + RIGHT, buff=0.02) for r in (10, 20)])
        cross = VGroup(Line(c + LEFT * 2.6, c + RIGHT * 2.6, color=P.MUTED, stroke_width=1),
                       Line(c + DOWN * 2.6, c + UP * 2.6, color=P.MUTED, stroke_width=1)
                       ).set_stroke(opacity=0.4)
        centre_lab = layout.label("pole of the planets' plane", font_size=14, color=P.MUTED)
        centre_lab.next_to(c + DOWN * 2.0, LEFT, buff=0.15)
        centre_arrow = Line(centre_lab.get_right() + RIGHT * 0.05, c + DOWN * 0.12,
                            color=P.MUTED, stroke_width=1.2)
        centre_lab = VGroup(centre_lab, centre_arrow)
        shift = ValueTracker(0.0)
        t = ValueTracker(0.0)
        free = [3.0, 5.0, 7.0, 9.0, 11.0, 6.0, 4.0]
        rates = [1.0, 0.8, 0.65, 0.55, 0.45, 0.7, 0.9]
        phase = np.linspace(0, 2 * np.pi, len(free), endpoint=False)

        def forced():
            return c + RIGHT * shift.get_value() * k

        def poles():
            g = VGroup()
            f = forced()
            for r, w, ph in zip(free, rates, phase):
                ang = ph + w * t.get_value()
                g.add(Circle(radius=r * k, color=P.GREEN, stroke_width=0.8).move_to(f)
                      .set_stroke(opacity=0.25))
                g.add(Dot(f + r * k * np.array([np.cos(ang), np.sin(ang), 0]), radius=0.05,
                          color=P.GREEN))
            g.add(Dot(f, radius=0.08, color=P.ORANGE).set_z_index(3))
            return g

        live = always_redraw(poles)
        mean_lab = always_redraw(lambda: layout.label("mean plane", font_size=14, color=P.ORANGE)
                                 .next_to(forced(), UP + RIGHT, buff=0.06))
        text = VGroup(
            layout.label("each dot: the pole of one distant orbit", font_size=16, color=P.GREEN),
            layout.label("the planets make each pole circle", font_size=16),
            layout.label("a forced pole (orange)", font_size=16),
            layout.label("average many orbits: you get", font_size=16),
            layout.label("the forced pole = the mean plane", font_size=16, color=P.ORANGE),
            layout.label("schematic: orbits illustrative", font_size=13, color=P.MUTED),
        ).arrange(DOWN, buff=0.16, aligned_edge=LEFT).move_to([3.4, 0.9, 0])
        cap = layout.caption("With only the known planets, the mean plane is their plane",
                             font_size=22)
        self.play(FadeIn(ring), FadeIn(ring_lab), FadeIn(cross), FadeIn(centre_lab), FadeIn(live),
                  FadeIn(mean_lab), FadeIn(text), FadeIn(cap), run_time=1.0)
        self.play(t.animate.set_value(5.0), run_time=4.0, rate_func=linear)
        timing.hold_to_read(self, cap, text[:5], settle=0.2)
        tilt = d["warp_bin_tilt_earth_deg"]
        cap2 = layout.caption(
            f"An Earth-mass planet at {earth['a_au']:.0f} AU moves the forced pole "
            f"{tilt:.0f}° at 80–200 AU (computed)", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), shift.animate.set_value(tilt),
                  t.animate.set_value(9.0), run_time=3.0, rate_func=linear)
        self.play(t.animate.set_value(12.0), run_time=2.2, rate_func=linear)
        timing.hold_to_read(self, cap2, settle=0.6)
        self.play(FadeOut(VGroup(ring, ring_lab, cross, centre_lab, live, mean_lab, text, cap2)),
                  run_time=0.7)

        # 2. computed: the forced tilt against distance, and the paper's bins
        ag = np.array(d["a_grid_au"])
        ax, labs = widgets.labeled_axes(
            [0, 400, 50], [0, 16, 4], x_label="semi-major axis  (AU)",
            y_label="tilt of the mean plane  (deg)", y_rotate=True, numbers=True,
            x_length=8.8, y_length=3.8, shift_down=0)
        ax.move_to([-0.9, 0.35, 0])
        labs[0].next_to(ax, DOWN, buff=0.5)
        labs[1].next_to(ax, LEFT, buff=0.5)
        bins = VGroup()
        bin_tags = VGroup()
        bin_notes = ["flat", f"warped, {100 * pub['confidence'][0]['confidence']:.0f}%",
                     "not significant"]
        for b, note in zip(d["bins"], bin_notes):
            col = P.GREEN if b["warped"] else P.MUTED
            bins.add(rect(ax.c2p(b["lo_au"], 0), ax.c2p(b["hi_au"], 16), col,
                          0.14 if b["warped"] else 0.08))
            bin_tags.add(layout.label(note, font_size=14, color=col)
                         .next_to(ax.c2p(0.5 * (b["lo_au"] + b["hi_au"]), 16), UP, buff=0.1))
        neptune = dataio.body("Neptune")
        nep = widgets.marker_line(ax, neptune["a_au"], (0, 16), "Neptune", color=P.FG, side=UP)
        meas = Line(ax.c2p(80, pub["warp_deg"]), ax.c2p(200, pub["warp_deg"]), color=P.GREEN,
                    stroke_width=5)
        meas_lab = layout.label(f"measured ≈ {pub['warp_deg']:.0f}°", font_size=15, color=P.GREEN)
        meas_lab.next_to(ax.c2p(200, pub["warp_deg"]), RIGHT, buff=0.1)
        cols = [P.MUTED, P.TEAL, P.ORANGE]
        curves, tags = VGroup(), VGroup()
        for prof, col in zip(d["profiles"], cols):
            curves.add(widgets.curve(ax, ag, prof["tilt_deg"], color=col))
            tags.add(layout.label(f"{prof['mass_earth']:g} M⊕", font_size=14, color=col)
                     .next_to(ax.c2p(ag[-1], prof["tilt_deg"][-1]), RIGHT, buff=0.1))
        py = Dot(ax.c2p(earth["a_au"], 0), radius=0.09, color=P.TEAL).set_z_index(4)
        cap3 = layout.caption("The paper's measurement: the mean plane, bin by bin", font_size=22)
        self.play(FadeIn(ax), FadeIn(labs), FadeIn(bins), FadeIn(bin_tags), Create(nep),
                  FadeIn(cap3), run_time=0.9)
        self.play(Create(meas), FadeIn(meas_lab), run_time=0.7)
        timing.hold_to_read(self, cap3, bin_tags, settle=0.6)
        cap4 = layout.caption(
            f"Computed: the forced plane for a planet at {earth['a_au']:.0f} AU (dot), "
            f"tilted {d['i_planet_deg']:.0f}°", font_size=22)
        self.play(FadeIn(py), FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        self.play(*[Create(cv) for cv in curves], FadeIn(tags), run_time=1.8)
        timing.hold_to_read(self, cap4, settle=0.6)
        avg = VGroup()
        for b, v in zip(d["bins"], earth["bin_tilt_deg"]):
            avg.add(Line(ax.c2p(b["lo_au"], v), ax.c2p(b["hi_au"], v), color=P.ORANGE,
                         stroke_width=2).set_stroke(opacity=0.9))
        cap4b = layout.caption(
            f"Averaged over each bin, 1 M⊕ gives {earth['bin_tilt_deg'][1]:.0f}° at 80–200 AU "
            f"and {earth['bin_tilt_deg'][0]:.1f}° at 50–80 AU", font_size=22)
        self.play(FadeIn(avg), FadeOut(cap4), FadeIn(cap4b), run_time=0.6)
        timing.hold_to_read(self, cap4b, settle=1.2)
        self.play(FadeOut(VGroup(ax, labs, bins, bin_tags, nep, meas, meas_lab, curves, tags, py,
                                 avg, cap4b)), run_time=0.7)

        # 3. computed: which planets warp 80-200 AU yet leave 50-80 AU flat?
        mp = d["map"]
        vmax = pub["warp_deg"]
        left = TiltMap(mp["mass_earth"], mp["a_au"], mp["warp_bin_tilt_deg"], vmax,
                       title="tilt at 80–200 AU  (measured ≈ 15°)")
        right = TiltMap(mp["mass_earth"], mp["a_au"], mp["inner_bin_tilt_deg"], vmax,
                        title="tilt at 50–80 AU  (measured ≈ 0°)")
        left.move_to([-3.05, 0.35, 0])
        right.move_to([3.35, 0.35, 0])

        def furniture(m):
            g = VGroup()
            for mass, name in zip(pub["mass_range_earth"], ("Mercury", "Earth")):
                at = m.p(mass, m.a0)
                g.add(Line(at, at + DOWN * 0.08, color=P.MUTED, stroke_width=1.2))
                g.add(layout.label(name, font_size=13, color=P.MUTED).next_to(at, DOWN, buff=0.1))
            for a in (100, 200, 300):
                at = m.p(10 ** m.l0, a)
                g.add(Line(at, at + LEFT * 0.08, color=P.MUTED, stroke_width=1.2))
                g.add(layout.label(f"{a}", font_size=13, color=P.MUTED).next_to(at, LEFT, buff=0.1))
            lo, hi = pub["mass_range_earth"]
            alo, ahi = pub["a_range_au"]
            box = Polygon(m.p(lo, alo), m.p(hi, alo), m.p(hi, ahi), m.p(lo, ahi),
                          color=P.TEAL, stroke_width=2.4)
            g.add(box)
            return g, box

        fl, box_l = furniture(left)
        fr, box_r = furniture(right)
        xlab = layout.label("planet mass", font_size=15).move_to([0.15, -1.95, 0])
        ylab = layout.label("planet distance (AU)", font_size=15).rotate(np.pi / 2)
        ylab.next_to(left, LEFT, buff=0.45)
        swatches = VGroup()
        for v in (0, 5, 10, 15):
            sq = rect([0, 0, 0], [0.3, 0.22, 0], P.ORANGE, 0.05 + 0.85 * min(v / vmax, 1.0))
            swatches.add(VGroup(sq, layout.label(f"{v}°", font_size=13, color=P.MUTED)
                                .next_to(sq, RIGHT, buff=0.06)))
        swatches.arrange(RIGHT, buff=0.2)
        key = VGroup(layout.label("tilt:", font_size=15, color=P.MUTED), swatches,
                     layout.label("   box: the paper's preferred planets", font_size=15,
                                  color=P.TEAL)).arrange(RIGHT, buff=0.15)
        key.move_to([0.15, -2.4, 0])
        cap5 = layout.caption("Computed for every planet: heavier or closer ones warp the inner belt too",
                              font_size=22)
        self.play(FadeIn(left[:2]), FadeIn(right[:2]), FadeIn(left[2]), FadeIn(right[2]),
                  FadeIn(fl[:-1]), FadeIn(fr[:-1]), FadeIn(xlab), FadeIn(ylab), FadeIn(cap5),
                  run_time=1.2)
        timing.hold_to_read(self, cap5, settle=0.6)
        self.play(Create(box_l), Create(box_r), FadeIn(key), run_time=0.8)
        timing.hold_to_read(self, key, settle=1.2)
        self.play(FadeOut(cap5), run_time=0.4)

        layout.show_takeaway(
            self, "Flat inside 80 AU, warped beyond: a small planet at 100–200 AU, not Planet Nine")
