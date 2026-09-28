"""Clement & Kaib (2020) -- orbital precession in the distant solar system.

The giant planets make every distant orbit precess, each at its own pace, so a
cluster of apsides left to itself scrambles within a few hundred Myr; something
must keep re-aligning it. The paper's N-body runs find that Planet Nine does,
but also that 17 detections are too few: random orbits drawn 17 at a time still
look partly clustered. Reproduced in p9-2020-clement-precession; the precession
periods, the drift of the clustering and the random-sample distribution are the
crate's own (anim.json -> papers -> p9-2020-clement-precession).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Arrow,
    Create,
    FadeIn,
    FadeOut,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2020-clement-precession"
CENTRE = np.array([-3.7, 0.0, 0.0])
ARROW_LEN = 2.3  # scene units for R-bar = 1


def display_orbits(objs, reach=2.75):
    """(a_display, e) per object: a range-compressed so every orbit fits within
    ``reach`` of the Sun at aphelion, e exact."""
    amax = max(o["a_au"] for o in objs)
    raw = [(np.sqrt(o["a_au"] / amax), o["e"]) for o in objs]
    k = reach / max(a * (1 + e) for a, e in raw)
    return [(k * a, e) for a, e in raw]


def mean_arrow(varpis_deg, color=P.FG):
    v = np.deg2rad(np.asarray(varpis_deg))
    mx, my = np.mean(np.cos(v)), np.mean(np.sin(v))
    tip = CENTRE + ARROW_LEN * np.array([mx, my, 0.0])
    return Arrow(CENTRE, tip, buff=0, color=color, stroke_width=6,
                 max_tip_length_to_length_ratio=0.3).set_z_index(6)


def unit_arrows(varpis_deg, color, length=1.5, width=2.5):
    g = VGroup()
    for w in np.deg2rad(varpis_deg):
        tip = CENTRE + length * np.array([np.cos(w), np.sin(w), 0.0])
        g.add(Arrow(CENTRE, tip, buff=0, color=color, stroke_width=width,
                    max_tip_length_to_length_ratio=0.12).set_opacity(0.8))
    return g


class ClementPrecession2020(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        varpi0 = np.array([o["varpi_deg"] for o in objs])
        periods = np.array([o["giant_period_gyr"] for o in objs])
        shapes = display_orbits(objs)
        drift_t = np.array(d["drift"]["t_gyr"])
        drift_r = np.array(d["drift"]["r_bar"])
        r_obs = d["r_bar_observed"]

        self.add(paper.scene_header(CRATE))

        # 1. the real sample, and what "clustered" means
        sun = orbits.sun().move_to(CENTRE)
        t = ValueTracker(0.0)

        def varpis_now():
            return varpi0 + 360.0 * t.get_value() / periods

        def swarm():
            g = VGroup()
            for (a, e), w in zip(shapes, varpis_now()):
                g.add(orbits.ellipse_orbit(a, e, color=P.GREEN, varpi=np.deg2rad(w),
                                           stroke_width=1.8, opacity=0.8).shift(CENTRE))
            return g

        live = always_redraw(swarm)
        self.play(FadeIn(sun), Create(live), run_time=1.6)
        cap = layout.caption(f"The {len(objs)} distant objects of the 2017 sample: "
                             "their perihelia bunch to one side", font_size=22)
        self.play(FadeIn(cap))
        timing.hold_to_read(self, cap, settle=0.6)

        arrow = always_redraw(lambda: mean_arrow(varpis_now()))
        rlab = always_redraw(lambda: layout.label(
            f"R̄ = {np.interp(t.get_value(), drift_t, drift_r):.2f}", font_size=22,
            color=P.FG, weight="BOLD").move_to([-6.2, 2.75, 0], aligned_edge=LEFT))
        cap2 = layout.caption("Average the perihelion directions: arrow length R̄ "
                              "(0 = random, 1 = all aligned)", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), FadeIn(arrow), FadeIn(rlab))
        timing.hold_to_read(self, cap2, settle=0.8)

        # 2. the giant planets turn every orbit at its own pace
        ax, labels = widgets.labeled_axes(
            [0, 2, 0.5], [0, 1, 0.25], x_label="time (Gyr)", y_label="clustering R̄",
            y_rotate=True, numbers=True, x_length=5.6, y_length=3.6, shift_down=0)
        VGroup(ax, labels).move_to([3.4, 0.25, 0])
        t_max = 2.0
        trace = always_redraw(lambda: widgets.curve(
            ax, drift_t[drift_t <= max(t.get_value(), 0.02) + 1e-9],
            drift_r[drift_t <= max(t.get_value(), 0.02) + 1e-9], color=P.ORANGE))
        clock = always_redraw(lambda: layout.label(
            f"t = {1000 * t.get_value():.0f} Myr", font_size=20, color=P.MUTED)
            .next_to(ax.c2p(2.0, 1.0), DOWN + LEFT, buff=0.1))
        cap3 = layout.caption(
            f"Neptune and the giants turn each orbit once every "
            f"{periods.min():.2f} to {periods.max():.1f} Gyr", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), Create(ax), FadeIn(labels))
        self.add(trace, clock)
        timing.hold_to_read(self, cap3, settle=0.3)
        self.play(t.animate.set_value(t_max), run_time=9.0, rate_func=linear)

        null_mean = d["r_bar_null_mean_sample"]
        scrambled = drift_t[np.argmax(drift_r < null_mean)]
        level = widgets.curve(ax, [0, t_max], [null_mean, null_mean], color=P.RED,
                              stroke_width=2).set_stroke(opacity=0.7)
        level_lab = layout.label(f"typical of {len(objs)} random orbits", font_size=15,
                                 color=P.RED).next_to(ax.c2p(t_max, null_mean), UP + LEFT,
                                                      buff=0.08)
        cap4 = layout.caption(
            f"Unheld, it sinks to the random level within ~{1000 * scrambled:.0f} Myr, "
            "a blink beside 4.5 Gyr", font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), Create(level), FadeIn(level_lab))
        timing.hold_to_read(self, cap4, settle=0.6)
        cap5 = layout.caption("Paper's N-body runs: with Planet Nine, >60% of "
                              "detectable survivors stay anti-aligned", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=0.8)
        for m in (trace, clock, live, arrow, rlab):
            m.clear_updaters()
        self.play(FadeOut(VGroup(ax, labels, sun, cap5, trace, clock, live, arrow, rlab,
                                 level, level_lab)),
                  run_time=0.8)

        # 3. but seventeen orbits are a small sample
        ex = d["examples_17"]
        n17 = d["null"]["n_paper"]
        sun2 = orbits.sun().move_to(CENTRE)
        spokes = unit_arrows(ex[0]["varpi_deg"], P.MUTED)
        marrow = mean_arrow(ex[0]["varpi_deg"])
        exlab = layout.label(f"R̄ = {ex[0]['r_bar']:.2f}", font_size=22, color=P.FG,
                             weight="BOLD").move_to(CENTRE + np.array([0.0, -2.3, 0.0]))
        cap6 = layout.caption(f"Now {n17} orbits pointing at random: is R̄ zero?",
                              font_size=22)
        self.play(FadeIn(sun2), Create(spokes), FadeIn(cap6), run_time=1.2)
        self.play(FadeIn(marrow), FadeIn(exlab))
        timing.hold_to_read(self, cap6, settle=0.3)
        for draw in ex[1:]:
            new_sp = unit_arrows(draw["varpi_deg"], P.MUTED)
            new_ar = mean_arrow(draw["varpi_deg"])
            new_lab = layout.label(f"R̄ = {draw['r_bar']:.2f}", font_size=22, color=P.FG,
                                   weight="BOLD").move_to(exlab)
            self.play(FadeOut(spokes), FadeOut(marrow), FadeOut(exlab),
                      FadeIn(new_sp), FadeIn(new_ar), FadeIn(new_lab), run_time=0.8)
            self.wait(0.6)
            spokes, marrow, exlab = new_sp, new_ar, new_lab

        edges = np.array(d["null"]["edges"])
        frac = np.array(d["null"]["fraction_paper"])
        top = 0.16
        ax2, lab2 = widgets.labeled_axes(
            [0, 0.8, 0.2], [0, top, 0.04], x_label="clustering R̄ of 17 random orbits",
            y_label="fraction of draws", y_rotate=True, numbers=True,
            x_length=5.6, y_length=3.5, shift_down=0)
        VGroup(ax2, lab2).move_to([3.4, 0.35, 0])
        bars = widgets.histogram(ax2, edges, frac, color=P.MUTED, opacity=0.6)
        mean = d["r_bar_null_mean_17"]
        p95 = d["r_bar_null_p95_17"]
        m1 = widgets.marker_line(ax2, mean, (0, top), f"mean {mean:.2f}", color=P.RED,
                                 side=UP)
        m2 = widgets.marker_line(ax2, p95, (0, top * 0.8), f"1 in 20 exceed {p95:.2f}",
                                 color=P.RED, side=UP)
        cap7 = layout.caption(f"{d['null']['n_draws']:,} random draws: "
                              "a small sample is never unclustered", font_size=22)
        self.play(FadeOut(cap6), FadeIn(cap7), Create(ax2), FadeIn(lab2))
        self.play(FadeIn(bars, lag_ratio=0.05), run_time=1.2)
        self.play(Create(m1), Create(m2))
        timing.hold_to_read(self, cap7, settle=1.0)
        self.play(FadeOut(cap7))

        layout.show_takeaway(
            self, f"Random 17-orbit samples still show R̄ ≈ {mean:.2f}: "
                  "more detections are needed.")

