"""Brown & Batygin (2021) -- the orbit of Planet Nine.

Two results: the clustering of the distant orbits survives an updated survey-bias
null, and inverting it through a grid of simulations gives a posterior for the
planet's mass and orbit. Reproduced in p9-2021-orbit: the clustering confidence
is recomputed from the vetted 10-object table, and the published posterior is
resampled by the workspace's emulator. Everything drawn is the crate's own
(anim.json -> papers -> p9-2021-orbit).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Arrow,
    Create,
    GrowArrow,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2021-orbit"


def fit_paths(paths, box):
    """Scale and origin that fit every (x, y) path in AU into `box`
    (x0, x1, y0, y1) in scene units."""
    pts = np.array([p for path in paths for p in path])
    lo, hi = pts.min(axis=0), pts.max(axis=0)
    s = min((box[1] - box[0]) / (hi[0] - lo[0]), (box[3] - box[2]) / (hi[1] - lo[1]))
    mid = 0.5 * (lo + hi)
    origin = np.array([0.5 * (box[0] + box[1]) - s * mid[0],
                       0.5 * (box[2] + box[3]) - s * mid[1], 0.0])
    return s, origin


def path_curve(path, s, origin, color, stroke_width=2.0, opacity=1.0):
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners([origin + s * np.array([x, y, 0.0]) for x, y in path])
    return m.set_stroke(opacity=opacity)


def unit(deg):
    t = np.deg2rad(deg)
    return np.array([np.cos(t), np.sin(t), 0.0])


def interval_text(name, v, unit, digits=0):
    f = f"{{:.{digits}f}}"
    if abs(v["plus"] - v["minus"]) < 1e-9:
        return f"{name} = {f.format(v['median'])} ± {f.format(v['plus'])} {unit}"
    return (f"{name} = {f.format(v['median'])}  (+{f.format(v['plus'])} "
            f"−{f.format(v['minus'])}) {unit}")


class Orbit2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        nul = d["null_r_bar"]
        post = d["posterior"]
        n = d["n_sample"]

        self.add(paper.scene_header(CRATE))

        # 1. the clustering, against a null that knows where surveys look
        top = 3.5
        ax, labels = widgets.labeled_axes(
            [0, 1, 0.2], [0, top, 1], x_label=f"alignment strength R of {n} longitudes",
            y_label="probability density", y_rotate=True, numbers=True,
            x_length=8.2, y_length=4.1, shift_down=-0.35)
        VGroup(ax, labels).shift(LEFT * 2.0)
        hist = widgets.histogram(ax, nul["edges"], nul["biased"], color=P.MUTED, opacity=0.6)
        obs = widgets.marker_line(ax, d["r_bar"], (0, top), f"observed  R = {d['r_bar']:.2f}",
                                  color=P.GREEN, font_size=16, side=UP)
        cap = layout.caption(
            "Random orbits, found by surveys that avoid the galactic plane, rarely align",
            font_size=22)
        self.play(Create(ax), FadeIn(labels), FadeIn(cap))
        self.play(FadeIn(hist, lag_ratio=0.05), run_time=1.2)
        self.play(Create(obs))
        conf = paper.result_readout("confidence the clustering is real",
                                    f"{100 * d['confidence_mc_bias']:.1f}%", color=P.GREEN)
        conf.scale(0.85).move_to([4.6, 1.0, 0])
        ref = layout.label(
            f"paper: {100 * d['paper_confidence']:.1f}% from {d['paper_n_sample']} objects\n"
            f"here: the {n} in the vetted table", font_size=17, color=P.FG, line_spacing=0.9)
        ref.next_to(conf, DOWN, buff=0.3)
        self.play(FadeIn(conf), FadeIn(ref))
        timing.hold_to_read(self, cap, conf, ref, settle=1.0)
        self.play(FadeOut(VGroup(ax, labels, hist, obs, conf, ref, cap)))

        # 2. the posterior
        samples = post["samples"]
        ax2, labels2 = widgets.labeled_axes(
            [200, 900, 100], [2, 14, 2], x_label="semi-major axis a (AU)",
            y_label="mass (Earth masses)", y_rotate=True, numbers=True,
            x_length=7.6, y_length=4.2, shift_down=-0.35)
        VGroup(ax2, labels2).shift(LEFT * 2.4)
        cloud = VGroup(*[
            Dot(ax2.c2p(s["a"], s["mass"]), radius=0.022, color=P.BLUE).set_opacity(0.45)
            for s in samples if 200 <= s["a"] <= 900 and 2 <= s["mass"] <= 14])
        a, m = post["a"], post["mass"]
        cross = VGroup(
            Line(ax2.c2p(a["median"] - a["minus"], m["median"]),
                 ax2.c2p(a["median"] + a["plus"], m["median"]), color=P.FG, stroke_width=3),
            Line(ax2.c2p(a["median"], m["median"] - m["minus"]),
                 ax2.c2p(a["median"], m["median"] + m["plus"]), color=P.FG, stroke_width=3),
            Dot(ax2.c2p(a["median"], m["median"]), radius=0.08, color=P.FG),
        ).set_z_index(4)
        rows = VGroup(
            layout.label(interval_text("mass", m, "M⊕", 1), font_size=20, color=P.BLUE,
                         weight="BOLD"),
            layout.label(interval_text("a", a, "AU"), font_size=20, color=P.BLUE, weight="BOLD"),
            layout.label(interval_text("perihelion", post["q"], "AU"), font_size=20,
                         color=P.BLUE, weight="BOLD"),
            layout.label(interval_text("inclination", post["i"], "°"), font_size=20,
                         color=P.BLUE, weight="BOLD"),
            layout.label(f"{len(samples)} draws from the published\nposterior (each dot one planet)",
                         font_size=15, color=P.FG, line_spacing=0.9),
        ).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
        rows.move_to([4.5, 0.4, 0])
        cap2 = layout.caption(
            "Simulations that best reproduce the real orbits pick out the planet", font_size=22)
        self.play(Create(ax2), FadeIn(labels2), FadeIn(cap2))
        self.play(FadeIn(cloud, lag_ratio=0.002), run_time=1.8)
        self.play(FadeIn(cross), FadeIn(rows, lag_ratio=0.15))
        timing.hold_to_read(self, cap2, rows, settle=1.2)
        self.play(FadeOut(VGroup(ax2, labels2, cloud, cross, rows, cap2)))

        # 3. the orbit among the orbits it shepherds
        etno_paths = [e["path"] for e in d["etnos"]]
        s, origin = fit_paths(etno_paths + [d["p9_path"]], (-6.4, 1.2, -2.45, 3.05))
        sun = orbits.sun(radius=0.07).move_to(origin)
        swarm = VGroup(*[path_curve(p, s, origin, P.GREEN, 1.8, 0.8) for p in etno_paths])
        p9 = path_curve(d["p9_path"], s, origin, P.BLUE, 3.5)
        bar_au = 500.0
        bar = Line(origin + np.array([0.0, 0.0, 0.0]), origin + np.array([s * bar_au, 0.0, 0.0]),
                   color=P.MUTED, stroke_width=2).move_to([4.2, -2.3, 0])
        bar_lab = layout.label(f"{bar_au:.0f} AU", font_size=14, color=P.MUTED)
        bar_lab.next_to(bar, UP, buff=0.08)
        legend = VGroup(
            layout.label(f"the {n} distant objects", font_size=19, color=P.GREEN),
            layout.label(f"perihelia toward ϖ = {d['mean_varpi_deg']:.0f}°", font_size=16,
                         color=P.FG),
            layout.label("Planet Nine, median orbit", font_size=19, color=P.BLUE, weight="BOLD"),
            layout.label(f"a = {d['a']:.0f} AU,  e = {d['e']:.2f},  i = {d['i']:.0f}°\n"
                         f"perihelion toward ϖ = {d['varpi_deg']:.0f}°", font_size=16,
                         color=P.FG, line_spacing=0.9),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        legend[2].shift(DOWN * 0.2)
        legend[3].shift(DOWN * 0.2)
        legend.move_to([4.4, 1.0, 0])
        cap3 = layout.caption("Seen from above: the planet's orbit points the opposite way",
                              font_size=22)
        self.play(FadeIn(sun), FadeIn(cap3))
        self.play(LaggedStart(*[Create(c) for c in swarm], lag_ratio=0.12), run_time=2.2)
        self.play(FadeIn(legend[0]), FadeIn(legend[1]), FadeIn(bar), FadeIn(bar_lab))
        self.play(Create(p9), run_time=1.6)
        arrows = VGroup(
            Arrow(origin, origin + 1.6 * unit(d["mean_varpi_deg"]), buff=0, color=P.GREEN,
                  stroke_width=5, max_tip_length_to_length_ratio=0.2),
            Arrow(origin, origin + 1.6 * unit(d["varpi_deg"]), buff=0, color=P.BLUE,
                  stroke_width=5, max_tip_length_to_length_ratio=0.2),
        ).set_z_index(6)
        self.play(FadeIn(legend[2]), FadeIn(legend[3]), GrowArrow(arrows[0]),
                  GrowArrow(arrows[1]))
        timing.hold_to_read(self, cap3, legend, settle=1.2)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, f"{d['mass']:.1f} Earth masses, {d['a']:.0f} AU out, tilted {d['i']:.0f}°: "
                  "closer and brighter than thought.")
