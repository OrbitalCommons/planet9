"""Siraj, Chyba & Tremaine (2026) -- measuring apsidal clustering.

A likelihood-ratio estimator: is the distribution of perihelion longitudes
better described by a bump (von Mises) than by a flat line? Each object casts a
vote, ln f_bump(ϖ) / f_flat, and the votes sum to ln Λ. Reproduced in
p9-2026-apsidal-clustering. The paper's 21- and 25-object tables are not in the
crate (its seeded stand-ins are tuned to the published answers, so they are not
drawn); the estimator is run instead on the workspace's vetted sample, before
and after the two 2025 discoveries. Samples, fits, votes and significances are
the crate's own (anim.json -> papers -> p9-2026-apsidal-clustering).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Create,
    DashedLine,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    ReplacementTransform,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2026-apsidal-clustering"

TOP = 0.5


def density_at(sample, deg):
    return float(np.interp(deg, sample["lon_deg"], sample["density"]))


def fit_curve(ax, sample, color):
    return widgets.curve(ax, sample["lon_deg"], sample["density"], color=color,
                         stroke_width=3.5)


def votes(ax, sample, flat, colors):
    """One stroke per object from the flat line to the bump at its ϖ: green
    where the bump is higher (a vote for clustering), red where it is lower."""
    g = VGroup()
    for deg, vote, col in zip(sample["varpi_deg"], sample["votes"], colors):
        y = density_at(sample, deg)
        c = P.GREEN if vote > 0 else P.RED
        g.add(VGroup(
            Line(ax.c2p(deg, flat), ax.c2p(deg, y), color=c, stroke_width=5),
            Line(ax.c2p(deg, 0.0), ax.c2p(deg, 0.035), color=col, stroke_width=4),
        ))
    return g


def verdict(sample, color, note):
    rows = VGroup(
        layout.label(f"n = {sample['n']}", font_size=20, color=color, weight="BOLD"),
        layout.label(f"sum of votes  ln Λ = {sample['lambda']:.2f}", font_size=17,
                     color=P.FG),
        layout.label(f"{sample['sigma_two_sided']:.1f}σ", font_size=34, color=color,
                     weight="BOLD"),
        layout.label(note, font_size=15, color=P.FG, line_spacing=0.9),
    ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
    return rows.move_to([3.95, 0.7, 0], aligned_edge=LEFT)


class ApsidalClustering2026(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        real, plus = d["real"], d["real_plus"]
        flat = d["uniform_density"]
        n0 = real["n"]

        self.add(paper.scene_header(CRATE))

        ax, labels = widgets.labeled_axes(
            [0, 360, 60], [0, TOP, 0.1], x_label="longitude of perihelion ϖ (degrees)",
            y_label="probability density (per radian)", y_rotate=True, numbers=True,
            x_length=8.0, y_length=4.1, shift_down=-0.35)
        VGroup(ax, labels).shift(LEFT * 2.1)
        flat_line = DashedLine(ax.c2p(0, flat), ax.c2p(360, flat), color=P.MUTED,
                               stroke_width=2.5)
        flat_lab = layout.label("flat: no clustering", font_size=15, color=P.FG)
        flat_lab.next_to(ax.c2p(200, flat), UP, buff=0.1)

        # 1. flat line or bump?
        ticks = VGroup(*[Line(ax.c2p(v, 0.0), ax.c2p(v, 0.035), color=P.GREEN, stroke_width=4)
                         for v in real["varpi_deg"]])
        bump = fit_curve(ax, real, P.GREEN)
        cap = layout.caption(
            f"The perihelia of {n0} distant objects: a flat line, or a bump?", font_size=22)
        self.play(Create(ax), FadeIn(labels), FadeIn(cap))
        self.play(FadeIn(ticks, lag_ratio=0.1), Create(flat_line), FadeIn(flat_lab))
        self.play(Create(bump), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.6)

        # 2. every object votes
        v_real = votes(ax, real, flat, [P.GREEN] * n0)
        verdict_real = verdict(real, P.GREEN,
                               f"paper, {d['paper_n_before']} stable objects:\n"
                               f"{d['paper_sigma21']:.1f}σ")
        cap2 = layout.caption(
            "Each object votes: bump above the flat line counts for clustering", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), FadeOut(ticks))
        self.play(LaggedStart(*[Create(x) for x in v_real], lag_ratio=0.15), run_time=1.8)
        self.play(FadeIn(verdict_real, lag_ratio=0.1))
        timing.hold_to_read(self, cap2, verdict_real, settle=1.0)

        # 3. two newcomers that sit where the bump is low
        new_cols = [P.GREEN] * n0 + [P.RED] * (plus["n"] - n0)
        newcomers = VGroup(*[
            Line(ax.c2p(v, 0.0), ax.c2p(v, 0.06), color=P.RED, stroke_width=6)
            for v in plus["varpi_deg"][n0:]])
        new_lab = layout.label("2017 OF201 and Ammonite", font_size=15, color=P.RED)
        new_lab.move_to(ax.c2p(288, 0.30))
        cap3 = layout.caption(
            "Add the two 2025 discoveries: they vote against, and the bump sags",
            font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), FadeIn(newcomers, lag_ratio=0.3),
                  FadeIn(new_lab))
        timing.hold_to_read(self, cap3, settle=0.4)
        v_plus = votes(ax, plus, flat, new_cols)
        verdict_plus = verdict(plus, P.RED,
                               f"paper, {d['paper_n_sample']} stable objects:\n"
                               f"{d['paper_sigma25']:.1f}σ")
        self.play(ReplacementTransform(bump, fit_curve(ax, plus, P.RED)),
                  ReplacementTransform(v_real, v_plus), FadeOut(newcomers), FadeOut(new_lab),
                  ReplacementTransform(verdict_real, verdict_plus), run_time=2.0)
        timing.hold_to_read(self, cap3, verdict_plus, settle=1.4)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, f"Off-cluster newcomers dilute it: {real['sigma_two_sided']:.1f}σ to "
                  f"{plus['sigma_two_sided']:.1f}σ here, {d['paper_sigma21']:.1f}σ to "
                  f"{d['paper_sigma25']:.1f}σ in the paper.")
