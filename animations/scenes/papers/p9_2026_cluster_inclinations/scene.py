"""Bansal et al. (2026) -- distant TNO inclinations constrain the birth cluster.

A violent stellar encounter in the Sun's birth cluster would have left the
distant orbits steeply tilted; Planet Nine might have cooled them since. The
paper measures how cold the observed high-perihelion objects really are,
correcting for where each was found, and integrates both kinds of primordial
population beside Planet Nine. Reproduced in p9-2026-cluster-inclinations:
the sample in its mean plane, the debiasing (Brown 2001) probabilities, the
Kuiper-statistic scan and the width histories of a short, small integration
of each population are the crate's own (anim.json -> papers ->
p9-2026-cluster-inclinations); the paper's 4 Gyr widths are labelled bands.
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
    Polygon,
    Rectangle,
    Scene,
    Transform,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2026-cluster-inclinations"


def band(ax, x0, x1, y0, y1, color, opacity):
    p0, p1 = np.array(ax.c2p(x0, y0)), np.array(ax.c2p(x1, y1))
    r = Rectangle(width=abs(p1[0] - p0[0]), height=abs(p1[1] - p0[1]), stroke_width=0)
    return r.set_fill(color, opacity=opacity).move_to((p0 + p1) / 2)


def pp_dots(ax, probs, color):
    n = len(probs)
    return VGroup(*[Dot(ax.c2p((k + 1) / (n + 1), p), radius=0.06, color=color)
                    for k, p in enumerate(probs)])


class ClusterInclinations2026(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        pub = d["published"]
        w_best, w_lo, w_hi = d["w_best_deg"], d["w_lo_deg"], d["w_hi_deg"]
        self.add(paper.scene_header(CRATE))

        # 1. the sample, and the catch: where each was found
        ax, labels = widgets.labeled_axes(
            [0, 25, 5], [0, 25, 5], x_label="latitude where it was found, |β| (deg)",
            y_label="inclination to the mean plane, i (deg)", y_rotate=True, numbers=True,
            x_length=4.6, y_length=4.6, shift_down=0)
        VGroup(ax, labels).move_to([-3.2, 0.1, 0])
        forbid = Polygon(ax.c2p(0, 0), ax.c2p(25, 25), ax.c2p(25, 0), stroke_width=0)
        forbid.set_fill(P.RED, opacity=0.10)
        forbid_l = layout.label("impossible:\ni < |β|", font_size=15, color=P.RED)
        forbid_l.move_to(ax.c2p(18, 6))
        dots = VGroup(*[Dot(ax.c2p(o["beta_deg"], o["i_deg"]), radius=0.07, color=P.GREEN)
                        .set_z_index(4) for o in objs])
        text = VGroup(
            layout.label(f"{d['n_observed']} distant objects that never", font_size=19),
            layout.label("come near Neptune (q = 40–80 AU)", font_size=19),
            layout.label("An orbit tilted by i spends most", font_size=19, color=P.MUTED),
            layout.label("of its time near latitude ±i, but", font_size=19, color=P.MUTED),
            layout.label("surveys look near the ecliptic:", font_size=19, color=P.MUTED),
            layout.label("low tilts are over-counted", font_size=19, color=P.MUTED),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT).move_to([3.4, 0.4, 0])
        text[2].shift(DOWN * 0.25)
        text[3:].shift(DOWN * 0.25)
        cap = layout.caption("How tilted are the distant orbits, really?", font_size=22)
        self.play(Create(ax), FadeIn(labels), FadeIn(cap))
        self.play(FadeIn(dots, lag_ratio=0.08), FadeIn(text[:2]), run_time=1.4)
        self.play(FadeIn(forbid), FadeIn(forbid_l), FadeIn(text[2:]))
        cap2 = layout.caption("So judge each tilt against the latitude where that object "
                              "was found (Brown 2001)", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        timing.hold_to_read(self, cap2, text, settle=0.6)
        self.play(FadeOut(VGroup(ax, labels, forbid, forbid_l, dots, text, cap2)))

        # 2. testing one width: are the probabilities uniform?
        ax2, lab2 = widgets.labeled_axes(
            [0, 1, 0.25], [0, 1, 0.25], x_label="expected if the width is right",
            y_label="each object's probability", y_rotate=True, numbers=True,
            x_length=4.4, y_length=4.4, shift_down=0)
        VGroup(ax2, lab2).move_to([-3.3, 0.15, 0])
        diag = DashedLine(ax2.c2p(0, 0), ax2.c2p(1, 1), color=P.MUTED, stroke_width=2)
        pts26 = pp_dots(ax2, d["probabilities_26"], P.RED)
        pts12 = pp_dots(ax2, d["probabilities_best"], P.GREEN)
        w26 = layout.label("try a width of 26°: points sag below the line", font_size=18,
                           color=P.RED).move_to([2.9, 1.6, 0])
        w12 = layout.label(f"try {w_best:.0f}°: they follow it", font_size=18,
                           color=P.GREEN).move_to([2.9, 0.9, 0])
        cap3 = layout.caption("Guess an intrinsic width w; if it is right, each object's "
                              "probability is uniform", font_size=22)
        self.play(Create(ax2), FadeIn(lab2), Create(diag), FadeIn(cap3))
        self.play(FadeIn(pts26, lag_ratio=0.05), FadeIn(w26), run_time=1.2)
        timing.hold_to_read(self, cap3, w26, settle=0.3)
        self.play(Transform(pts26, pts12), FadeIn(w12), run_time=1.4)
        timing.hold_to_read(self, w12, settle=0.8)
        self.play(FadeOut(VGroup(ax2, lab2, diag, pts26, w26, w12, cap3)))

        # 3. the scan over widths
        sc = d["scan"]
        ws = np.array(sc["w_deg"])
        pv = np.array(sc["p_value"])
        ax3, lab3 = widgets.labeled_axes(
            [0, 45, 5], [0, 0.3, 0.1], x_label="intrinsic inclination width w (deg)",
            y_label="how well it fits (Kuiper p)", y_rotate=True, numbers=True,
            x_length=8.6, y_length=4.0, shift_down=0)
        VGroup(ax3, lab3).move_to([-1.0, 0.3, 0])
        one_sig = band(ax3, w_lo, w_hi, 0, 0.3, P.GREEN, 0.14)
        curve = widgets.curve(ax3, ws, pv, color=P.GREEN, stroke_width=3.5)
        thr = DashedLine(ax3.c2p(0, d["p_one_sigma"]), ax3.c2p(45, d["p_one_sigma"]),
                         color=P.MUTED, stroke_width=1.5)
        thr_l = layout.label("1σ", font_size=15, color=P.MUTED).next_to(
            ax3.c2p(45, d["p_one_sigma"]), RIGHT, buff=0.1)
        best = widgets.marker_line(ax3, w_best, (0, 0.3),
                                   f"w = {w_best:.0f}° ({w_lo:.0f}–{w_hi:.0f}°)",
                                   color=P.GREEN, side=UP)
        hot = widgets.marker_line(ax3, 26.0, (0, 0.2), "26°: cluster-stirred", color=P.RED,
                                  side=UP)
        verdict = VGroup(
            layout.label(f"p(26°) = {d['p_at_26']:.4f}", font_size=19, color=P.RED),
            layout.label(f"paper {pub['p_reject_w26']:.4f}: ~3σ", font_size=16, color=P.FG),
        ).arrange(DOWN, buff=0.08, aligned_edge=LEFT).move_to([4.9, 1.2, 0])
        cap4 = layout.caption("Scan every width from 3° to 45°", font_size=22)
        self.play(Create(ax3), FadeIn(lab3), FadeIn(cap4))
        self.play(Create(curve), Create(thr), FadeIn(thr_l), run_time=1.6)
        cap5 = layout.caption(f"Best width {w_best:.0f}° (paper {pub['w_obs_deg']:.0f}°); "
                              "a cluster-stirred 26° is ruled out", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), run_time=0.6)
        self.play(FadeIn(one_sig), Create(best))
        self.play(Create(hot), FadeIn(verdict))
        timing.hold_to_read(self, cap5, verdict, settle=0.8)
        self.play(FadeOut(VGroup(ax3, lab3, one_sig, curve, thr, thr_l, best, hot, verdict,
                                 cap5)))

        # 4. can Planet Nine cool a stirred population?
        hot_pop, cold_pop = d["cluster_influenced"], d["cluster_free"]
        t_end = d["t_myr"]
        ax4, lab4 = widgets.labeled_axes(
            [0, t_end, 20], [0, 40, 10], x_label="time beside Planet Nine (Myr)",
            y_label="inclination width w (deg)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=0)
        VGroup(ax4, lab4).move_to([-2.0, 0.25, 0])
        obs = band(ax4, 0, t_end, w_lo, w_hi, P.GREEN, 0.16)
        obs_l = layout.label(f"observed {w_best:.0f}°", font_size=16, color=P.GREEN).next_to(
            ax4.c2p(t_end, w_best), RIGHT, buff=0.1)
        c_hot = widgets.curve(ax4, hot_pop["t_myr"], hot_pop["w_deg"], color=P.RED)
        c_cold = widgets.curve(ax4, cold_pop["t_myr"], cold_pop["w_deg"], color=P.TEAL)
        l_hot = layout.label("cluster-stirred start", font_size=16, color=P.RED).next_to(
            ax4.c2p(2, hot_pop["w_deg"][0] + 6), RIGHT, buff=0.0)
        l_cold = layout.label("quiet start", font_size=16, color=P.TEAL).next_to(
            ax4.c2p(0, cold_pop["w_deg"][0]), DOWN + RIGHT, buff=0.08)
        paper_key = VGroup(
            layout.label("paper, after 4 Gyr:", font_size=16, color=P.MUTED),
            layout.label(f"stirred {pub['w_cluster_influenced_deg'][0]:.0f}–"
                         f"{pub['w_cluster_influenced_deg'][1]:.1f}°", font_size=17,
                         color=P.RED),
            layout.label(f"quiet {pub['w_cluster_free_deg'][0]:.0f}–"
                         f"{pub['w_cluster_free_deg'][1]:.1f}°", font_size=17, color=P.TEAL),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT).move_to([4.6, 1.5, 0])
        run_note = layout.label(f"here: {hot_pop['n_particles']} + {cold_pop['n_particles']} "
                                f"particles, {t_end:.0f} Myr", font_size=14, color=P.MUTED)
        run_note.next_to(paper_key, DOWN, buff=0.3, aligned_edge=LEFT)
        cap6 = layout.caption("Evolve both populations with the giant planets, Planet Nine, "
                              "passing stars and the tide", font_size=22)
        self.play(Create(ax4), FadeIn(lab4), FadeIn(obs), FadeIn(obs_l), FadeIn(cap6))
        self.play(Create(c_hot), Create(c_cold), FadeIn(l_hot), FadeIn(l_cold), run_time=2.4)
        self.play(FadeIn(paper_key), FadeIn(run_note))
        timing.hold_to_read(self, cap6, settle=0.3)
        cap7 = layout.caption("Planet Nine never cools the stirred one; the quiet one warms, "
                              "but stays far below 26°", font_size=22)
        self.play(FadeOut(cap6), FadeIn(cap7), run_time=0.6)
        timing.hold_to_read(self, cap7, paper_key, settle=0.8)
        self.play(FadeOut(cap7))

        layout.show_takeaway(
            self, f"Cold {w_best:.0f}° orbits, even with Planet Nine: the birth cluster "
                  "was gentle.")
