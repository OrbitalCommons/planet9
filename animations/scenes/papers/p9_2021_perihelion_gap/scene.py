"""Oldroyd & Trujillo (2021) -- the perihelion gap at 50-65 AU.

Plot the distant, eccentric orbits by perihelion and a lane between 50 and
65 AU is nearly empty. A single smooth population rarely leaves such a lane;
two populations -- Neptune-coupled extreme TNOs below, detached inner Oort
cloud objects above -- do, and a distant planet that cycles perihelia through
the lane quickly keeps it thin. Reproduced in p9-2021-perihelion-gap: the
sample, the smooth-population expectation and the Monte Carlo odds are the
crate's own (anim.json -> papers -> p9-2021-perihelion-gap).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Rectangle,
    Scene,
    VGroup,
    Transform,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2021-perihelion-gap"


def lane(ax, lo, hi, x0, x1, color=P.RED, opacity=0.12):
    p0, p1 = np.array(ax.c2p(x0, lo)), np.array(ax.c2p(x1, hi))
    r = Rectangle(width=p1[0] - p0[0], height=p1[1] - p0[1], stroke_width=0)
    return r.set_fill(color, opacity=opacity).move_to((p0 + p1) / 2)


class PerihelionGap2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        g_lo, g_hi = d["gap_au"]
        objs = d["objects"]
        pe = d["paper_epoch"]
        self.add(paper.scene_header(CRATE))

        # 1. the sample, and the empty lane
        ax = widgets.axes([2.1, 3.5, 0.1], [30, 90, 10], x_length=8.4, y_length=4.6,
                          shift_down=0)
        ax.move_to([-1.3, 0.2, 0])
        ax.get_y_axis().add_numbers(font_size=16)
        xt = VGroup(*[layout.label(str(v), font_size=15, color=P.FG)
                      .next_to(ax.c2p(np.log10(v), 30), DOWN, buff=0.12)
                      for v in (150, 300, 1000, 3000)])
        xl = layout.label("semi-major axis a (AU, log scale)", font_size=18).next_to(
            xt, DOWN, buff=0.12).set_x(ax.get_center()[0])
        yl = layout.label("perihelion q (AU)", font_size=15).rotate(np.pi / 2).next_to(
            ax, LEFT, buff=0.4)
        band = lane(ax, g_lo, g_hi, 2.1, 3.5)
        band_l = layout.label(f"the gap: q = {g_lo:.0f}–{g_hi:.0f} AU", font_size=18,
                              color=P.RED).next_to(ax.c2p(3.5, (g_lo + g_hi) / 2), RIGHT,
                                                   buff=0.15)
        old = [o for o in objs if not o["post_paper"]]
        new = [o for o in objs if o["post_paper"]]
        dots = VGroup(*[Dot(ax.c2p(np.log10(o["a_au"]), o["q_au"]), radius=0.07,
                            color=P.GREEN) for o in old])
        self.play(Create(ax), FadeIn(xt), FadeIn(xl), FadeIn(yl))
        cap = layout.caption(f"{len(old)} known distant orbits (e > {d['e_floor']:.2f}) "
                             "at the time of the paper", font_size=22)
        self.play(FadeIn(dots, lag_ratio=0.08), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.3)
        low_l = layout.label("Neptune-coupled\nextreme TNOs", font_size=16, color=P.GREEN)
        low_l.next_to(ax.c2p(3.5, 41), RIGHT, buff=0.15)
        high_l = layout.label("detached\ninner Oort cloud", font_size=16, color=P.GREEN)
        high_l.next_to(ax.c2p(3.5, 77), RIGHT, buff=0.15)
        cap2 = layout.caption("Below it, orbits Neptune still touches; above it, "
                              "orbits cut loose. Between: almost nothing", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        self.play(FadeIn(band), FadeIn(band_l), FadeIn(low_l), FadeIn(high_l))
        timing.hold_to_read(self, cap2, settle=0.8)
        self.play(FadeOut(VGroup(ax, xt, xl, yl, band, band_l, low_l, high_l, dots, cap2)))

        # 2. could one smooth population do that?
        edges = np.array(pe["edges"])
        counts = np.array(pe["counts"], dtype=float)
        expect = np.array(pe["null_expected"])
        ax2, lab2 = widgets.labeled_axes(
            [30, 90, 10], [0, 8, 2], x_label="perihelion q (AU)",
            y_label="objects per 5 AU", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=0)
        VGroup(ax2, lab2).move_to([-2.3, 0.25, 0])
        band2 = lane(ax2, 0, 8, g_lo, g_hi)
        bars = widgets.histogram(ax2, edges, counts, color=P.GREEN, opacity=0.7)
        xs = np.repeat(edges, 2)[1:-1]
        ys = np.repeat(expect, 2)
        smooth = widgets.curve(ax2, xs, ys, color=P.TEAL, stroke_width=3)
        exp_gap = float(sum(e for e, lo in zip(expect, edges[:-1]) if g_lo <= lo < g_hi))
        key = VGroup(
            layout.label(f"observed: {pe['n_in_gap']} in the gap", font_size=18,
                         color=P.GREEN),
            layout.label(f"one smooth population: ~{exp_gap:.1f}", font_size=18,
                         color=P.TEAL),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([4.4, 1.9, 0])
        self.play(Create(ax2), FadeIn(lab2), FadeIn(band2))
        self.play(FadeIn(bars, lag_ratio=0.1), FadeIn(key[0]), run_time=1.2)
        cap3 = layout.caption("Fit one smooth, falling distribution to the same objects",
                              font_size=22)
        self.play(FadeIn(cap3))
        self.play(Create(smooth), FadeIn(key[1]), run_time=1.4)
        timing.hold_to_read(self, cap3, key, settle=0.4)
        odds = paper.result_readout("a smooth population leaves the gap this empty",
                                    f"{100 * d['p_paper_epoch']:.0f}% of the time",
                                    color=P.RED).scale(0.8)
        odds.move_to([4.1, -0.2, 0])
        cap4 = layout.caption("Paper: two separate populations, very unlikely to be one",
                              font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        self.play(FadeIn(odds))
        timing.hold_to_read(self, cap4, odds, settle=0.8)

        # 3. the paper's planet prediction, then a new object in the lane
        pred = VGroup(
            layout.label("With a distant planet, orbits cycle", font_size=17),
            layout.label("through the gap fast and linger above:", font_size=17),
            layout.label(f"gap ≈ {100 * d['published_gap_relative_abundance']:.0f}% as full "
                         "as q = 65–100 AU", font_size=17, color=P.TEAL),
            layout.label("(paper's prediction)", font_size=14, color=P.MUTED),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT).move_to([4.4, -0.4, 0])
        self.play(FadeOut(odds), FadeIn(pred))
        timing.hold_to_read(self, pred, settle=0.6)
        nb = new[0]
        today = d["today"]
        bars_now = widgets.histogram(ax2, edges, np.array(today["counts"], dtype=float),
                                     color=P.GREEN, opacity=0.7)
        drop = Dot(ax2.c2p(nb["q_au"], 7.5), radius=0.09, color=P.ORANGE)
        drop_l = layout.label(f"{nb['name']}: q = {nb['q_au']:.0f} AU, found after the paper",
                              font_size=16, color=P.ORANGE).next_to(drop, UP, buff=0.1)
        cap5 = layout.caption(f"Now the gap holds {today['n_in_gap']}: a smooth population "
                              f"odds {100 * d['p_today']:.0f}%. Thin, not empty",
                              font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), FadeIn(drop), FadeIn(drop_l), run_time=0.8)
        self.play(drop.animate.move_to(ax2.c2p(nb["q_au"], 0.5)), run_time=1.0)
        seen_now = layout.label(f"observed today: {today['n_in_gap']} in the gap", font_size=18,
                                color=P.GREEN).move_to(key[0], aligned_edge=LEFT)
        self.play(Transform(bars, bars_now), FadeOut(drop), Transform(key[0], seen_now))
        timing.hold_to_read(self, cap5, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "The 50–65 AU lane is thin, as a distant planet would keep it.")
