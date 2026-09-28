"""Belyakov & Batygin (2025) -- perturbation theory of scattered-disk stability.

Batygin et al. (2021) explained the scattered disk's chaos with one chain of
2:j resonances with Neptune, the leading (quadrupole) term of Neptune's
potential. Closer to Neptune that chain alone is too sparse to overlap.
Carrying the expansion to octupole order and beyond adds 1:j, 3:j, 4:j...
chains whose islands interleave with the 2:j ones; their mutual overlap is
what makes the nearer orbits chaotic. Reproduced in p9-2025-perturbation: the
resonance chains, their overlap and the onset perihelia are the crate's own
(anim.json -> papers -> p9-2025-perturbation); the paper's boundary fit is
drawn as a labelled comparison.
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
    RoundedRectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-perturbation"


class Perturbation2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        lo, hi = d["comb_window_au"]
        comb = min(d["combs"], key=lambda c: c["q_au"])
        self.add(paper.scene_header(CRATE))

        # 1. the resonance chains near a = 200 AU
        x0, x1 = -5.0, 3.6

        def xa(a):
            return x0 + (x1 - x0) * (a - lo) / (hi - lo)

        axis = VGroup(Line([x0, -2.0, 0], [x1, -2.0, 0], color=P.MUTED, stroke_width=2))
        for a in range(int(lo), int(hi) + 1, 5):
            axis.add(Line([xa(a), -2.07, 0], [xa(a), -1.93, 0], color=P.MUTED, stroke_width=2))
            axis.add(layout.label(str(a), font_size=15, color=P.MUTED).move_to([xa(a), -2.3, 0]))
        axis.add(layout.label("semi-major axis a (AU)", font_size=17).move_to([(x0 + x1) / 2, -2.7, 0]))

        def box(r, y, col, h=0.42):
            a_l = max(r["a_au"] - r["width_au"], lo)
            a_r = min(r["a_au"] + r["width_au"], hi)
            b = RoundedRectangle(width=max(xa(a_r) - xa(a_l), 0.05), height=h, corner_radius=0.1,
                                 color=col, stroke_width=1.3)
            return b.set_fill(col, opacity=0.3).move_to([(xa(a_l) + xa(a_r)) / 2, y, 0])

        ys = {1: 2.3, 2: 1.55, 3: 0.8, 4: 0.05}
        rows, tags = {}, {}
        for m, y in ys.items():
            col = P.ORANGE if m == 2 else P.FG
            rows[m] = VGroup(*[box(r, y, col) for r in comb["resonances"] if r["m"] == m])
            tags[m] = layout.label(f"{m}:j", font_size=18, color=col).next_to([x0, y, 0], LEFT,
                                                                             buff=0.35)
        union_y = -1.1
        union = VGroup(*[box(r, union_y, P.ORANGE if r["m"] == 2 else P.FG, h=0.5)
                         for r in comb["resonances"]])
        union_t = layout.label("all\ntogether", font_size=16).next_to([x0, union_y, 0], LEFT,
                                                                   buff=0.35)
        q = comb["q_au"]
        self.play(FadeIn(axis))
        cap = layout.caption(f"An orbit with perihelion {q:.0f} AU near a = 200 AU. The 2021 "
                             "theory: Neptune's 2:j resonances only", font_size=22)
        self.play(FadeIn(cap), FadeIn(tags[2]), Create(rows[2], lag_ratio=0.2), run_time=1.4)
        k2 = VGroup(
            layout.label("islands apart", font_size=17, color=P.ORANGE),
            layout.label(f"overlap K = {comb['overlap_2j_only']:.2f}", font_size=17,
                         color=P.ORANGE),
        ).arrange(DOWN, buff=0.06, aligned_edge=LEFT).next_to([x1, ys[2], 0], RIGHT, buff=0.3)
        self.play(FadeIn(k2))
        timing.hold_to_read(self, cap, k2, settle=0.4)
        cap2 = layout.caption("Expand Neptune's pull further: octupole 1:j and 3:j chains, "
                              "then 4:j, fill the gaps", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        for m in (1, 3, 4):
            self.play(FadeIn(tags[m]), Create(rows[m], lag_ratio=0.15), run_time=0.9)
        timing.hold_to_read(self, cap2, settle=0.3)
        k_all = VGroup(
            layout.label("islands overlap", font_size=17, color=P.RED),
            layout.label(f"K = {comb['overlap_all']:.2f} > 1: chaos", font_size=17,
                         color=P.RED),
        ).arrange(DOWN, buff=0.06, aligned_edge=LEFT).next_to([x1, union_y, 0], RIGHT, buff=0.3)
        cap3 = layout.caption("Stacked together the islands overlap: this orbit is chaotic "
                              "after all", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), run_time=0.6)
        self.play(FadeIn(union_t), FadeIn(union, lag_ratio=0.05), run_time=1.2)
        self.play(FadeIn(k_all))
        timing.hold_to_read(self, cap3, k_all, settle=0.8)
        self.play(FadeOut(VGroup(axis, *rows.values(), *tags.values(), union, union_t, k2, k_all,
                                 cap3)))

        # 2. where chaos begins, chain by chain
        on = d["onset"]
        a_on = np.array([o["a_au"] for o in on])
        ax, labels = widgets.labeled_axes(
            [100, 450, 50], [20, 60, 10], x_label="semi-major axis a (AU)",
            y_label="perihelion where chaos begins (AU)", y_rotate=True, numbers=True,
            x_length=8.4, y_length=4.2, shift_down=0)
        VGroup(ax, labels).move_to([-1.3, 0.25, 0])

        def series(key):
            pts = [(o["a_au"], o[key]) for o in on if o[key] is not None]
            return np.array([p[0] for p in pts]), np.array([p[1] for p in pts])

        a2, q2 = series("q_quadrupole_chain_au")
        aa, qa = series("q_all_chains_au")
        ap, qp = series("q_published_fit_au")
        c2 = widgets.curve(ax, a2, q2, color=P.ORANGE, stroke_width=3)
        ca = widgets.curve(ax, aa, qa, color=P.FG, stroke_width=3.5)
        cp = VGroup(*[DashedLine(ax.c2p(ap[k], qp[k]), ax.c2p(ap[k + 1], qp[k + 1]),
                                 color=P.TEAL, stroke_width=2.5) for k in range(len(ap) - 1)])
        key = VGroup(
            layout.label("2:j chain alone", font_size=17, color=P.ORANGE),
            layout.label("all chains (reproduced)", font_size=17, color=P.FG),
            layout.label("paper's boundary fit", font_size=17, color=P.TEAL),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([4.9, 1.6, 0])
        self.play(Create(ax), FadeIn(labels))
        cap4 = layout.caption(f"The 2:j chain alone overlaps only beyond ~{a2[0]:.0f} AU",
                              font_size=22)
        self.play(FadeIn(cap4), Create(c2), FadeIn(key[0]), run_time=1.2)
        timing.hold_to_read(self, cap4, settle=0.3)
        cap5 = layout.caption(f"With every chain, chaos reaches in to a = {aa[0]:.0f} AU "
                              "and sits close to the paper's fit", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), run_time=0.6)
        self.play(Create(ca), FadeIn(key[1]), run_time=1.4)
        self.play(Create(cp), FadeIn(key[2]), run_time=1.0)
        a_h = d["a_headline_au"]
        q_h, q_pub = d["q_onset_headline_au"], d["q_published_headline_au"]
        dot = Dot(ax.c2p(a_h, q_h), radius=0.09, color=P.FG).set_z_index(5)
        lab = layout.label(f"a = {a_h:.0f} AU: {q_h:.1f} AU (paper {q_pub:.1f})", font_size=17)
        lab.next_to(dot, DOWN + RIGHT, buff=0.12)
        self.play(FadeIn(dot), FadeIn(lab))
        timing.hold_to_read(self, cap5, lab, settle=0.8)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Near Neptune, chaos comes from many resonance chains crossing, not one.")
