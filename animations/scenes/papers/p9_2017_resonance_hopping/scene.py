"""Becker et al. (2017) -- distant objects hop between Planet Nine resonances.

Becker et al. integrate eight distant TNOs under a Monte Carlo suite of Planet
Nines and find that some stay in one mean-motion resonance while others hop
between neighbouring resonances and still keep their anti-aligned orbits. The
reproduction crate (p9-2017-resonance-hopping) gives the analytic reason: each
resonance is an island of finite width in semimajor axis, the islands crowd and
widen toward Planet Nine, and past a computed semimajor axis neighbours overlap
(Chirikov K >= 1), so an object there cannot stay in one. Every width, overlap
parameter, threshold and classification drawn here is the crate's
(anim.json -> papers -> p9-2017-resonance-hopping).
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
    LaggedStart,
    Line,
    Rectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2017-resonance-hopping"

A_LO, A_HI = 250.0, 720.0   # semimajor-axis window shared by both panels (AU)
STRIP_Y = 2.05              # baseline of the resonance-island strip
BAND_H = 0.42


def island(ax, a, half_width, color, opacity=0.55, a_max=A_HI):
    """A resonance island drawn to scale: centre a, full width 2 * half_width,
    clipped at ``a_max`` (Planet Nine's own orbit)."""
    x0 = ax.c2p(max(a - half_width, A_LO), 0)[0]
    x1 = ax.c2p(min(a + half_width, a_max), 0)[0]
    r = Rectangle(width=max(x1 - x0, 0.02), height=BAND_H, stroke_width=0)
    r.set_fill(color, opacity=opacity)
    r.move_to([(x0 + x1) / 2, STRIP_Y + BAND_H / 2, 0])
    return r


def x_of(ax, a):
    return ax.c2p(a, 0)[0]


class ResonanceHopping2017(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        a9 = d["a9_au"]
        a_hop = d["a_hop"]
        prof_a = np.array(d["profile"]["a_mid_au"])
        prof_k = np.array(d["profile"]["k"])

        self.add(paper.scene_header(CRATE))

        # Shared semimajor-axis frame: overlap parameter K below, islands above.
        ax = widgets.axes([A_LO, A_HI, 50], [-1, 1.5, 0.5], x_length=12.0, y_length=2.7,
                          shift_down=0.0)
        ax.move_to([0.25, -0.75, 0])
        ax.get_x_axis().add_numbers(range(300, 701, 100), font_size=16)
        x_lab = layout.label("semimajor axis of the small body  (AU)", font_size=16)
        x_lab.next_to(ax, DOWN, buff=0.2)
        base = Line([x_of(ax, A_LO), STRIP_Y, 0], [x_of(ax, A_HI), STRIP_Y, 0],
                    color=P.MUTED, stroke_width=2)
        p9 = Dot([x_of(ax, a9), STRIP_Y, 0], radius=0.09, color=P.BLUE)
        p9_lab = layout.label(f"Planet Nine  a = {a9:.0f} AU", font_size=14, color=P.BLUE)
        p9_lab.next_to(p9, DOWN, buff=0.14).shift(LEFT * 0.55)
        self.play(Create(base), FadeIn(p9), FadeIn(p9_lab))
        p9.set_z_index(3)

        # 1. the familiar resonances are narrow, isolated islands
        simple = [r for r in d["simple"] if r["a_au"] >= A_LO + 15]
        isles = VGroup(*[island(ax, r["a_au"], r["half_width_au"], P.ORANGE) for r in simple])
        names = VGroup(*[
            layout.label(f"{r['p']}:{r['q']}", font_size=14, color=P.ORANGE)
            .next_to(isles[k], UP, buff=0.1)
            for k, r in enumerate(simple)
        ])
        cap = layout.caption("Each resonance is an island in semimajor axis, drawn here to scale",
                             font_size=22)
        self.play(LaggedStart(*[FadeIn(i) for i in isles], lag_ratio=0.15),
                  FadeIn(names), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.6)
        w21 = next(r for r in simple if (r["p"], r["q"]) == (2, 1))
        gaps = np.diff(sorted(r["a_au"] for r in simple))
        cap2 = layout.caption(
            f"The 2:1 is only ±{w21['half_width_au']:.1f} AU wide, with {gaps.min():.0f}-"
            f"{gaps.max():.0f} AU gaps between islands: an object can sit in one", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.6)

        # 2. toward Planet Nine the first-order islands crowd together and widen
        chain = [r for r in d["first_order"] if (r["p"], r["q"]) != (2, 1) and r["a_au"] <= A_HI]
        c_isles = VGroup(*[island(ax, r["a_au"], r["half_width_au"], P.ORANGE, opacity=0.45,
                                  a_max=a9) for r in chain])
        c_names = VGroup(*[
            layout.label(f"{r['p']}:{r['q']}", font_size=14, color=P.ORANGE)
            .next_to(c_isles[k], UP, buff=0.1)
            for k, r in enumerate(chain[:3])
        ])
        dots_lab = layout.label("… 40:39", font_size=14, color=P.ORANGE)
        dots_lab.next_to(c_isles[-1], UP, buff=0.1).shift(LEFT * 0.3)

        y_labels = VGroup(*[
            layout.label(t, font_size=14, color=P.FG).next_to(ax.c2p(A_LO, v), LEFT, buff=0.12)
            for t, v in (("0.1", -1), ("1", 0), ("10", 1))
        ])
        y_title = layout.label("overlap K", font_size=15).rotate(np.pi / 2)
        y_title.next_to(y_labels, LEFT, buff=0.14)
        keep = (prof_a >= A_LO) & (prof_k <= 10 ** 1.5)
        kcurve = widgets.curve(ax, prof_a[keep], np.log10(prof_k[keep]), color=P.ORANGE)
        k1 = DashedLine(ax.c2p(A_LO, 0), ax.c2p(A_HI, 0), color=P.MUTED, stroke_width=2)
        k1_lab = layout.label("K = 1: neighbouring islands touch", font_size=14, color=P.FG)
        k1_lab.next_to(ax.c2p(A_LO, 0), UP + RIGHT, buff=0.08).shift(RIGHT * 0.1)
        cap3 = layout.caption("Closer to Planet Nine the islands crowd together and grow wider",
                              font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), Create(ax), FadeIn(y_labels), FadeIn(y_title),
                  FadeIn(x_lab), Create(k1), FadeIn(k1_lab))
        self.play(LaggedStart(*[FadeIn(i) for i in c_isles], lag_ratio=0.08),
                  FadeIn(c_names), FadeIn(dots_lab), Create(kcurve), run_time=3.0)
        timing.hold_to_read(self, cap3, settle=0.4)

        # the overlap edge
        hop_line = DashedLine([x_of(ax, a_hop), ax.c2p(0, -1)[1], 0],
                              [x_of(ax, a_hop), STRIP_Y + BAND_H + 0.05, 0],
                              color=P.RED, stroke_width=2.5)
        zone = Rectangle(width=x_of(ax, a9) - x_of(ax, a_hop), height=BAND_H + 0.16,
                         stroke_width=0).set_fill(P.RED, opacity=0.14)
        zone.move_to([(x_of(ax, a_hop) + x_of(ax, a9)) / 2, STRIP_Y + BAND_H / 2, 0])
        hop_dot = Dot(ax.c2p(a_hop, np.log10(np.interp(a_hop, prof_a, prof_k))),
                      radius=0.07, color=P.RED)
        hop_lab = layout.label(f"overlap from {a_hop:.0f} AU:\nobjects hop", font_size=15,
                               color=P.RED)
        hop_lab.next_to(ax.c2p(a_hop, 0.9), LEFT, buff=0.25)
        n1_lab = layout.label(f"the whole n:1 chain stays below K = {d['simple_k_max']:.2f}",
                              font_size=14, color=P.FG)
        n1_lab.next_to(ax.c2p(A_LO, -1), UP + RIGHT, buff=0.12).shift(UP * 0.2)
        cap4 = layout.caption("Once neighbours overlap (K > 1), an object cannot stay in any one",
                              font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), Create(hop_line), FadeIn(zone), FadeIn(hop_dot),
                  FadeIn(hop_lab), FadeIn(n1_lab))
        timing.hold_to_read(self, cap4, hop_lab, settle=1.0)

        stage1 = VGroup(ax, y_labels, y_title, x_lab, k1, k1_lab, kcurve, hop_line, zone, hop_dot,
                        hop_lab, n1_lab, base, p9, p9_lab, isles, names, c_isles, c_names,
                        dots_lab)
        self.play(FadeOut(stage1), FadeOut(cap4))

        # 3. which of the observed distant objects reach the overlap zone
        rows = sorted(d["sample"], key=lambda s: s["big_q_au"])
        ax2 = widgets.axes([0, 1000, 100], [0, len(rows) + 1, 1], x_length=9.4, y_length=4.5,
                           shift_down=0.0)
        ax2.move_to([0.55, 0.72, 0])
        ax2.get_x_axis().add_numbers(range(0, 1001, 200), font_size=16)
        ax2.get_y_axis().set_opacity(0)
        x2_lab = layout.label("distance from the Sun (AU):  perihelion ── aphelion, dot = a",
                              font_size=15)
        x2_lab.next_to(ax2, DOWN, buff=0.38)
        z0, z1 = ax2.c2p(a_hop, 0), ax2.c2p(a9, len(rows))
        zone2 = Rectangle(width=z1[0] - z0[0], height=z1[1] - z0[1], stroke_width=0)
        zone2.set_fill(P.RED, opacity=0.14).move_to([(z0[0] + z1[0]) / 2, (z0[1] + z1[1]) / 2, 0])
        zone2_lab = layout.label(f"overlap zone {a_hop:.0f}-{a9:.0f} AU", font_size=14,
                                 color=P.RED)
        zone2_lab.next_to(zone2, UP, buff=0.1).align_to(zone2, RIGHT)
        p9_line = DashedLine(ax2.c2p(a9, 0), ax2.c2p(a9, len(rows)), color=P.BLUE,
                             stroke_width=2)
        p9_lab2 = layout.label("Planet Nine's a", font_size=14, color=P.BLUE)
        p9_lab2.next_to(p9_line.get_end(), UP + RIGHT, buff=0.1)

        verdict = {"Sedna": "stays", "2012 VP113": "stays", "2007 TG422": "hops",
                   "2013 RF98": "hops"}
        bars = VGroup()
        tags = VGroup()
        for k, s in enumerate(rows):
            y = k + 0.5
            hop = s["state"] == "Hopping"
            col = P.ORANGE if hop else P.GREEN
            ln = Line(ax2.c2p(s["q_au"], y), ax2.c2p(s["big_q_au"], y), color=col, stroke_width=4)
            dot = Dot(ax2.c2p(s["a_au"], y), radius=0.06, color=col)
            name = layout.label(s["name"], font_size=14, color=P.FG)
            name.next_to(ax2.c2p(0, y), LEFT, buff=0.15)
            bars.add(VGroup(name, ln, dot))
            if s["name"] in verdict:
                agree = verdict[s["name"]] == ("hops" if hop else "stays")
                t = layout.label(f"paper: {verdict[s['name']]}", font_size=13,
                                 color=P.FG if agree else P.RED)
                t.next_to(ax2.c2p(s["big_q_au"], y), RIGHT, buff=0.15)
                tags.add(t)
        key = VGroup(
            layout.label("reaches the zone: hops", font_size=14, color=P.ORANGE),
            layout.label("stays clear: locked", font_size=14, color=P.GREEN),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        key.next_to(ax2.c2p(0, len(rows) + 0.5), RIGHT, buff=0.1)
        cap5 = layout.caption(
            f"Crate's test on the observed objects: {d['n_hopping']} of {d['n_sample']} "
            "reach the overlap zone at aphelion", font_size=22)
        self.play(Create(ax2), FadeIn(x2_lab), FadeIn(zone2), FadeIn(zone2_lab),
                  Create(p9_line), FadeIn(p9_lab2), FadeIn(cap5))
        self.play(LaggedStart(*[FadeIn(b) for b in bars], lag_ratio=0.12), FadeIn(key),
                  run_time=2.0)
        timing.hold_to_read(self, cap5, settle=0.6)
        cap6 = layout.caption(
            "Becker's hoppers TG422 and RF98 match; Sedna, which the paper finds stays put, does not",
            font_size=20)
        self.play(FadeOut(cap5), FadeIn(cap6), FadeIn(tags))
        timing.hold_to_read(self, cap6, settle=1.2)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, f"Past ~{a_hop:.0f} AU Planet Nine's resonances overlap: objects hop, not lock.")
