"""Bailey, Brown & Batygin (2018) -- can resonances locate Planet Nine?

Earlier work read Planet Nine's semimajor axis off the observed objects by
assuming each sits in a simple N/1 or N/2 resonance with it. Bailey et al.
simulate a scattered disk under eccentric Planet Nines and find most resonant
objects occupy high-order ratios instead; once every ratio is allowed, the
implied semimajor axis spreads into a plateau. The resonance forest, the
reduced-scale N-body census and both implied-a9 distributions shown here are the
crate's (anim.json -> papers -> p9-2018-resonance).
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
    Scene,
    Transform,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2018-resonance"

A_LO, A_HI = 140.0, 610.0   # semimajor-axis window of the resonance strip (AU)
X_LO, X_HI = -6.1, 6.1      # its screen extent
STRIP_Y = -0.2              # strip baseline


def sx(a):
    return X_LO + (a - A_LO) / (A_HI - A_LO) * (X_HI - X_LO)


def tick(a, h, color, width, opacity=1.0):
    return Line([sx(a), STRIP_Y, 0], [sx(a), STRIP_Y + h, 0], color=color,
                stroke_width=width).set_opacity(opacity)


class Resonance2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        a_p9 = d["a_p9"]
        full = [r for r in d["catalog_full"] if A_LO <= r["a_au"] <= A_HI]
        simple = [r for r in full if r["simple"]]
        other = [r for r in full if not r["simple"]]

        self.add(paper.scene_header(CRATE))

        # 1. the resonance forest: simple ratios are sparse, all ratios are everywhere
        base = Line([X_LO, STRIP_Y, 0], [X_HI, STRIP_Y, 0], color=P.MUTED, stroke_width=2)
        nums = VGroup(*[
            layout.label(f"{a}", font_size=14, color=P.FG).next_to([sx(a), STRIP_Y, 0], DOWN,
                                                                   buff=0.12)
            for a in range(200, 501, 100)
        ])
        a_lab = layout.label("semimajor axis of a scattered-disk object (AU)", font_size=16)
        a_lab.next_to(nums, DOWN, buff=0.15)
        p9 = Dot([sx(a_p9), STRIP_Y, 0], radius=0.1, color=P.BLUE).set_z_index(3)
        p9_lab = layout.label(f"Planet Nine\na = {a_p9:.0f} AU", font_size=14, color=P.BLUE)
        p9_lab.next_to(p9, DOWN, buff=0.15)
        self.play(Create(base), FadeIn(nums), FadeIn(a_lab), FadeIn(p9), FadeIn(p9_lab))

        s_ticks = VGroup(*[tick(r["a_au"], 1.1, P.ORANGE, 3) for r in simple])
        s_names = VGroup(*[
            layout.label(f"{r['p']}:{r['q']}", font_size=14, color=P.ORANGE)
            .next_to([sx(r["a_au"]), STRIP_Y + 1.1, 0], UP, buff=0.08)
            for r in simple if (r["p"], r["q"]) in ((2, 1), (3, 1), (5, 2), (3, 2), (4, 1))
        ])
        cap = layout.caption(
            f"Earlier fits assumed N/1 or N/2 resonances: {len(simple)} of them in this range",
            font_size=22)
        self.play(LaggedStart(*[Create(t) for t in s_ticks], lag_ratio=0.08), FadeIn(s_names),
                  FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.6)

        o_ticks = VGroup(*[tick(r["a_au"], 0.7, P.FG, 1.2, opacity=0.45) for r in other])
        cap2 = layout.caption(
            f"Allow every ratio p:q with p ≤ 35, q ≤ 20: {len(full)} of them, one every "
            f"~{(A_HI - A_LO) / len(full):.1f} AU", font_size=22)
        self.play(LaggedStart(*[Create(t) for t in o_ticks], lag_ratio=0.01), FadeOut(cap),
                  FadeIn(cap2), run_time=2.2)
        timing.hold_to_read(self, cap2, settle=0.8)

        # 2. which resonances simulated objects actually occupy
        res = sorted(d["resonant"], key=lambda r: r["a_res_au"])
        stack = {}
        pdots = VGroup()
        for r in res:
            col = P.ORANGE if r["simple"] else P.TEAL
            key = round(sx(r["a_res_au"]) / 0.18)
            level = stack.get(key, 0)
            stack[key] = level + 1
            pdots.add(Dot([sx(r["a_res_au"]), STRIP_Y + 1.75 + 0.2 * level, 0], radius=0.075,
                          color=col))
        n_res, n_s, n_tot = d["n_resonant"], d["n_simple"], d["n_total"]
        legend = VGroup(
            layout.label(f"in N/1 or N/2: {n_s}", font_size=16, color=P.ORANGE),
            layout.label(f"in another ratio: {n_res - n_s}", font_size=16, color=P.TEAL),
        ).arrange(RIGHT, buff=0.6)
        legend.move_to([0, 2.75, 0])
        e9s = [row["e9"] for row in d["per_e9"]]
        cap3 = layout.caption(
            f"N-body: {n_tot} particles, e₉ = {min(e9s):.1f}-{max(e9s):.1f}, {d['t_kyr']:.0f} kyr: "
            f"{n_res} hold a resonant angle", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3),
                  LaggedStart(*[FadeIn(p, shift=DOWN * 0.3) for p in pdots], lag_ratio=0.06),
                  run_time=2.0)
        self.play(FadeIn(legend))
        timing.hold_to_read(self, cap3, legend, settle=0.8)

        eq = layout.equation(
            rf"P(\text{{all 6 in }} N/1,\,N/2) = ({n_s}/{n_res})^6 = "
            rf"{self._sci(d['p_all6'])}", color=P.FG, scale=0.8)
        eq.move_to([0, -2.1, 0])
        note = layout.label("reduced-scale run (paper: 4 Gyr); paper's bound: below 5%",
                            font_size=16, color=P.FG)
        note.next_to(eq, DOWN, buff=0.15)
        cap4 = layout.caption("So the six clustered objects are unlikely to all sit in simple "
                              "resonances", font_size=22)
        self.play(FadeOut(a_lab), FadeOut(cap3), FadeIn(cap4), FadeIn(eq), FadeIn(note))
        timing.hold_to_read(self, cap4, eq, note, settle=1.0)
        stage = VGroup(base, nums, p9, p9_lab, s_ticks, s_names, o_ticks, pdots, legend, eq, note)
        self.play(FadeOut(stage), FadeOut(cap4))

        # 3. what that does to the semimajor axis inferred from the six objects
        a9 = np.array(d["a9_dist"]["a9_au"])
        f5 = np.array(d["a9_dist"]["f5"]) * 1e3
        fu = np.array(d["a9_dist"]["full"]) * 1e3
        top = float(np.ceil(f5.max() * 1.25))
        ax, labels = widgets.labeled_axes(
            [300, 1000, 100], [0, top, top / 4], x_label="Planet Nine semimajor axis a₉ (AU)",
            y_label="relative likelihood", y_rotate=True, numbers=False, x_length=10.0,
            y_length=4.0, shift_down=-0.35)
        ax.get_x_axis().add_numbers(range(300, 1001, 100), font_size=16)
        labels[0].shift(DOWN * 0.3)
        c_f5 = widgets.curve(ax, a9, f5, color=P.ORANGE)
        c_fu = widgets.curve(ax, a9, fu, color=P.TEAL)
        peak = widgets.marker_line(ax, d["f5_peak_a9"], (0, top),
                                   f"p, q ≤ 5 only: peak at {d['f5_peak_a9']:.0f} AU",
                                   color=P.ORANGE, side=RIGHT)
        cap5 = layout.caption("Read a₉ off the six objects allowing only ratios with p, q ≤ 5: "
                              "a sharp peak", font_size=22)
        kepler = layout.equation(r"a_9 = a_{\rm obj}\,(p/q)^{2/3}", color=P.FG, scale=0.75)
        kepler.move_to(ax.c2p(430, 0.86 * top))
        self.play(Create(ax), FadeIn(labels), FadeIn(cap5), FadeIn(kepler))
        self.play(Create(c_f5), run_time=1.5)
        self.play(Create(peak))
        timing.hold_to_read(self, cap5, settle=0.6)

        ghost = c_f5.copy().set_stroke(opacity=0.35)
        self.add(ghost)
        cap6 = layout.caption("Count every ratio and the peak dissolves into a plateau",
                              font_size=22)
        self.play(Transform(c_f5, c_fu), FadeOut(cap5), FadeIn(cap6), run_time=2.0)
        ratios = VGroup(
            layout.label(f"peak / mean, p, q ≤ 5: {d['f5_peak_to_mean']:.2f}",
                         font_size=16, color=P.ORANGE),
            layout.label(f"peak / mean, every ratio: {d['full_peak_to_mean']:.2f}",
                         font_size=16, color=P.TEAL),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        ratios.move_to(ax.c2p(870, 0.8 * top))
        self.play(FadeIn(ratios))
        timing.hold_to_read(self, cap6, ratios, settle=1.2)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "Resonances are too densely packed to tell us where Planet Nine is.")

    @staticmethod
    def _sci(x):
        m, e = f"{x:.1e}".split("e")
        return rf"{m}\times 10^{{{int(e)}}}"
