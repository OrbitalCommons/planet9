"""Khain, Becker & Adams (2020) -- the resonance hopping effect.

In simulations, distant objects jump abruptly from one resonance with Planet
Nine to another. The paper traces the jumps to Neptune: its kicks at each
perihelion passage make the semi-major axis random-walk across the ladder of
Planet Nine resonances, and the anti-aligned objects survive whether or not
they sit in one. Reproduced in p9-2020-resonance-hopping: the resonance
ladder, the random walk (Neptune's diffusion coefficient) and the classified
synthetic belt are the crate's own (anim.json -> papers ->
p9-2020-resonance-hopping).
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
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2020-resonance-hopping"
CLASS_COLOR = {"resonant": P.BLUE, "hopping": P.ORANGE, "non_resonant": P.FG}
CLASS_NAME = {"resonant": "stuck in one resonance", "hopping": "hopping between resonances",
              "non_resonant": "in no resonance"}


class ResonanceHopping2020(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        res = d["simple_resonances"]
        a_lo, a_hi = d["a_range_au"]
        walk = d["walk"]
        a_walk = np.array(walk["a_au"])
        t_walk = np.arange(len(a_walk)) * walk["dt_myr"]

        self.add(paper.scene_header(CRATE))

        # 1. the ladder of Planet Nine resonances
        ax = widgets.axes([150, 500, 50], [0, 1, 1], x_length=10.6, y_length=1.2,
                          shift_down=0)
        ax.move_to([-0.6, 1.6, 0])
        ax.get_y_axis().set_opacity(0)
        ax.get_x_axis().add_numbers(font_size=16)
        xl = layout.label("semi-major axis a (AU)", font_size=18).next_to(ax, DOWN, buff=0.45)
        rungs = VGroup()
        tags = VGroup()
        for r in res:
            x = r["a_au"]
            rungs.add(Line(ax.c2p(x, 0), ax.c2p(x, 1), color=P.BLUE, stroke_width=3))
            if r["label"] in ("4:1", "3:1", "5:2", "2:1", "3:2"):
                tags.add(layout.label(r["label"], font_size=16, color=P.BLUE)
                         .next_to(ax.c2p(x, 1), UP, buff=0.08))
        p9 = layout.label(f"Planet Nine\n{d['p9']['a_au']:.0f} AU →", font_size=16,
                          color=P.BLUE).next_to(ax.c2p(500, 0.5), RIGHT, buff=0.1)
        self.play(Create(ax), FadeIn(xl))
        self.play(Create(rungs, lag_ratio=0.12), FadeIn(tags), FadeIn(p9), run_time=1.6)
        spacing = d["median_simple_spacing_au"]
        cap = layout.caption(
            f"Planet Nine's main resonances: rungs about {spacing:.0f} AU apart, "
            f"plus {d['n_spectrum'] - len(res)} weaker ones between", font_size=22)
        self.play(FadeIn(cap))
        timing.hold_to_read(self, cap, settle=0.6)

        # 2. Neptune's kicks make a random-walk
        t = ValueTracker(0.0)
        n = len(a_walk)

        def upto():
            return max(2, int(round(t.get_value() * (n - 1))) + 1)

        ax2 = widgets.axes([0, float(t_walk[-1]), 0.5], [280, 400, 40], x_length=8.2,
                           y_length=2.7, shift_down=0)
        ax2.move_to([0.3, -1.25, 0])
        ax2.add_coordinates()
        l2x = layout.label("time (Myr)", font_size=16).next_to(ax2, DOWN, buff=0.1)
        l2y = layout.label("a (AU)", font_size=16).next_to(ax2, LEFT, buff=0.12)
        guide = VGroup()
        for r in res:
            if 280 <= r["a_au"] <= 400:
                guide.add(DashedLine(ax2.c2p(0, r["a_au"]), ax2.c2p(t_walk[-1], r["a_au"]),
                                     color=P.BLUE, stroke_width=1.5).set_opacity(0.6))
                guide.add(layout.label(r["label"], font_size=14, color=P.BLUE)
                          .next_to(ax2.c2p(t_walk[-1], r["a_au"]), RIGHT, buff=0.08))
        trace = always_redraw(lambda: widgets.curve(
            ax2, t_walk[:upto()], a_walk[:upto()], color=P.GREEN, stroke_width=2.2))
        dot = always_redraw(lambda: Dot(ax.c2p(a_walk[upto() - 1], 0.5), radius=0.1,
                                        color=P.GREEN).set_z_index(5))
        cap2 = layout.caption(
            f"Each pass by Neptune (perihelion {walk['q_au']:.0f} AU) nudges a: "
            "a random walk up and down the ladder", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        self.play(Create(ax2), FadeIn(l2x), FadeIn(l2y),
                  FadeIn(guide))
        self.add(trace, dot)
        self.play(t.animate.set_value(1.0), run_time=9.0, rate_func=linear)
        cross = d["kick_crossing_time_myr"]
        cap3 = layout.caption(
            f"It drifts one rung ({spacing:.0f} AU) in ~{cross:.1f} Myr: "
            "the resonance 'hops' are Neptune's doing", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.8)
        trace.clear_updaters()
        dot.clear_updaters()
        self.play(FadeOut(VGroup(ax, xl, rungs, tags, p9, ax2, l2x, l2y, guide, trace, dot,
                                 cap3)))

        # 3. who is in a resonance at all
        bins = d["bins"]
        edges = [b["a_lo"] for b in bins] + [bins[-1]["a_hi"]]
        ax3, lab3 = widgets.labeled_axes(
            [150, 500, 50], [0, 1, 0.25], x_label="semi-major axis a (AU)",
            y_label="fraction of belt", y_rotate=True, numbers=True,
            x_length=7.6, y_length=4.2, shift_down=0)
        VGroup(ax3, lab3).move_to([-2.0, 0.2, 0])
        base = np.zeros(len(bins))
        stacks = VGroup()
        for cls in ("resonant", "hopping", "non_resonant"):
            vals = np.array([b[cls] for b in bins])
            stacks.add(widgets.histogram(ax3, edges, vals, color=CLASS_COLOR[cls],
                                         opacity=0.75 if cls != "non_resonant" else 0.35,
                                         base=base.copy()))
            base += vals
        key = VGroup()
        for cls in ("non_resonant", "hopping", "resonant"):
            key.add(VGroup(
                layout.label(f"{100 * d[cls]:.0f}%", font_size=24, color=CLASS_COLOR[cls],
                             weight="BOLD"),
                layout.label(CLASS_NAME[cls], font_size=17, color=CLASS_COLOR[cls]),
            ).arrange(RIGHT, buff=0.2))
        key.arrange(DOWN, buff=0.3, aligned_edge=LEFT).move_to([4.6, 0.6, 0])
        head = layout.label(f"{d['n_population']:,} synthetic distant orbits",
                            font_size=17, color=P.MUTED).next_to(key, UP, buff=0.35,
                                                                 aligned_edge=LEFT)
        cap4 = layout.caption("Hopping takes over toward Planet Nine, where "
                              "the rungs crowd and overlap", font_size=22)
        self.play(Create(ax3), FadeIn(lab3), FadeIn(cap4))
        self.play(FadeIn(stacks, lag_ratio=0.3), run_time=1.6)
        self.play(FadeIn(head), FadeIn(key, lag_ratio=0.2))
        timing.hold_to_read(self, cap4, key, settle=0.6)
        cap5 = layout.caption("Paper: anti-aligned objects survive without any "
                              "resonance, so the herding is secular", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), key[0].animate.scale(1.12, about_edge=LEFT))
        timing.hold_to_read(self, cap5, settle=0.8)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Neptune drives the hops; Planet Nine's alignment needs no resonance.")
