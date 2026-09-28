"""Brown & Batygin (2019) -- orbital clustering in the distant solar system.

Fourteen distant orbits, and a bias model that covers where the perihelia point
and how the orbital planes tilt at the same time. Each orbit becomes two small
vectors (Poincare variables): one for its perihelion direction, one for its
orbital pole. Average them over the sample: random orientations cancel, aligned
ones do not. The null is thousands of fake samples drawn where the surveys
could have found them. The vectors, the null clouds, the probabilities and the
detection threshold are the crate's own (anim.json -> papers ->
p9-2019-clustering).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arrow,
    Circle,
    Create,
    DashedLine,
    DashedVMobject,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Polygon,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2019-clustering"

# Published in the paper's abstract; drawn only as a labelled comparison.
PAPER_P_COMBINED = 0.002


class Plane:
    """A square panel for one pair of Poincare variables, origin at the centre."""

    def __init__(self, centre, half_width, extent, title, axes_names):
        self.c = np.array(centre, dtype=float)
        self.k = half_width / extent
        self.extent = extent
        box = Polygon(*[self.c + half_width * np.array(v) for v in
                        ([-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0])],
                      color=P.MUTED, stroke_width=1.2)
        box.set_fill("#16171f", opacity=1.0)
        cross = VGroup(
            Line(self.c + LEFT * half_width, self.c + RIGHT * half_width, color=P.MUTED,
                 stroke_width=0.8).set_stroke(opacity=0.5),
            Line(self.c + DOWN * half_width, self.c + UP * half_width, color=P.MUTED,
                 stroke_width=0.8).set_stroke(opacity=0.5))
        head = layout.label(title, font_size=19, color=P.FG)
        head.next_to(box, UP, buff=0.12)
        xn = layout.label(axes_names[0], font_size=15, color=P.MUTED)
        xn.next_to(self.c + RIGHT * half_width, DOWN, buff=0.08).shift(LEFT * 0.2)
        yn = layout.label(axes_names[1], font_size=15, color=P.MUTED)
        yn.next_to(self.c + UP * half_width, RIGHT, buff=0.08).shift(DOWN * 0.15)
        self.frame = VGroup(box, cross, head, xn, yn)

    def p(self, u, v):
        return self.c + self.k * np.array([u, v, 0.0])

    def arrow(self, u, v, color, width=2.5, opacity=1.0):
        a = Arrow(self.c, self.p(u, v), buff=0, color=color, stroke_width=width,
                  max_tip_length_to_length_ratio=0.12, max_stroke_width_to_length_ratio=8)
        return a.set_opacity(opacity)

    def ring(self, radius, color):
        c = Circle(radius=radius * self.k, color=color, stroke_width=2).move_to(self.c)
        return DashedVMobject(c, num_dashes=48)


class Clustering2019(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        mean = d["mean"]
        null = d["null_means"]

        self.add(paper.scene_header(CRATE))

        # 1. every orbit as two arrows, and their average
        left = Plane((-3.4, 0.0, 0), 2.45, 1.3, "perihelion direction",
                     ("toward ϖ = 0°", "ϖ = 90°"))
        right = Plane((3.4, 0.0, 0), 2.45, 0.4, "tilt of the orbital plane",
                      ("toward Ω = 0°", "Ω = 90°"))
        old = [o for o in objs if o["in_2017"]]
        new = [o for o in objs if not o["in_2017"]]
        arrows_l = VGroup(*[left.arrow(o["x"], o["y"], P.GREEN, 2.2, 0.55) for o in old],
                          *[left.arrow(o["x"], o["y"], P.GREEN, 3.2) for o in new])
        arrows_r = VGroup(*[right.arrow(o["p"], o["q_var"], P.GREEN, 2.2, 0.55) for o in old],
                          *[right.arrow(o["p"], o["q_var"], P.GREEN, 3.2) for o in new])
        cap = layout.caption(
            f"{len(objs)} orbits ({len(new)} new since 2017), each an arrow: "
            f"longer when more eccentric or more tilted", font_size=22)
        self.play(FadeIn(left.frame), FadeIn(right.frame))
        self.play(LaggedStart(*[Create(a) for a in arrows_l], lag_ratio=0.08),
                  LaggedStart(*[Create(a) for a in arrows_r], lag_ratio=0.08),
                  FadeIn(cap), run_time=2.6)
        timing.hold_to_read(self, cap, settle=0.8)

        m_l = Dot(left.p(mean["x"], mean["y"]), radius=0.09, color=P.GREEN).set_z_index(4)
        m_r = Dot(right.p(mean["p"], mean["q_var"]), radius=0.09, color=P.GREEN).set_z_index(4)
        ring_l = left.ring(d["observed_perihelion"], P.GREEN)
        ring_r = right.ring(d["observed_pole"], P.GREEN)
        cap2 = layout.caption(
            "Average the arrows: random directions cancel, these leave a net pull",
            font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(arrows_l.animate.set_opacity(0.18), arrows_r.animate.set_opacity(0.18),
                  FadeIn(m_l, scale=2), FadeIn(m_r, scale=2), FadeIn(cap2), run_time=1.2)
        self.play(Create(ring_l), Create(ring_r), run_time=1.0)
        timing.hold_to_read(self, cap2, settle=0.8)

        # 2. the same average for thousands of fake samples the surveys could have found
        cloud_l = VGroup(*[Dot(left.p(s["x"], s["y"]), radius=0.018, color=P.FG)
                           .set_opacity(0.55) for s in null])
        cloud_r = VGroup(*[Dot(right.p(s["p"], s["q_var"]), radius=0.018, color=P.FG)
                           .set_opacity(0.55) for s in null])
        cap3 = layout.caption(
            f"{len(null):,} fake samples of {d['n_sample']}, placed where the surveys "
            f"could have seen them", font_size=22)
        self.play(FadeOut(cap2), FadeOut(arrows_l), FadeOut(arrows_r), run_time=0.5)
        self.play(FadeIn(cloud_l, lag_ratio=0.002), FadeIn(cloud_r, lag_ratio=0.002),
                  FadeIn(cap3), run_time=2.4)
        self.add(m_l, m_r)
        timing.hold_to_read(self, cap3, settle=0.8)

        row_l = layout.label(f"as far out as the real sample:  {100 * d['p_perihelion']:.1f}%",
                             font_size=16, color=P.FG)
        row_l.next_to(left.frame[0], DOWN, buff=0.12)
        row_r = layout.label(f"as far out as the real sample:  {100 * d['p_pole']:.1f}%",
                             font_size=16, color=P.FG)
        row_r.next_to(right.frame[0], DOWN, buff=0.12)
        cap4 = layout.caption(
            f"Both at once, under the bias:  {100 * d['p_combined']:.2f}%   "
            f"(paper: {100 * PAPER_P_COMBINED:.1f}%)", font_size=24)
        self.play(FadeOut(cap3), run_time=0.4)
        self.play(FadeIn(row_l), FadeIn(row_r))
        timing.hold_to_read(self, row_l, settle=0.6)
        self.play(FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=1.2)
        self.play(FadeOut(VGroup(left.frame, right.frame, cloud_l, cloud_r, m_l, m_r, ring_l,
                                 ring_r, row_l, row_r, cap4)))

        # 3. why OSSOS on its own could not tell
        th = d["threshold"]
        n_max = th["n"][-1]
        ax, labels = widgets.labeled_axes(
            [0, n_max, 5], [0, 1.0, 0.2], x_label="number of distant orbits in the sample",
            y_label="alignment of the perihelia", y_rotate=True, numbers=True,
            x_length=9.0, y_length=4.2, shift_down=-0.35)
        plot = VGroup(ax, labels).shift(LEFT * 1.4)
        top = [ax.c2p(n, r) for n, r in zip(th["n"], th["r_bar"]) if r <= 1.0]
        shade = Polygon(*top, ax.c2p(th["n"][-1], 1.0), ax.c2p(th["n"][0], 1.0), stroke_width=0)
        shade.set_fill(P.TEAL, opacity=0.12)
        line = widgets.curve(ax, [n for n, r in zip(th["n"], th["r_bar"]) if r <= 1.0],
                             [r for r in th["r_bar"] if r <= 1.0], color=P.TEAL)
        line_lab = layout.label("detectable at 95% above this line", font_size=15,
                                color=P.TEAL)
        line_lab.next_to(ax.c2p(18, 0.85), UP, buff=0.05)
        oss = d["ossos"]
        p_oss = Dot(ax.c2p(oss["n"], oss["r_bar_varpi"]), radius=0.09, color=P.GREEN)
        need = DashedLine(ax.c2p(oss["n"], oss["r_bar_varpi"]),
                          ax.c2p(oss["n"], oss["r_bar_needed"]), color=P.RED, stroke_width=2.5)
        oss_lab = VGroup(
            layout.label(f"OSSOS alone: {oss['n']} orbits", font_size=15, color=P.GREEN),
            layout.label(f"alignment {oss['r_bar_varpi']:.2f}, needs {oss['r_bar_needed']:.2f}",
                         font_size=14, color=P.RED),
        ).arrange(DOWN, buff=0.06, aligned_edge=LEFT)
        oss_lab.next_to(p_oss, RIGHT, buff=0.2).shift(UP * 0.1)
        p_all = Dot(ax.c2p(d["n_sample"], d["r_bar_varpi"]), radius=0.09, color=P.GREEN)
        all_lab = VGroup(
            layout.label(f"all {d['n_sample']} orbits", font_size=15, color=P.GREEN),
            layout.label(f"alignment {d['r_bar_varpi']:.2f}, needs {d['r_bar_needed']:.2f}",
                         font_size=14, color=P.FG),
        ).arrange(DOWN, buff=0.06, aligned_edge=LEFT)
        all_lab.next_to(p_all, DOWN + RIGHT, buff=0.12)
        cap5 = layout.caption(
            f"With {oss['n']} orbits, only an alignment above {oss['r_bar_needed']:.2f} "
            f"would stand out: OSSOS alone could not see it", font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(shade), Create(line), FadeIn(line_lab), run_time=1.4)
        self.play(FadeIn(p_oss), Create(need), FadeIn(oss_lab), FadeIn(cap5))
        timing.hold_to_read(self, cap5, oss_lab, settle=1.0)
        cap6 = layout.caption(
            "Even 14 perihelia alone sit near the line; adding the orbital planes decides it",
            font_size=22)
        self.play(FadeOut(cap5), run_time=0.4)
        self.play(FadeIn(p_all), FadeIn(all_lab), FadeIn(cap6))
        timing.hold_to_read(self, cap6, all_lab, settle=1.2)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, f"With both biases modelled, {d['n_sample']} orbits this aligned: "
                  f"a {100 * d['p_combined']:.1f}% chance.")
