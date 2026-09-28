"""Napier et al. (2021) -- no evidence for orbital clustering in the ETNOs.

The distant objects were found by surveys that looked in particular directions
at particular times of year. Napier et al. fold each survey's selection into the
null and find the observed alignment unremarkable. Reproduced in
p9-2021-napier-critique with a stand-in selection function for the surveys'
pointing histories; the sample, the selection curve, both null distributions and
every p-value shown are the crate's own
(anim.json -> papers -> p9-2021-napier-critique).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Arrow,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Rotate,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2021-napier-critique"


class Dial(VGroup):
    """A compass of orbital longitude: 0° to the right, increasing anticlockwise."""

    def __init__(self, radius=2.0, centre=(0.0, 0.0, 0.0), title=None):
        super().__init__()
        self.radius = radius
        self.centre = np.array(centre, dtype=float)
        self.add(Circle(radius=radius, color=P.MUTED, stroke_width=1.5).move_to(self.centre))
        for deg in range(0, 360, 30):
            a, b = self.p(deg, radius), self.p(deg, radius + (0.12 if deg % 90 == 0 else 0.06))
            self.add(Line(a, b, color=P.MUTED, stroke_width=1.2))
        for deg in (0, 90, 180, 270):
            lab = layout.label(f"{deg}°", font_size=13, color=P.MUTED)
            self.add(lab.move_to(self.p(deg, radius + 0.38)))
        if title:
            t = layout.label(title, font_size=15, color=P.FG)
            self.add(t.move_to(self.centre + DOWN * (radius + 0.85)))

    def p(self, deg, r=None):
        r = self.radius if r is None else r
        t = np.deg2rad(deg)
        return self.centre + r * np.array([np.cos(t), np.sin(t), 0.0])

    def dots(self, degs, color=P.GREEN, radius=0.075):
        return VGroup(*[Dot(self.p(d), radius=radius, color=color).set_z_index(3) for d in degs])

    def polar(self, degs, values, vmax, r_max, color, opacity=0.22):
        """A closed polar curve r = r_max * value / vmax, filled."""
        m = VMobject(color=color, stroke_width=2.2)
        m.set_points_as_corners([self.p(d, r_max * v / vmax) for d, v in zip(degs, values)])
        return m.set_fill(color, opacity=opacity)


def pct(p):
    return f"{100 * p:.1f}%" if p < 0.1 else f"{100 * p:.0f}%"


class NapierCritique2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        sel = d["selection"]
        nul = d["null_r_bar"]
        varpis = [o["varpi_deg"] for o in d["objects"]]
        n = d["n_sample"]

        self.add(paper.scene_header(CRATE))

        # 1. the alignment, taken at face value
        dial = Dial(radius=2.0, centre=(-3.3, 0.5, 0.0), title="longitude of perihelion ϖ")
        seen = dial.dots(varpis)
        mean = Arrow(dial.centre, dial.p(d["mean_varpi_deg"], dial.radius * d["r_bar"]),
                     buff=0, color=P.GREEN, stroke_width=4, max_tip_length_to_length_ratio=0.18)
        cap = layout.caption(f"{n} distant objects (the paper used {d['paper_n_sample']}): "
                             "perihelia bunched on one side", font_size=22)
        naive = VGroup(
            layout.label(f"alignment strength  R = {d['r_bar']:.2f}", font_size=20, color=P.GREEN),
            layout.label(f"if any direction were equally likely:  p = {pct(d['rayleigh_p'])}",
                         font_size=20, color=P.FG),
        ).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
        naive.move_to([3.2, 1.6, 0])
        self.play(FadeIn(dial), FadeIn(cap))
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in seen], lag_ratio=0.1),
                  run_time=1.2)
        self.play(Create(mean), FadeIn(naive))
        timing.hold_to_read(self, cap, naive, settle=0.8)

        # 2. but the surveys did not look everywhere
        w_max = max(sel["weight"])
        lobe = dial.polar(sel["lon_deg"], sel["weight"], w_max, 1.9, P.PURPLE)
        cap2 = layout.caption("Surveys find objects near perihelion, and they favoured one side",
                              font_size=22)
        note = VGroup(
            layout.label("survey sensitivity to ϖ", font_size=20, color=P.PURPLE, weight="BOLD"),
            layout.label(f"{sel['contrast']:.0f} times higher toward "
                         f"ϖ = {sel['phi1_deg']:.0f}° than away from it", font_size=18,
                         color=P.FG),
            layout.label("stand-in for the DES, OSSOS and Sheppard-Trujillo\n"
                         "pointing histories the paper simulates", font_size=15, color=P.FG,
                         line_spacing=0.9),
        ).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
        note.next_to(naive, DOWN, buff=0.6, aligned_edge=LEFT)
        self.play(FadeIn(lobe), FadeOut(cap), FadeIn(cap2), run_time=1.2)
        self.play(FadeIn(note))
        timing.hold_to_read(self, cap2, note, settle=1.0)
        self.play(FadeOut(VGroup(dial, seen, mean, lobe, naive, note, cap2)))

        # 3. how aligned is a uniform population, seen through the surveys?
        top = 3.5
        ax, labels = widgets.labeled_axes(
            [0, 1, 0.2], [0, top, 1], x_label=f"alignment strength R of {n} longitudes",
            y_label="probability density", y_rotate=True, numbers=True,
            x_length=9.6, y_length=4.1, shift_down=-0.35)
        flat = widgets.histogram(ax, nul["edges"], nul["flat"], color=P.MUTED, opacity=0.55)
        biased = widgets.histogram(ax, nul["edges"], nul["selection"], color=P.PURPLE,
                                   opacity=0.5)
        obs = widgets.marker_line(ax, d["r_bar"], (0, top), f"observed  R = {d['r_bar']:.2f}",
                                  color=P.GREEN, font_size=16, side=UP)
        flat_lab = layout.label(f"uniform population, uniform survey\n"
                                f"{pct(d['rayleigh_p'])} are this aligned", font_size=16,
                                color=P.FG, line_spacing=0.9)
        flat_lab.move_to(ax.c2p(0.7, 3.1), aligned_edge=LEFT)
        sel_lab = layout.label(f"uniform population, seen by the surveys\n"
                               f"{pct(d['consistency_p'])} are this aligned", font_size=16,
                               color=P.PURPLE, weight="BOLD", line_spacing=0.9)
        sel_lab.move_to(ax.c2p(0.7, 2.45), aligned_edge=LEFT)
        cap3 = layout.caption("Draw uniform populations and measure how aligned they look",
                              font_size=22)
        self.play(Create(ax), FadeIn(labels), FadeIn(cap3))
        self.play(FadeIn(flat, lag_ratio=0.05), FadeIn(flat_lab), run_time=1.2)
        self.play(Create(obs))
        timing.hold_to_read(self, cap3, flat_lab, settle=0.6)
        band = d["paper_band"]
        cap4 = layout.caption(
            f"Through the surveys the alignment is ordinary  (paper: {pct(band[0])} to "
            f"{pct(band[1])})", font_size=22)
        self.play(FadeIn(biased, lag_ratio=0.05), FadeIn(sel_lab), FadeOut(cap3), FadeIn(cap4),
                  run_time=1.4)
        timing.hold_to_read(self, cap4, sel_lab, settle=1.2)
        self.play(FadeOut(VGroup(ax, labels, flat, biased, obs, flat_lab, sel_lab, cap4)))

        # 4. the verdict depends on where the surveys looked
        dial2 = Dial(radius=1.8, centre=(-3.3, 0.5, 0.0), title="longitude of perihelion ϖ")
        seen2 = dial2.dots(varpis)
        lobe2 = dial2.polar(sel["lon_deg"], sel["weight"], w_max, 1.7, P.PURPLE)
        on = layout.label(
            f"survey bias points at the cluster:  p = {pct(d['p_directional_aligned'])}",
            font_size=20, color=P.PURPLE)
        off = layout.label(
            f"same bias, turned 90° away:  p = {pct(d['p_directional_rotated'])}",
            font_size=20, color=P.ORANGE)
        rows = VGroup(on, off).arrange(DOWN, buff=0.4, aligned_edge=LEFT).move_to([3.2, 0.7, 0])
        off.set_opacity(0)
        cap5 = layout.caption("So the answer rests on knowing exactly where each survey looked",
                              font_size=22)
        self.play(FadeIn(dial2), FadeIn(seen2), FadeIn(lobe2), FadeIn(on), FadeIn(cap5))
        self.wait(0.8)
        self.play(Rotate(lobe2, np.pi / 2, about_point=dial2.centre), run_time=1.6)
        self.play(lobe2.animate.set_color(P.ORANGE), off.animate.set_opacity(1))
        timing.hold_to_read(self, cap5, rows, settle=1.2)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, f"Seen through the surveys, random orbits align this well "
                  f"{pct(d['consistency_p'])} of the time.")
