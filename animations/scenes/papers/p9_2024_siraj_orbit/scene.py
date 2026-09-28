"""Siraj, Chyba & Tremaine (2024) -- orbit of a possible Planet X.

An independent fit of the unseen planet from the apsidal confinement of distant
objects on long-term stable orbits. Reproduced in p9-2024-siraj-orbit as a
mass-distance posterior: the confinement strength fixes m / a^3 and a prior on
the semi-major axis breaks the degeneracy. The posterior map, the ridges, the
best fit from the workspace's 10-object sample and the tensions shown are the
crate's own (anim.json -> papers -> p9-2024-siraj-orbit).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    Create,
    DashedVMobject,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2024-siraj-orbit"


def density_cells(ax, xs, ys, values, color, max_opacity=0.8, floor=0.04):
    """A gridded map (row-major, y rows by x columns) as translucent cells."""
    vals = np.asarray(values, dtype=float).reshape(len(ys), len(xs))
    dx, dy = xs[1] - xs[0], ys[1] - ys[0]
    g = VGroup()
    for iy, y in enumerate(ys):
        for ix, x in enumerate(xs):
            f = vals[iy, ix]
            if f < floor:
                continue
            a, b = ax.c2p(x - dx / 2, y - dy / 2), ax.c2p(x + dx / 2, y + dy / 2)
            cell = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0)
            g.add(cell.set_fill(color, opacity=max_opacity * min(f, 1.0)))
    return g


def clipped_curve(ax, xs, ys, y_max, color, stroke_width=3.0):
    pts = [(x, y) for x, y in zip(xs, ys) if y <= y_max]
    return widgets.curve(ax, [p[0] for p in pts], [p[1] for p in pts], color=color,
                         stroke_width=stroke_width)


def error_cross(ax, x, sx, y, sy, color):
    return VGroup(
        Line(ax.c2p(x - sx, y), ax.c2p(x + sx, y), color=color, stroke_width=3),
        Line(ax.c2p(x, y - sy), ax.c2p(x, y + sy), color=color, stroke_width=3),
        Dot(ax.c2p(x, y), radius=0.08, color=color),
    ).set_z_index(4)


class SirajOrbit2024(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        bb = d["bb21"]
        grid, ridge = d["grid"], d["ridge"]

        self.add(paper.scene_header(CRATE))

        # 1. the mass-distance plane
        m_top = 14
        ax, labels = widgets.labeled_axes(
            [150, 650, 100], [0, m_top, 2], x_label="semi-major axis of the planet (AU)",
            y_label="mass (Earth masses)", y_rotate=True, numbers=True,
            x_length=7.0, y_length=4.3, shift_down=-0.35)
        VGroup(ax, labels).shift(LEFT * 2.3)
        line = clipped_curve(ax, ridge["a"], ridge["mass"], m_top, P.ORANGE)
        cap = layout.caption(
            "How tightly the orbits are confined fixes one thing: mass over distance cubed",
            font_size=22)
        ridge_lab = layout.label(
            f"confinement of the {d['n_sample']}-object sample\n"
            f"alignment strength R = {d['r_bar']:.2f}", font_size=17, color=P.ORANGE,
            line_spacing=0.9)
        ridge_lab.move_to([4.5, 2.2, 0])
        self.play(Create(ax), FadeIn(labels), FadeIn(cap))
        self.play(Create(line), FadeIn(ridge_lab), run_time=1.4)
        timing.hold_to_read(self, cap, ridge_lab, settle=0.6)

        cells = density_cells(ax, grid["a"], grid["mass"], grid["density"], P.TEAL)
        here = Dot(ax.c2p(d["map_a"], d["map_mass"]), radius=0.09, color=P.GREEN).set_z_index(5)
        here_lab = layout.label(
            f"best fit here\n{d['map_mass']:.1f} M⊕ at {d['map_a']:.0f} AU", font_size=17,
            color=P.GREEN, weight="BOLD", line_spacing=0.9)
        here_lab.next_to(ridge_lab, DOWN, buff=0.4, aligned_edge=LEFT)
        cap2 = layout.caption(
            f"A distance prior of {d['a_prior']:.0f} ± {d['a_prior_sigma']:.0f} AU, taken from "
            "the paper, picks the point on the ridge", font_size=22)
        self.play(FadeIn(cells, lag_ratio=0.002), FadeOut(cap), FadeIn(cap2), run_time=1.4)
        self.play(FadeIn(here), FadeIn(here_lab))
        timing.hold_to_read(self, cap2, here_lab, settle=0.8)

        paper_line = DashedVMobject(
            clipped_curve(ax, ridge["a"], ridge["mass_paper_sample"], m_top, P.BLUE, 2.2),
            num_dashes=40)
        pub = error_cross(ax, d["a"], d["a_sigma"], d["mass"], d["mass_sigma"], P.BLUE)
        old = error_cross(ax, bb["a"], bb["a_sigma"], bb["mass"], bb["mass_sigma"], P.FG)
        pub_lab = layout.label(
            f"paper, stable sample (dashed ridge)\n{d['mass']:.1f} M⊕ "
            f"at {d['a']:.0f} ± {d['a_sigma']:.0f} AU", font_size=17, color=P.BLUE,
            weight="BOLD", line_spacing=0.9)
        pub_lab.next_to(here_lab, DOWN, buff=0.4, aligned_edge=LEFT)
        old_lab = layout.label(
            f"Brown & Batygin (2021)\n{bb['mass']:.1f} M⊕ at {bb['a']:.0f} AU", font_size=17,
            color=P.FG, line_spacing=0.9)
        old_lab.next_to(pub_lab, DOWN, buff=0.4, aligned_edge=LEFT)
        cap3 = layout.caption(
            "The paper's larger sample is less tightly aligned, so its planet is lighter",
            font_size=22)
        self.play(Create(paper_line), FadeIn(pub), FadeIn(pub_lab), FadeOut(cap2), FadeIn(cap3),
                  run_time=1.4)
        self.play(FadeIn(old), FadeIn(old_lab))
        timing.hold_to_read(self, cap3, pub_lab, old_lab, settle=1.2)
        self.play(FadeOut(VGroup(ax, labels, line, ridge_lab, cells, here, here_lab, paper_line,
                                 pub, old, pub_lab, old_lab, cap3)))

        # 2. two different planets
        s = 0.0058
        c_top = np.array([-3.5, 0.45, 0.0])
        sun = orbits.sun(radius=0.07).move_to(c_top)
        nep = orbits.ellipse_orbit(30.0 * s, 0.0, color=P.MUTED, stroke_width=1.2).shift(c_top)
        o_old = orbits.ellipse_orbit(bb["a"] * s, bb["e"], color=P.FG, varpi=np.pi,
                                     stroke_width=2.2, opacity=0.8).shift(c_top)
        o_new = orbits.ellipse_orbit(d["a"] * s, d["e"], color=P.BLUE, varpi=np.pi,
                                     stroke_width=3.5).shift(c_top)
        top_title = layout.label("from above", font_size=17, color=P.FG)
        top_title.move_to([-3.0, 2.85, 0])

        s2 = 0.0125
        c_side = np.array([1.2, -0.9, 0.0])
        ecl = Line(c_side, c_side + RIGHT * 5.4, color=P.ORANGE, stroke_width=1.6)
        ecl_lab = layout.label("plane of the known planets", font_size=13, color=P.ORANGE)
        ecl_lab.next_to(ecl, DOWN, buff=0.12)

        def tilted(a, i_deg, color, width):
            t = np.deg2rad(i_deg)
            tip = c_side + a * s2 * np.array([np.cos(t), np.sin(t), 0.0])
            return VGroup(Line(c_side, tip, color=color, stroke_width=width),
                          Dot(tip, radius=0.08, color=color))

        side_old = tilted(bb["a"], bb["i"], P.FG, 2.2)
        side_new = tilted(d["a"], d["i"], P.BLUE, 3.5)
        side_sun = orbits.sun(radius=0.07).move_to(c_side)
        side_title = layout.label("from the side", font_size=17, color=P.FG)
        side_title.move_to([3.9, 2.85, 0])
        old_tag = layout.label(
            f"Brown & Batygin (2021)\n{bb['mass']:.1f} M⊕,  a = {bb['a']:.0f} AU,  "
            f"i = {bb['i']:.0f}°", font_size=17, color=P.FG, line_spacing=0.9)
        old_tag.move_to([3.9, 2.0, 0])
        new_tag = layout.label(
            f"Siraj, Chyba & Tremaine (2024)\n{d['mass']:.1f} M⊕,  a = {d['a']:.0f} AU,  "
            f"e = {d['e']:.2f},  i = {d['i']:.1f}°", font_size=17, color=P.BLUE, weight="BOLD",
            line_spacing=0.9)
        new_tag.move_to([3.9, 1.05, 0])
        cap4 = layout.caption("Two fits to the same kind of evidence, two different planets",
                              font_size=22)
        self.play(FadeIn(VGroup(sun, nep, top_title, side_sun, ecl, ecl_lab, side_title)),
                  FadeIn(cap4))
        self.play(Create(o_old), Create(side_old[0]), FadeIn(side_old[1]), FadeIn(old_tag),
                  run_time=1.4)
        self.play(Create(o_new), Create(side_new[0]), FadeIn(side_new[1]), FadeIn(new_tag),
                  run_time=1.4)
        timing.hold_to_read(self, cap4, old_tag, new_tag, settle=1.0)

        cap5 = layout.caption(
            f"Paper: only {100 * d['paper_bb21_overlap']:.2f}% of the 2021 predicted orbits "
            "come within 1σ of the new best fit", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=1.2)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, f"A lighter, closer, flatter planet: {d['mass']:.1f} Earth masses at "
                  f"{d['a']:.0f} AU, tilted {d['i']:.0f}°.")
