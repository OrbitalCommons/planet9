"""Pichierri & Batygin (2025) -- measuring the clustering and diffusion of
trans-Neptunian objects.

Clones of each distant object are integrated to measure how fast its semi-major
axis wanders under Neptune's kicks; the clustering of perihelia and of orbital
poles is then measured by stability class. Reproduced in p9-2025-clustering at
reduced scale on the workspace's 10-object sample; the diffusion coefficients,
classes, von Mises fit and poles shown are the crate's own
(anim.json -> papers -> p9-2025-clustering).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Circle,
    Create,
    DashedLine,
    DashedVMobject,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    MathTex,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-clustering"

CLASS_COLOUR = {"stable": P.GREEN, "metastable": P.ORANGE, "unstable": P.RED}


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

    def bump(self, degs, values, vmax, r_base, r_max, color, opacity=0.2, stroke_width=2.4):
        """A density drawn outward from a base circle of radius `r_base`."""
        m = VMobject(color=color, stroke_width=stroke_width)
        m.set_points_as_corners([self.p(d, r_base + (r_max - r_base) * v / vmax)
                                 for d, v in zip(degs, values)])
        return m.set_fill(color, opacity=opacity)


def log_axes(x_range, decades, x_label, y_label, x_length, y_length, shift):
    """Axes whose vertical coordinate counts decades above 10^decades[0].
    Returns (frame, axes, y) where y(value) is the axis coordinate of a value."""
    lo, hi = decades
    ax = widgets.axes(x_range, [0, hi - lo, 1], x_length=x_length, y_length=y_length,
                      shift_down=0)
    ax.get_x_axis().add_numbers(font_size=16)
    extras = VGroup()
    for k in range(lo, hi + 1):
        t = MathTex(rf"10^{{{k}}}", color=P.MUTED).scale(0.5)
        extras.add(t.next_to(ax.c2p(x_range[0], k - lo), LEFT, buff=0.15))
    xl = layout.label(x_label, font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.12)
    yl = layout.label(y_label, font_size=15, color=P.FG).rotate(np.pi / 2)
    yl.next_to(extras, LEFT, buff=0.15)
    extras.add(xl, yl)

    def y(value):
        return float(np.clip(np.log10(max(value, 1e-300)), lo, hi)) - lo

    return VGroup(ax, extras).shift(shift), ax, y


class Clustering2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        cfg = d["config"]
        law = d["diffusion_law"]
        n = d["n_sample"]

        self.add(paper.scene_header(CRATE))

        # 1. how fast does each orbit wander?
        frame, ax, ylog = log_axes([30, 85, 10], (-7, -2), "perihelion distance q (AU)",
                                   "diffusion of semi-major axis (AU² / yr)", 7.8, 4.3,
                                   np.array([-1.5, 0.45, 0.0]))
        theory = widgets.curve(ax, law["q_au"], [ylog(v) for v in law["d"]], color=P.ORANGE,
                               stroke_width=3)
        theory_lab = layout.label("analytical law for\nNeptune's kicks", font_size=16,
                                  color=P.ORANGE, line_spacing=0.9)
        theory_lab.move_to(ax.c2p(77, ylog(1.6e-4)))
        crit = ylog(0.5 * d["d_crit"])
        gate = DashedLine(ax.c2p(30, crit), ax.c2p(85, crit), color=P.RED, stroke_width=2)
        gate_lab = layout.label("stable below this line", font_size=16, color=P.RED)
        gate_lab.next_to(ax.c2p(85, crit), UP, buff=0.1, aligned_edge=RIGHT)
        cap = layout.caption("Neptune kicks a distant orbit at every perihelion passage",
                             font_size=22)
        self.play(FadeIn(frame), FadeIn(cap))
        self.play(Create(theory), FadeIn(theory_lab), run_time=1.4)
        self.play(Create(gate), FadeIn(gate_lab))
        timing.hold_to_read(self, cap, theory_lab, gate_lab, settle=0.5)

        marks = VGroup()
        for o in objs:
            col = CLASS_COLOUR[o["class"]]
            marks.add(VGroup(
                Line(ax.c2p(o["q"], ylog(o["d_mean"] - o["d_std"])),
                     ax.c2p(o["q"], ylog(o["d_mean"] + o["d_std"])), color=col, stroke_width=2),
                Dot(ax.c2p(o["q"], ylog(o["d_mean"])), radius=0.075, color=col).set_z_index(3)))
        fastest = max(objs, key=lambda o: o["d_mean"])
        tag = layout.label(fastest["name"], font_size=15, color=P.GREEN)
        tag.next_to(ax.c2p(fastest["q"], ylog(fastest["d_mean"])), RIGHT, buff=0.2)
        tag.shift(DOWN * 0.25)
        counts = VGroup(
            layout.label(f"{d['n_stable']} stable", font_size=22, color=P.GREEN, weight="BOLD"),
            layout.label(f"{d['n_metastable']} metastable", font_size=22, color=P.ORANGE,
                         weight="BOLD"),
            layout.label(f"{d['n_unstable']} unstable", font_size=22, color=P.RED,
                         weight="BOLD"),
            layout.label(f"here: {cfg['n_clones']} clones, {cfg['t_total_yr'] / 1e6:.0f} Myr\n"
                         f"paper: {cfg['paper_n_clones']} clones, "
                         f"{cfg['paper_t_total_yr'] / 1e9:.0f} Gyr", font_size=15,
                         color=P.FG, line_spacing=0.9),
        ).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
        counts.move_to([5.2, 0.6, 0])
        cap2 = layout.caption(
            f"Measured on clones of {n} real objects: all wander slowly enough to be stable",
            font_size=22)
        self.play(LaggedStart(*[FadeIn(x) for x in marks], lag_ratio=0.12), FadeIn(tag),
                  FadeOut(cap), FadeIn(cap2), run_time=1.8)
        self.play(FadeIn(counts, lag_ratio=0.15))
        timing.hold_to_read(self, cap2, counts, settle=1.0)
        self.play(FadeOut(VGroup(frame, theory, theory_lab, gate, gate_lab, marks, tag, counts,
                                 cap2)))

        # 2. where the long-lived orbits point
        fit = d["fit"]
        vmax = max(max(fit["density"]), max(fit["paper_density"]))
        dial = Dial(radius=2.0, centre=(-3.3, 0.5, 0.0), title="longitude of perihelion ϖ")
        dots = VGroup(*[Dot(dial.p(o["varpi_deg"]), radius=0.075,
                            color=CLASS_COLOUR[o["class"]]).set_z_index(3) for o in objs])
        base = Circle(radius=0.55, color=P.MUTED, stroke_width=1).move_to(dial.centre)
        base.set_stroke(opacity=0.6)
        mine = dial.bump(fit["lon_deg"], fit["density"], vmax, 0.55, 1.85, P.GREEN)
        theirs = DashedVMobject(
            dial.bump(fit["lon_deg"], fit["paper_density"], vmax, 0.55, 1.85, P.FG, opacity=0.0,
                      stroke_width=2.0), num_dashes=70)
        rows = VGroup(
            layout.label(f"fitted here:  centre ϖ = {d['mean_varpi_deg']:.0f}°,  "
                         f"spread {d['spread_rad']:.1f} rad", font_size=19, color=P.GREEN,
                         weight="BOLD"),
            layout.label(f"paper (dashed):  centre ϖ = {d['paper_mean_varpi_deg']:.0f}°,  "
                         f"spread {d['paper_spread_rad']:.1f} rad", font_size=19, color=P.FG),
            layout.label(f"chance of this alignment, allowing for\nwhere surveys look:  "
                         f"p = {100 * d['bias_mc_p']:.1f}%", font_size=17, color=P.FG,
                         line_spacing=0.9),
            layout.label("paper: unstable orbits instead split\ninto two groups, near "
                         + " and ".join(f"{m:.0f}°" for m in d["paper_unstable_modes_deg"]),
                         font_size=17, color=P.RED, line_spacing=0.9),
        ).arrange(DOWN, buff=0.35, aligned_edge=LEFT)
        modes = VGroup(*[
            Line(dial.p(m, dial.radius - 0.25), dial.p(m, dial.radius + 0.25), color=P.RED,
                 stroke_width=5) for m in d["paper_unstable_modes_deg"]])
        rows.move_to([3.3, 0.7, 0])
        cap3 = layout.caption(
            f"The {d['n_kept']} stable orbits point one way: a von Mises bump fitted to them",
            font_size=22)
        self.play(FadeIn(dial), FadeIn(cap3))
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in dots], lag_ratio=0.1),
                  run_time=1.2)
        self.play(FadeIn(base), FadeIn(mine), FadeIn(rows[0]))
        self.play(Create(theirs), FadeIn(rows[1]))
        self.play(FadeIn(rows[2]))
        self.play(FadeIn(rows[3]), Create(modes))
        timing.hold_to_read(self, cap3, rows, settle=1.2)
        self.play(FadeOut(VGroup(dial, dots, base, mine, theirs, rows, modes, cap3)))

        # 3. and how their planes tilt
        mp, lp = d["mean_pole"], d["laplace_pole"]
        ax2, labels2 = widgets.labeled_axes(
            [-30, 30, 10], [-30, 30, 10], x_label="i cos Ω  (degrees)",
            y_label="i sin Ω  (degrees)", y_rotate=True, numbers=False,
            x_length=4.5, y_length=4.5, shift_down=-0.65)
        VGroup(ax2, labels2).shift(LEFT * 3.2)
        rings = VGroup(*[
            Circle(radius=np.linalg.norm(ax2.c2p(r, 0) - ax2.c2p(0, 0)), color=P.MUTED,
                   stroke_width=1).set_stroke(opacity=0.5).move_to(ax2.c2p(0, 0))
            for r in (10, 20, 30)])
        ring_labs = VGroup(*[
            layout.label(f"{r}°", font_size=12, color=P.MUTED)
            .move_to(ax2.c2p(r * 0.72, -r * 0.72) + np.array([0.18, -0.05, 0]))
            for r in (10, 20, 30)])
        poles = VGroup(*[Dot(ax2.c2p(o["pole_x_deg"], o["pole_y_deg"]), radius=0.07,
                             color=CLASS_COLOUR[o["class"]]) for o in objs])
        planets = Dot(ax2.c2p(lp["x_deg"], lp["y_deg"]), radius=0.09, color=P.ORANGE)
        mean = Dot(ax2.c2p(mp["x_deg"], mp["y_deg"]), radius=0.12, color=P.FG).set_z_index(4)
        link = Line(planets.get_center(), mean.get_center(), color=P.FG, stroke_width=2)
        key = VGroup(
            layout.label(f"orbital poles of the {n} objects", font_size=18, color=P.GREEN),
            layout.label("pole of the giant planets' plane", font_size=18, color=P.ORANGE),
            layout.label(f"their mean pole, tilted {mp['offset_deg']:.0f}° from it",
                         font_size=18, color=P.FG, weight="BOLD"),
        ).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
        key.move_to([3.3, 0.7, 0])
        cap4 = layout.caption("Each orbit's pole, seen from above the solar system",
                              font_size=22)
        self.play(FadeIn(VGroup(ax2, labels2, rings, ring_labs)), FadeIn(cap4))
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in poles], lag_ratio=0.1),
                  FadeIn(key[0]), run_time=1.2)
        self.play(FadeIn(planets), FadeIn(key[1]))
        self.play(Create(link), FadeIn(mean), FadeIn(key[2]))
        timing.hold_to_read(self, cap4, key, settle=1.2)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, f"Long-lived orbits point toward ϖ = {d['mean_varpi_deg']:.0f}°, "
                  f"their planes tilted {mp['offset_deg']:.0f}° together.")
