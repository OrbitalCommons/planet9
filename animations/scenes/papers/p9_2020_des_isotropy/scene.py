"""Bernardinelli et al. (2020) -- is the DES sample of extreme TNOs isotropic?

DES watched one patch of southern sky, so the orbits it can find are not a fair
sample of directions. The paper compares the detected angles with an isotropic
population seen through the DES selection function: twelve Kuiper tests (three
angles, four sample definitions). Reproduced in p9-2020-des-isotropy; the
footprint, perihelion directions, null distribution and every p-value shown are
the crate's own (anim.json -> papers -> p9-2020-des-isotropy).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    AnnularSector,
    Circle,
    Create,
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
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2020-des-isotropy"

ANGLE_NAMES = {"Omega": "node Ω", "omega": "argument ω", "varpi": "longitude ϖ"}


class Dial(VGroup):
    """A compass of orbital longitude: 0° to the right, increasing anticlockwise."""

    def __init__(self, radius=1.9, centre=(0.0, 0.0, 0.0), title=None):
        super().__init__()
        self.radius = radius
        self.centre = np.array(centre, dtype=float)
        ring = Circle(radius=radius, color=P.MUTED, stroke_width=1.5).move_to(self.centre)
        self.add(ring)
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

    def dots(self, degs, color=P.GREEN, radius=0.075, r=None):
        return VGroup(*[Dot(self.p(d, r), radius=radius, color=color).set_z_index(3)
                        for d in degs])

    def sectors(self, centres, width, values, color, r_inner, r_outer, opacity=0.55):
        """Radial histogram: one wedge per bin, length proportional to value."""
        vmax = max(values)
        g = VGroup()
        for c, v in zip(centres, values):
            if v <= 0.002 * vmax:
                continue
            g.add(AnnularSector(
                inner_radius=r_inner, outer_radius=r_inner + (r_outer - r_inner) * v / vmax,
                angle=np.deg2rad(width), start_angle=np.deg2rad(c - width / 2),
                color=color, fill_opacity=opacity, stroke_width=0.6,
                arc_center=self.centre))
        return g


def footprint(m, bands):
    g = VGroup()
    for b in bands:
        g.add(m.box(b["ra_start_deg"], b["ra_end_deg"], b["dec_min_deg"], b["dec_max_deg"],
                    color=P.PURPLE, opacity=0.28, stroke_width=0))
    return g


def p_text(p):
    return f"{p:.2f}" if p >= 0.095 else f"{p:.3f}"


def battery_grid(cells, cases, alpha, title, centre):
    """Three angles by four sample definitions; a cell lights up when its
    Kuiper test rejects isotropy at `alpha`."""
    cw, ch = 1.2, 0.8
    angles = ["Omega", "omega", "varpi"]
    g = VGroup()
    x0 = centre[0] - cw * (len(cases) - 1) / 2 + 0.7
    y0 = centre[1] + ch
    for j, case in enumerate(cases):
        head = layout.label(f"a > {case['a_min']:.0f}\nq > {case['q_min']:.0f}", font_size=14,
                            color=P.MUTED, line_spacing=0.8)
        g.add(head.move_to([x0 + j * cw, y0 + ch * 0.95, 0]))
    for i, ang in enumerate(angles):
        lab = layout.label(ANGLE_NAMES[ang], font_size=15, color=P.FG)
        lab.move_to([x0 - cw * 0.5 - 0.12, y0 - i * ch, 0], aligned_edge=RIGHT)
        g.add(lab)
        for j, case in enumerate(cases):
            cell = next(c for c in cells if c["case"] == case["case"] and c["angle"] == ang)
            hit = cell["kuiper_p"] < alpha
            box = Rectangle(width=cw - 0.08, height=ch - 0.08, stroke_width=1.2,
                            color=P.ORANGE if hit else P.MUTED)
            box.set_fill(P.ORANGE if hit else P.MUTED, opacity=0.35 if hit else 0.08)
            box.move_to([x0 + j * cw, y0 - i * ch, 0])
            txt = layout.label(p_text(cell["kuiper_p"]), font_size=18,
                               color=P.FG, weight="BOLD" if hit else "NORMAL")
            g.add(VGroup(box, txt.move_to(box)))
    head = layout.label(title, font_size=20, color=P.FG, weight="BOLD")
    head.move_to([x0 + cw * (len(cases) - 1) / 2, y0 + ch * 1.9, 0])
    g.add(head)
    return g


class DesIsotropy2020(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = [o for o in d["objects"] if o["y4_discovery"]]
        null = d["null_varpi"]

        self.add(paper.scene_header(CRATE))

        # 1. one patch of sky
        m = sky.SkyMap(width=12.4, dec_range=(-80, 40), centre=(0.0, 0.65, 0.0), ra_centre=0.0)
        ecl, gal = m.reference_curves()
        foot = footprint(m, d["footprint"])
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), FadeIn(foot), run_time=1.2)
        key = m.legend([("ecliptic", P.ORANGE), ("galactic plane ±10°", P.PURPLE)])
        foot_lab = layout.label(f"DES footprint  {d['footprint_deg2']:,.0f} deg²", font_size=15,
                                color=P.PURPLE, weight="BOLD")
        foot_lab.move_to(m.p(30, -72))
        self.play(FadeIn(key), FadeIn(foot_lab))
        cap = layout.caption("The Dark Energy Survey watched one patch of the southern sky",
                             font_size=22)
        self.play(FadeIn(cap))
        timing.hold_to_read(self, cap, settle=0.6)

        peri = m.dots(objs, color=P.GREEN, radius=0.07, opacity=1.0)
        cap2 = layout.caption(
            f"Where each of its {len(objs)} extreme orbits (a > 150 AU) comes closest to the Sun",
            font_size=22)
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in peri], lag_ratio=0.15),
                  FadeOut(cap), FadeIn(cap2), run_time=1.6)
        timing.hold_to_read(self, cap2, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, foot, key, foot_lab, peri, cap2)))

        # 2. what an isotropic population looks like through that window
        dial = Dial(radius=2.0, centre=(-3.3, 0.55, 0.0), title="longitude of perihelion ϖ")
        seen = dial.dots([o["varpi_deg"] for o in objs], color=P.GREEN)
        cap3 = layout.caption(f"Their {len(objs)} longitudes of perihelion look bunched",
                              font_size=22)
        flat_p = layout.label(
            f"against a uniform sky:  p = {p_text(d['p_varpi_flat'])}", font_size=22, color=P.FG)
        flat_p.move_to([3.0, 1.6, 0])
        self.play(FadeIn(dial), FadeIn(cap3))
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in seen], lag_ratio=0.12),
                  run_time=1.2)
        self.play(FadeIn(flat_p))
        timing.hold_to_read(self, cap3, flat_p, settle=0.8)

        wedges = dial.sectors(null["centre_deg"], null["bin_deg"], null["selection"],
                              P.PURPLE, 0.25, 1.85)
        cap4 = layout.caption(
            "Random orbits seen through the DES window already bunch near ϖ ≈ 0°",
            font_size=22)
        sel_p = layout.label(
            f"through the DES footprint:  p = {p_text(d['p_varpi_selection'])}", font_size=22,
            color=P.PURPLE, weight="BOLD")
        sel_p.next_to(flat_p, DOWN, buff=0.4, aligned_edge=LEFT)
        note = layout.label(
            f"Kuiper test on ϖ, n = {len(objs)}\n"
            "purple: isotropic orbits whose perihelion\nfalls inside the footprint\n"
            f"p ≈ 1: the real {len(objs)} are no more\nbunched than the survey forces",
            font_size=17, color=P.FG, line_spacing=0.9)
        note.next_to(sel_p, DOWN, buff=0.5, aligned_edge=LEFT)
        self.play(FadeIn(wedges, lag_ratio=0.1), FadeOut(cap3), FadeIn(cap4), run_time=1.4)
        self.play(FadeIn(sel_p), FadeIn(note))
        timing.hold_to_read(self, cap4, sel_p, note, settle=1.2)
        self.play(FadeOut(VGroup(dial, seen, wedges, flat_p, sel_p, note, cap4)))

        # 3. the twelve tests
        left = battery_grid(d["battery_flat"], d["cases"], d["alpha"],
                            "against a uniform sky", (-3.6, 0.6, 0))
        right = battery_grid(d["battery_selection"], d["cases"], d["alpha"],
                             "through the DES selection", (3.5, 0.6, 0))
        n_flat, n_sel, n = d["n_significant_flat"], d["n_significant_selection"], d["n_tests"]
        tally_l = layout.label(f"{n_flat} of {n} reject isotropy (p < {d['alpha']})",
                               font_size=19, color=P.ORANGE)
        tally_l.move_to([-2.9, -1.75, 0])
        tally_r = layout.label(
            f"{n_sel} of {n} reject isotropy   (paper: {d['paper_n_significant']} of "
            f"{d['paper_n_tests']})", font_size=19, color=P.ORANGE, weight="BOLD")
        tally_r.move_to([3.9, -1.75, 0])
        cap5 = layout.caption("Kuiper p-values: three angles, four definitions of 'extreme'",
                              font_size=22)
        self.play(FadeIn(left, lag_ratio=0.05), FadeIn(cap5), run_time=1.2)
        self.play(FadeIn(tally_l))
        timing.hold_to_read(self, cap5, tally_l, settle=0.8)
        self.play(FadeIn(right, lag_ratio=0.05), run_time=1.2)
        self.play(FadeIn(tally_r))
        timing.hold_to_read(self, tally_r, settle=1.4)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, f"Through its own footprint DES looks isotropic: "
                  f"smallest p = {d['min_p_selection']:.2f}.")
