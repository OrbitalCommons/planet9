"""Bernardinelli et al. (2021) -- the six-year DES search for outer solar system
objects.

A catalogue whose value is its known selection: one footprint, one depth, one
detection efficiency curve. Reproduced in p9-2021-des-catalog, which encodes
eight of the sixteen extreme objects; the footprint, efficiency curve, object
magnitudes and isotropy p-values shown are the crate's own
(anim.json -> papers -> p9-2021-des-catalog).
"""
from manim import (
    DOWN,
    LEFT,
    UP,
    Create,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Rectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2021-des-catalog"

ANGLE_NAMES = {"Omega": "node Ω", "omega": "argument ω", "varpi": "longitude ϖ"}


def footprint(m, bands):
    g = VGroup()
    for b in bands:
        g.add(m.box(b["ra_start_deg"], b["ra_end_deg"], b["dec_min_deg"], b["dec_max_deg"],
                    color=P.PURPLE, opacity=0.28, stroke_width=0))
    return g


def p_text(p):
    return f"{p:.2f}" if p >= 0.095 else f"{p:.3f}"


def isotropy_table(flat, sel, alpha, centre):
    """Three angles; p against a uniform sky and through the DES selection."""
    cw, ch = 3.0, 0.95
    g = VGroup()
    x_lab, x_a, x_b = centre[0] - 3.1, centre[0], centre[0] + cw
    y0 = centre[1] + ch
    for x, text, col in ((x_a, "uniform sky", P.FG), (x_b, "DES selection", P.PURPLE)):
        g.add(layout.label(text, font_size=21, color=col, weight="BOLD")
              .move_to([x, y0 + ch * 0.85, 0]))
    for i, (f, s) in enumerate(zip(flat, sel)):
        y = y0 - i * ch
        lab = layout.label(ANGLE_NAMES[f["angle"]], font_size=21, color=P.FG)
        g.add(lab.move_to([x_lab, y, 0]))
        for x, cell in ((x_a, f), (x_b, s)):
            hit = cell["mc_p"] < alpha
            box = Rectangle(width=cw - 0.12, height=ch - 0.1, stroke_width=1.2,
                            color=P.ORANGE if hit else P.MUTED)
            box.set_fill(P.ORANGE if hit else P.MUTED, opacity=0.35 if hit else 0.08)
            box.move_to([x, y, 0])
            txt = layout.label(f"p = {p_text(cell['mc_p'])}", font_size=22, color=P.FG,
                               weight="BOLD" if hit else "NORMAL")
            g.add(VGroup(box, txt.move_to(box)))
    return g


class DesCatalog2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        eff = d["efficiency"]

        self.add(paper.scene_header(CRATE))

        # 1. the survey and its haul
        m = sky.SkyMap(width=9.3, dec_range=(-80, 40), centre=(-1.6, 0.6, 0.0), ra_centre=0.0)
        ecl, gal = m.reference_curves()
        foot = footprint(m, d["footprint"])
        key = m.legend([("ecliptic", P.ORANGE), ("galactic plane ±10°", P.PURPLE)])
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), FadeIn(foot), FadeIn(key), run_time=1.2)
        haul = VGroup(
            layout.label(f"{d['footprint_deg2']:,.0f} deg²", font_size=24, color=P.PURPLE,
                         weight="BOLD"),
            layout.label(f"{100 * d['sky_fraction']:.0f}% of the sky, six years", font_size=16,
                         color=P.FG),
            layout.label(f"{d['n_tnos']} objects", font_size=24, color=P.GREEN, weight="BOLD"),
            layout.label(f"{d['n_new_tnos']} of them new", font_size=16, color=P.FG),
            layout.label(f"{d['n_extreme']} extreme", font_size=24, color=P.GREEN, weight="BOLD"),
            layout.label("a > 150 AU, q > 30 AU", font_size=16, color=P.FG),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        for k in (2, 4):
            VGroup(*haul[k:]).shift(DOWN * 0.22)
        haul.move_to([5.2, 0.6, 0])
        cap = layout.caption("Six years of the Dark Energy Survey, searched as one data set",
                             font_size=22)
        self.play(FadeIn(haul[0]), FadeIn(haul[1]), FadeIn(cap))
        self.play(FadeIn(VGroup(*haul[2:4])))
        timing.hold_to_read(self, cap, VGroup(*haul[:4]), settle=0.6)

        peri = m.dots(objs, color=P.GREEN, radius=0.07, opacity=1.0)
        cap2 = layout.caption(
            f"Where {len(objs)} of its {d['n_extreme']} extreme orbits come closest to the Sun",
            font_size=22)
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in peri], lag_ratio=0.15),
                  FadeIn(VGroup(*haul[4:])), FadeOut(cap), FadeIn(cap2), run_time=1.6)
        timing.hold_to_read(self, cap2, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, foot, key, peri, haul, cap2)))

        # 2. a known detection efficiency
        ax, labels = widgets.labeled_axes(
            [19, 27, 1], [0, 1, 0.25], x_label="apparent magnitude r  (fainter →)",
            y_label="fraction detected", y_rotate=True, numbers=True,
            x_length=10.0, y_length=4.0, shift_down=-0.4)
        curve = widgets.curve(ax, eff["r_mag"], eff["fraction"], color=P.PURPLE,
                              stroke_width=3.5)
        depth = widgets.marker_line(ax, d["depth_r"], (0, 1.0),
                                    f"half are found at r = {d['depth_r']:.1f}", font_size=16,
                                    side=UP)
        rug = VGroup(*[
            Line(ax.c2p(o["r_at_perihelion"], 0.0), ax.c2p(o["r_at_perihelion"], 0.09),
                 color=P.GREEN, stroke_width=3.5) for o in objs])
        rug_lab = layout.label("each extreme object's brightness at perihelion", font_size=16,
                               color=P.GREEN)
        rug_lab.move_to(ax.c2p(21.4, 0.2))
        cap3 = layout.caption("What makes it a test: the chance of finding each object is known",
                              font_size=22)
        self.play(Create(ax), FadeIn(labels), FadeIn(cap3))
        self.play(Create(curve), run_time=1.4)
        self.play(Create(depth))
        self.play(FadeIn(rug, lag_ratio=0.1), FadeIn(rug_lab))
        timing.hold_to_read(self, cap3, rug_lab, settle=1.2)
        self.play(FadeOut(VGroup(ax, labels, curve, depth, rug, rug_lab, cap3)))

        # 3. do the extreme objects point anywhere in particular?
        table = isotropy_table(d["battery_flat"], d["battery_selection"], 0.05, (0.0, 0.7, 0))
        cap4 = layout.caption(
            f"The angles of the {len(objs)} encoded extreme objects, tested for isotropy",
            font_size=22)
        verdict = layout.label(
            f"through the survey's own selection the smallest p is "
            f"{p_text(d['min_p_selection'])}", font_size=21, color=P.PURPLE, weight="BOLD")
        verdict.move_to([0.0, -2.0, 0])
        self.play(FadeIn(table, lag_ratio=0.05), FadeIn(cap4), run_time=1.4)
        self.play(FadeIn(verdict))
        timing.hold_to_read(self, cap4, verdict, settle=1.4)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, f"{d['n_tnos']} objects with a known selection: "
                  "the extreme ones look isotropic.")
