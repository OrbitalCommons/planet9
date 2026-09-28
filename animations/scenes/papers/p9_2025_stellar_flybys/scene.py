"""Hu, Huang, Gladman & Zhu (2025) -- early stellar flybys are unlikely.

A passing star can lift perihelia and make sednoids, but it must also leave
them at low inclination (i < 30 deg) and apsidally aligned. Only two families
of encounter orientation do that: a star passing in the plane of the disc, or
one crossing it perpendicularly with its periastron on the ecliptic. The scene
throws the crate's random orientations onto an equal-area map and keeps the
ones inside the crate's acceptance bands, then shows how often the birth
cluster delivers any pass inside 1000 AU, and multiplies the factors out.
Data: anim.json -> papers -> p9-2025-stellar-flybys.
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
    Polygon,
    Rectangle,
    Scene,
    TransformFromCopy,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-stellar-flybys"


def pct(x):
    return f"{100 * x:.1f}%" if x < 0.1 else f"{100 * x:.0f}%"


class StellarFlybys2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        cl = d["cluster"]

        self.add(paper.scene_header(CRATE))

        # 1. every direction a star could come from, on an equal-area map
        W, H = 7.4, 4.4
        x0, y0 = -4.55, -2.0

        def p(w_deg, cos_i):
            return np.array([x0 + W * w_deg / 360.0, y0 + H * (cos_i + 1) / 2, 0.0])

        box = Rectangle(width=W, height=H, color=P.MUTED, stroke_width=1.3)
        box.move_to(p(180, 0))
        xt = VGroup(*[layout.label(f"{w}°", font_size=14, color=P.MUTED)
                      .next_to(p(w, -1), DOWN, buff=0.1) for w in (0, 90, 180, 270, 360)])
        xl = layout.label("where the star's closest approach lies  (argument of periastron)",
                          font_size=15).next_to(p(180, -1), DOWN, buff=0.4)
        yt = VGroup(
            layout.label("in the plane", font_size=14, color=P.MUTED).next_to(p(0, 1), LEFT, buff=0.1),
            layout.label("perpendicular", font_size=14, color=P.MUTED).next_to(p(0, 0), LEFT, buff=0.1),
            layout.label("in the plane,", font_size=14, color=P.MUTED).next_to(p(0, -0.84), LEFT, buff=0.1),
            layout.label("reversed", font_size=14, color=P.MUTED).next_to(p(0, -0.93), LEFT, buff=0.1),
        )
        yt[3].align_to(yt[2], RIGHT).next_to(yt[2], DOWN, buff=0.06, aligned_edge=RIGHT)
        yl = layout.label("tilt of the star's path", font_size=15).rotate(np.pi / 2)
        yl.next_to(yt, LEFT, buff=0.12)

        # the crate's acceptance bands (tolerance-wide; widened to be visible)
        tol = d["coplanar_tol_deg"]
        c_top = np.cos(np.radians(tol))
        s_half = np.sin(np.radians(d["symmetric_tol_deg"]))
        w_half = d["symmetric_tol_deg"]
        bands = VGroup()
        for lo, hi in ((c_top, 1.0), (-1.0, -c_top)):
            a, b = p(0, lo), p(360, hi)
            bands.add(Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0)
                      .set_fill(P.TEAL, opacity=0.35))
        for wc in (0.0, 180.0, 360.0):
            lo_w, hi_w = max(0.0, wc - w_half), min(360.0, wc + w_half)
            a, b = p(lo_w, -s_half), p(hi_w, s_half)
            bands.add(Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0)
                      .set_fill(P.TEAL, opacity=0.35))

        dots_ok, dots_bad = VGroup(), VGroup()
        for o in d["orientations"]:
            dot = Dot(p(o["arg_peri_deg"], o["cos_i"]), radius=0.028)
            if o["band"] == "rejected":
                dots_bad.add(dot.set_color(P.MUTED))
            else:
                dots_ok.add(dot.set_color(P.GREEN).scale(1.4))

        side = VGroup(
            layout.label("a sednoid-making star must", font_size=17, weight="BOLD"),
            layout.label(f"keep their tilts below {d['inclination_limit_deg']:.0f}°", font_size=16),
            layout.label("and keep their perihelia aligned", font_size=16),
            layout.label("only two geometries do:", font_size=16, color=P.TEAL),
            layout.label("• pass in the disc's plane", font_size=16, color=P.TEAL),
            layout.label("• cross it square-on, with", font_size=16, color=P.TEAL),
            layout.label("  closest approach on the ecliptic", font_size=16, color=P.TEAL),
        ).arrange(DOWN, buff=0.15, aligned_edge=LEFT).move_to([5.0, 0.9, 0])

        cap = layout.caption(f"{d['n_orientations']} random directions a star could arrive from",
                             font_size=22)
        self.play(FadeIn(box), FadeIn(xt), FadeIn(xl), FadeIn(yt), FadeIn(yl), FadeIn(cap),
                  run_time=0.9)
        self.play(LaggedStart(FadeIn(dots_bad, lag_ratio=0.002), FadeIn(dots_ok, lag_ratio=0.02),
                              lag_ratio=0.3), run_time=2.0)
        timing.hold_to_read(self, cap, settle=0.4)
        cap2 = layout.caption("Keep only the ones that leave the sednoids flat and aligned",
                              font_size=22)
        self.play(FadeIn(side), FadeIn(bands), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        self.play(dots_bad.animate.set_opacity(0.18), run_time=1.0)
        tally = paper.result_readout(
            "orientations that work",
            f"{d['n_accepted']} of {d['n_orientations']}  ({pct(d['sampled_fraction'])})",
            color=P.TEAL).scale(0.6)
        tally.next_to(side, DOWN, buff=0.45)
        self.play(FadeIn(tally), run_time=0.5)
        timing.hold_to_read(self, cap2, side, settle=1.0)
        self.play(FadeOut(VGroup(box, xt, xl, yt, yl, bands, dots_ok, dots_bad, side, cap2,
                                 tally)), run_time=0.8)

        # 2. how often does the birth cluster deliver a pass inside 1000 AU?
        q = np.array(d["q_grid_au"])
        ax, labs = widgets.labeled_axes(
            [0, 3000, 500], [0, 1, 0.25], x_label="closest approach of the passing star  (AU)",
            y_label="chance of at least one pass", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.0, shift_down=0)
        ax.move_to([-1.9, 0.2, 0])
        labs[0].next_to(ax, DOWN, buff=0.5)
        labs[1].next_to(ax, LEFT, buff=0.55)
        cols = [P.MUTED, P.GREEN, P.FG]
        curves, tags = VGroup(), VGroup()
        for c, col in zip(d["p_vs_q"], cols):
            curves.add(widgets.curve(ax, q, c["p"], color=col,
                                     stroke_width=3.2 if col == P.GREEN else 2.0))
            tags.add(layout.label(f"{c['density_per_pc3']:.0f} stars/pc³", font_size=14, color=col)
                     .next_to(ax.c2p(q[-1], c["p"][-1]), RIGHT, buff=0.1))
        qmax = d["q_star_max_au"]
        mark = widgets.marker_line(ax, qmax, (0, 1), f"{qmax:.0f} AU", color=P.PURPLE, side=UP)
        hit = Dot(ax.c2p(qmax, d["p_encounter"]), radius=0.08, color=P.GREEN)
        hit_lab = layout.label(pct(d["p_encounter"]), font_size=17, color=P.GREEN)
        hit_lab.next_to(hit, LEFT, buff=0.12).shift(UP * 0.15)
        nep, sedna = dataio.body("Neptune"), dataio.body("Sedna")
        scale_note = VGroup(
            layout.label("for scale:", font_size=15, color=P.MUTED),
            layout.label(f"Neptune orbits at {nep['a_au']:.0f} AU", font_size=15, color=P.MUTED),
            layout.label(f"Sedna swings out to {sedna['Q_au']:.0f} AU", font_size=15,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT).move_to([5.3, 0.3, 0])
        cluster = VGroup(
            layout.label("birth cluster:", font_size=15, color=P.MUTED),
            layout.label(f"{cl['density_per_pc3']:.0f} stars/pc³, {cl['velocity_dispersion_kms']:.0f} km/s,",
                         font_size=15, color=P.GREEN),
            layout.label(f"for {cl['residence_myr']:.0f} Myr", font_size=15, color=P.GREEN),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT).next_to(scale_note, DOWN, buff=0.35,
                                                             aligned_edge=LEFT)
        cap3 = layout.caption("Computed: the chance the Sun's birth cluster sent any star that close",
                              font_size=22)
        self.play(FadeIn(ax), FadeIn(labs), FadeIn(cap3), FadeIn(scale_note), FadeIn(cluster),
                  run_time=0.9)
        self.play(*[Create(c) for c in curves], FadeIn(tags), run_time=1.6)
        self.play(Create(mark), FadeIn(hit), FadeIn(hit_lab), run_time=0.8)
        timing.hold_to_read(self, cap3, cluster, settle=1.0)
        self.play(FadeOut(VGroup(ax, labs, curves, tags, mark, scale_note, cluster, cap3)),
                  hit.animate.move_to([-6.0, 1.6, 0]),
                  hit_lab.animate.move_to([-6.0, 1.6, 0]).set_opacity(0), run_time=0.8)
        self.remove(hit, hit_lab)

        # 3. multiply it out, on a log scale
        steps = [
            ("any early star", 1.0, P.MUTED),
            (f"passes inside {qmax:.0f} AU", d["p_encounter"], P.GREEN),
            ("... from a direction that works", d["p_after_geometry"], P.TEAL),
            ("... and makes enough sednoids", d["p_total"], P.ORANGE),
        ]
        lx0, lx1 = -1.2, 5.6
        lo_log = -3.0

        def bx(v):
            return lx0 + (lx1 - lx0) * (np.log10(v) - lo_log) / (0 - lo_log)

        rows = VGroup()
        for k, (name, v, col) in enumerate(steps):
            y = 1.9 - 1.05 * k
            bar = Polygon([lx0, y - 0.26, 0], [bx(v), y - 0.26, 0], [bx(v), y + 0.26, 0],
                          [lx0, y + 0.26, 0], stroke_width=0).set_fill(col, opacity=0.75)
            lab = layout.label(name, font_size=18, color=P.FG).next_to([lx0, y, 0], LEFT, buff=0.25)
            val = layout.label(pct(v), font_size=18, color=col, weight="BOLD")
            val.next_to([bx(v), y, 0], RIGHT, buff=0.15)
            rows.add(VGroup(bar, lab, val))
        axis = Line([lx0, -1.75, 0], [lx1, -1.75, 0], color=P.MUTED, stroke_width=1.3)
        ticks = VGroup()
        for v in (0.001, 0.01, 0.1, 1.0):
            ticks.add(Line([bx(v), -1.75, 0], [bx(v), -1.83, 0], color=P.MUTED, stroke_width=1.3))
            ticks.add(layout.label(pct(v), font_size=14, color=P.MUTED)
                      .next_to([bx(v), -1.83, 0], DOWN, buff=0.08))
        ceiling = DashedLine([bx(d["published_max"]), -1.75, 0], [bx(d["published_max"]), 2.4, 0],
                             color=P.RED, stroke_width=2)
        ceil_lab = layout.label(f"paper: at most {pct(d['published_max'])}", font_size=15,
                                color=P.RED).next_to(ceiling.get_end(), UP, buff=0.06)
        notes = VGroup(
            layout.label(f"× {pct(d['f_geometry'])} of orientations", font_size=14, color=P.TEAL),
            layout.label(f"× {100 * d['f_success']:.0f}% assumed success", font_size=14,
                         color=P.ORANGE),
        )
        notes[0].next_to(rows[2][1], DOWN, buff=0.08, aligned_edge=RIGHT)
        notes[1].next_to(rows[3][1], DOWN, buff=0.08, aligned_edge=RIGHT)
        cap4 = layout.caption("Multiply the chances: each step is a cut", font_size=22)
        self.play(FadeIn(axis), FadeIn(ticks), FadeIn(cap4), run_time=0.6)
        for k, row in enumerate(rows):
            if k == 0:
                self.play(FadeIn(row), run_time=0.6)
            else:
                self.play(TransformFromCopy(rows[k - 1][0], row[0]), FadeIn(row[1:]),
                          *( [FadeIn(notes[k - 2])] if k >= 2 else []), run_time=0.9)
        self.play(Create(ceiling), FadeIn(ceil_lab), run_time=0.6)
        timing.hold_to_read(self, cap4, rows, settle=1.2)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, f"A flyby that fits everything: about {pct(d['p_total'])} here, at most 5% in the paper.")
