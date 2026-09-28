"""Iorio (2026) -- has Saturn left room for Planet Nine and its smaller cousins?

Cassini pinned Saturn's orbit to metres. A distant planet's tide would make
that orbit precess, faster the closer the planet is, so the uncertainties of
Saturn's precessions bound where Planet Nine, Planet X and Planet Y can be.
The paper uses all of Saturn's orbital rates at once; the reproduction crate
carries the perihelion precession alone.

Everything drawn comes from anim.json -> papers -> p9-2026-iorio-precession:
the bound in the mass-distance plane, each candidate's reach from perihelion to
aphelion, and Saturn's induced precession against Planet Nine's true anomaly;
the paper's bound and verdicts appear as labelled comparisons.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    DashedVMobject,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2026-iorio-precession"


class Plot(VGroup):
    """Axes anchored at their lower-left corner, with optional log scales and
    hand-placed tick labels. ``c2p`` takes data values."""

    def __init__(self, x, y, x_ticks, y_ticks, x_label, y_label, size=(10.0, 4.2),
                 centre=(0.0, 0.2), xlog=False, ylog=False, tick_size=15):
        super().__init__()
        self.xlog, self.ylog = xlog, ylog
        self.x, self.y = x, y
        self.w, self.h = size
        self.corner = np.array([centre[0] - self.w / 2, centre[1] - self.h / 2, 0.0])
        self.add(Line(self.c2p(x[0], y[0]), self.c2p(x[1], y[0]), color=P.MUTED, stroke_width=2),
                 Line(self.c2p(x[0], y[0]), self.c2p(x[0], y[1]), color=P.MUTED, stroke_width=2))
        xt, yt = VGroup(), VGroup()
        for v, text in x_ticks.items():
            p = self.c2p(v, y[0])
            self.add(Line(p, p + DOWN * 0.08, color=P.MUTED, stroke_width=2))
            xt.add(layout.label(text, font_size=tick_size, color=P.MUTED).next_to(p, DOWN, buff=0.14))
        for v, text in y_ticks.items():
            p = self.c2p(x[0], v)
            self.add(Line(p, p + LEFT * 0.08, color=P.MUTED, stroke_width=2))
            yt.add(layout.label(text, font_size=tick_size, color=P.MUTED).next_to(p, LEFT, buff=0.14))
        self.add(xt, yt)
        self.add(layout.label(x_label, font_size=18).next_to(xt, DOWN, buff=0.14)
                 .set_x(self.corner[0] + self.w / 2))
        self.add(layout.label(y_label, font_size=15).rotate(np.pi / 2).next_to(yt, LEFT, buff=0.14)
                 .set_y(self.corner[1] + self.h / 2))

    @staticmethod
    def _f(v, log):
        return np.log10(v) if log else v

    def c2p(self, x, y):
        fx = (self._f(x, self.xlog) - self._f(self.x[0], self.xlog)) / (
            self._f(self.x[1], self.xlog) - self._f(self.x[0], self.xlog))
        fy = (self._f(y, self.ylog) - self._f(self.y[0], self.ylog)) / (
            self._f(self.y[1], self.ylog) - self._f(self.y[0], self.ylog))
        return self.corner + np.array([fx * self.w, fy * self.h, 0.0])


class Iorio2026(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        cases, pub, edge = d["cases"], d["published"], d["boundary"]
        bound = d["bound_mas_cy"]

        self.add(paper.scene_header(CRATE))

        # 1. the bound in the mass-distance plane, and what each candidate spans
        x0, x1, y0, y1 = 30.0, 1000.0, 0.03, 15.0
        ax = Plot(
            (x0, x1), (y0, y1),
            {30: "30", 100: "100", 300: "300", 1000: "1000"},
            {0.03: "0.03", 0.1: "0.1", 0.3: "0.3", 1: "1", 3: "3", 10: "10"},
            "distance from the Sun  (AU)", "mass  (M⊕)",
            size=(8.0, 4.2), centre=(-2.1, 0.45), xlog=True, ylog=True)
        mass = np.array(edge["mass_earth"])
        dist = np.array(edge["distance_au"])
        keep = (mass >= y0) & (mass <= y1)
        wall = widgets.curve(ax, dist[keep], mass[keep], color=P.RED, stroke_width=3.0)
        shade = Polygon(ax.c2p(x0, mass[keep][0]),
                        *[ax.c2p(x, m) for x, m in zip(dist[keep], mass[keep])],
                        ax.c2p(x0, mass[keep][-1]), stroke_width=0).set_fill(P.RED, opacity=0.14)
        paper_wall = DashedVMobject(
            widgets.curve(ax, np.array(edge["published_bound_distance_au"])[keep], mass[keep],
                          color=P.FG, stroke_width=2.0), num_dashes=40)
        key = VGroup(
            VGroup(Line(LEFT * 0.2, RIGHT * 0.2, color=P.RED, stroke_width=3),
                   layout.label(
                       f"reproduced: Saturn's perihelion\nturns {bound:g} mas per century",
                       font_size=14, color=P.RED, line_spacing=0.9)).arrange(RIGHT, buff=0.15),
            VGroup(DashedLine(LEFT * 0.2, RIGHT * 0.2, color=P.FG, stroke_width=2),
                   layout.label(
                       f"the paper's bound:\n{pub['bound_mas_cy']:g} mas per century",
                       font_size=14, line_spacing=0.9)).arrange(RIGHT, buff=0.15),
            VGroup(Line(LEFT * 0.2, RIGHT * 0.2, color=P.TEAL, stroke_width=5),
                   layout.label("part of an orbit\nSaturn allows", font_size=14, color=P.TEAL,
                                line_spacing=0.9)).arrange(RIGHT, buff=0.15),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        key.move_to([4.7, -0.75, 0])

        cap = layout.caption(
            "Closer than this line, a planet would turn Saturn's orbit faster than Cassini allows",
            font_size=22)
        self.play(FadeIn(ax), run_time=0.9)
        self.play(FadeIn(shade), Create(wall), Create(paper_wall), FadeIn(key[0]),
                  FadeIn(key[1]), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, key[0], key[1], settle=0.8)

        spans, names = VGroup(), VGroup()
        survivors = {}
        for c in cases:
            q, big_q, r_crit, m = (c["perihelion_au"], c["aphelion_au"],
                                   c["critical_distance_au"], c["mass_earth"])
            cut = min(max(r_crit, q), big_q)
            if cut > q:
                spans.add(Line(ax.c2p(q, m), ax.c2p(cut, m), color=P.RED, stroke_width=6))
            if cut < big_q:
                survivors[c["name"]] = Line(ax.c2p(cut, m), ax.c2p(big_q, m), color=P.TEAL,
                                            stroke_width=6)
                spans.add(survivors[c["name"]])
            far = c["published_min_distance_au"]
            if far is not None and far < big_q:
                p = ax.c2p(far, m)
                spans.add(Line(p + DOWN * 0.11, p + UP * 0.11, color=P.FG, stroke_width=3))
            says = f"paper: {c['published_verdict']}"
            if big_q > 2.0 * q:
                lab = layout.label(f"{c['name']}  ({says})", font_size=15)
                lab.next_to(ax.c2p(big_q, m), RIGHT, buff=0.12)
            elif r_crit <= q:
                lab = layout.label(f"{c['name']}\n({says})", font_size=15, line_spacing=0.9)
                lab.next_to(ax.c2p(big_q, m), RIGHT, buff=0.12)
            else:
                lab = layout.label(f"{c['name']}\n({says})", font_size=15, line_spacing=0.9)
                lab.next_to(ax.c2p(q, m), LEFT, buff=0.12)
            names.add(lab)
        cap2 = layout.caption(
            "Each candidate from perihelion to aphelion: only the far part of an orbit survives",
            font_size=22)
        self.play(FadeIn(spans, lag_ratio=0.1), FadeIn(names), FadeIn(key[2]), FadeOut(cap),
                  FadeIn(cap2), run_time=1.6)
        timing.hold_to_read(self, cap2, names, settle=1.0)

        cap3 = layout.caption(
            f"Paper: Planet X is ruled out; Planet Y only as a Mercury-mass body beyond "
            f"{pub['planet_y_mercury_min_a_au']:.0f} AU", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), run_time=0.6)
        timing.hold_to_read(self, cap3, settle=0.6)

        # at the paper's tighter bound the reproduction closes Planet X too
        px = next(c for c in cases if c["name"].startswith("Planet X"))
        closes = px["critical_distance_paper_bound_au"] >= px["aphelion_au"]
        verdict = "closes here too" if closes else "still stays open here"
        cap3b = layout.caption(
            f"At the paper's {pub['bound_mas_cy']:g} bound, Planet X's last arc {verdict} "
            f"({px['critical_distance_paper_bound_au']:.1f} vs {px['aphelion_au']:.1f} AU)",
            font_size=22)
        anims = [FadeOut(cap3), FadeIn(cap3b)]
        if closes and px["name"] in survivors:
            anims.append(survivors[px["name"]].animate.set_color(P.RED))
        self.play(*anims, run_time=1.0)
        timing.hold_to_read(self, cap3b, settle=0.8)
        self.play(FadeOut(VGroup(ax, wall, shade, paper_wall, key, spans, names, cap3b)))

        # 2. Planet Nine: Saturn's precession against where the planet is
        ax2 = Plot(
            (0, 360), (0.05, 10),
            {v: f"{v}°" for v in range(0, 361, 60)},
            {0.1: "0.1", 0.3: "0.3", 1: "1", 3: "3", 10: "10"},
            "true anomaly of Planet Nine  (0° = perihelion)",
            "Saturn's induced precession  (mas per century)",
            size=(9.2, 4.1), centre=(-0.75, 0.55), ylog=True)
        nines = [c for c in cases if c["published_allowed_deg"] is not None]
        shades = (P.BLUE, P.ORANGE)
        curves = VGroup()
        for c, color in zip(nines, shades):
            curves.add(widgets.curve(ax2, c["f_deg"], c["rate_mas_cy"], color=color,
                                     stroke_width=3.2))
        limit = Line(ax2.c2p(0, bound), ax2.c2p(360, bound), color=P.RED, stroke_width=2.5)
        limit_p = DashedLine(ax2.c2p(0, pub["bound_mas_cy"]), ax2.c2p(360, pub["bound_mas_cy"]),
                             color=P.FG, stroke_width=2)
        limits = VGroup(
            limit, limit_p,
            layout.label(f"reproduced bound: {bound:g}", font_size=15, color=P.RED)
            .next_to(ax2.c2p(360, bound), RIGHT, buff=0.12).shift(UP * 0.06),
            layout.label(f"the paper's bound: {pub['bound_mas_cy']:g}", font_size=15)
            .next_to(ax2.c2p(360, pub["bound_mas_cy"]), RIGHT, buff=0.12).shift(DOWN * 0.06),
        )
        heights = ((6.5, 4.3), (2.7, 1.8))
        bars = VGroup()
        for c, color, (h_mine, h_paper) in zip(nines, shades, heights):
            lo, hi = c["allowed_from_deg"], c["allowed_to_deg"]
            plo, phi = c["published_allowed_deg"]
            bars.add(Line(ax2.c2p(lo, h_mine), ax2.c2p(hi, h_mine), color=color, stroke_width=4))
            bars.add(layout.label(
                f"{c['mass_earth']:g} M⊕, reproduced: {lo:.0f}° to {hi:.0f}°", font_size=15,
                color=color).next_to(ax2.c2p(180, h_mine), UP, buff=0.05))
            bars.add(Line(ax2.c2p(plo, h_paper), ax2.c2p(phi, h_paper), color=P.FG,
                          stroke_width=4))
            bars.add(layout.label(f"paper: {plo:.0f}° to {phi:.0f}°", font_size=15)
                     .next_to(ax2.c2p(phi, h_paper), RIGHT, buff=0.1))
        nine = nines[0]
        cap4 = layout.caption(
            f"Planet Nine, {nine['a_au']:.0f} AU and e = {nine['e']:.2f}: "
            "slow enough only on the far side of its orbit", font_size=22)
        self.play(FadeIn(ax2), run_time=0.9)
        self.play(Create(curves), FadeIn(limits), FadeIn(cap4), run_time=1.8)
        self.play(FadeIn(bars), run_time=0.8)
        timing.hold_to_read(self, cap4, bars, settle=1.2)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, "Saturn allows Planet Nine only near aphelion; Planets X and Y fare worse.")
