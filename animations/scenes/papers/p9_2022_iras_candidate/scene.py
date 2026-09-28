"""Rowan-Robinson (2021) -- a search for Planet 9 in the IRAS data.

IRAS scanned each patch of sky on passes weeks to months apart. A body a few
hundred AU away shifts by arcminutes between passes because the Earth itself
moves, so Planet Nine would appear as a 60 µm source seen on some passes and
a lone detection nearby on another. Of several hundred such pairings one
survives: a 0.57 Jy source that moved 20 arcmin in 12 weeks. Reproduced in
p9-2022-iras-candidate: the chance-pairing estimate, the flux-distance curves
and the predicted population are the crate's own
(anim.json -> papers -> p9-2022-iras-candidate).
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
    Dot,
    Ellipse,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    Star,
    SurroundingRectangle,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2022-iras-candidate"


class Plot(VGroup):
    """Axes with the origin at the lower-left corner of the data range, optional
    log10 scales, and tick labels and axis titles in the house style."""

    def __init__(self, x_range, y_range, x_ticks, y_ticks, x_label, y_label, width=10.0,
                 height=4.2, centre=(0.35, 0.55), x_log=False, y_log=False, x_fmt="{:g}",
                 y_fmt="{:g}", **kwargs):
        super().__init__(**kwargs)
        self.x_log, self.y_log = x_log, y_log
        self.x0, self.x1 = (self._s(v, x_log) for v in x_range)
        self.y0, self.y1 = (self._s(v, y_log) for v in y_range)
        self.left, self.bottom = centre[0] - width / 2, centre[1] - height / 2
        self.w, self.h = width, height
        frame = VGroup(
            Line([self.left, self.bottom, 0], [self.left + width, self.bottom, 0]),
            Line([self.left, self.bottom, 0], [self.left, self.bottom + height, 0]),
        ).set_stroke(P.MUTED, width=2)
        marks = VGroup()
        for v in x_ticks:
            at = self.p(v, y_range[0])
            marks.add(Line(at, at + DOWN * 0.08, color=P.MUTED, stroke_width=1.5))
            marks.add(layout.label(self._text(x_fmt, v), font_size=15, color=P.MUTED)
                      .next_to(at, DOWN, buff=0.14))
        widest = 0.0
        for v in y_ticks:
            at = self.p(x_range[0], v)
            lab = layout.label(self._text(y_fmt, v), font_size=15, color=P.MUTED)
            lab.next_to(at, LEFT, buff=0.14)
            widest = max(widest, lab.width)
            marks.add(Line(at, at + LEFT * 0.08, color=P.MUTED, stroke_width=1.5), lab)
        xl = layout.label(x_label, font_size=19, color=P.FG)
        xl.move_to([centre[0], self.bottom - 0.66, 0])
        yl = layout.label(y_label, font_size=17, color=P.FG).rotate(np.pi / 2)
        yl.move_to([self.left - widest - 0.45, centre[1], 0])
        self.add(frame, marks, xl, yl)

    @staticmethod
    def _s(v, log):
        return float(np.log10(v)) if log else float(v)

    @staticmethod
    def _text(fmt, v):
        return fmt(v) if callable(fmt) else fmt.format(v)

    def p(self, x, y):
        """Scene point of the data point (x, y), clamped to the plot area."""
        u = (self._s(x, self.x_log) - self.x0) / (self.x1 - self.x0)
        v = (self._s(y, self.y_log) - self.y0) / (self.y1 - self.y0)
        return np.array([self.left + self.w * min(max(u, 0.0), 1.0),
                         self.bottom + self.h * min(max(v, 0.0), 1.0), 0.0])

    def inside(self, x, y):
        sx, sy = self._s(x, self.x_log), self._s(y, self.y_log)
        return (min(self.x0, self.x1) <= sx <= max(self.x0, self.x1)
                and min(self.y0, self.y1) <= sy <= max(self.y0, self.y1))

    def curve(self, xs, ys, color, stroke_width=3.0):
        """Polyline through the data points that fall inside the plot area."""
        runs, cur = [], []
        for x, y in zip(xs, ys):
            if self.inside(x, y):
                cur.append(self.p(x, y))
            elif cur:
                runs.append(cur)
                cur = []
        if cur:
            runs.append(cur)
        g = VGroup()
        for pts in runs:
            if len(pts) > 1:
                g.add(VMobject(color=color, stroke_width=stroke_width)
                      .set_points_as_corners(pts))
        return g

    def hline(self, y, color, x_from=None, x_to=None, stroke_width=1.6):
        a = self.p(self._x(x_from, self.x0), y)
        b = self.p(self._x(x_to, self.x1), y)
        return DashedLine(a, b, color=color, stroke_width=stroke_width)

    def vline(self, x, color, y_from=None, y_to=None, stroke_width=1.6):
        a = self.p(x, self._y(y_from, self.y0))
        b = self.p(x, self._y(y_to, self.y1))
        return DashedLine(a, b, color=color, stroke_width=stroke_width)

    def band(self, x_from, x_to, y_from, y_to, color, opacity=0.16):
        a, b = self.p(x_from, y_from), self.p(x_to, y_to)
        r = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=color, stroke_width=0)
        return r.set_fill(color, opacity=opacity)

    def _x(self, v, scaled):
        return (10 ** scaled if self.x_log else scaled) if v is None else v

    def _y(self, v, scaled):
        return (10 ** scaled if self.y_log else scaled) if v is None else v


def _readout(title, value, note, colour):
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


class IrasCandidate2022(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        c = d["candidate"]
        self.add(paper.scene_header(CRATE))

        # 1. why a distant planet shows up as a pair of mismatched sources
        per_arcmin = 0.105
        a = d["parallax_semi_major_arcmin"]
        b = a * np.sin(np.radians(c["ecliptic_lat_deg"]))
        centre = np.array([-2.4, 0.45, 0])
        ell = Ellipse(width=2 * a * per_arcmin, height=2 * b * per_arcmin, color=P.MUTED,
                      stroke_width=1.6).move_to(centre)
        ell_lab = layout.label(f"its yearly parallax loop at {d['published_distance_au']:.0f} AU",
                               font_size=16, color=P.MUTED)
        ell_lab.next_to(ell, DOWN, buff=0.45)
        # sweep between the passes: the chord of the loop equals the measured motion
        sweep = 2 * np.arcsin(min(1.0, c["motion_arcmin"] / (2 * a)))
        th0 = np.radians(200.0)

        def on_loop(th):
            return centre + per_arcmin * np.array([a * np.cos(th), b * np.sin(th), 0])

        p1, p3 = on_loop(th0), on_loop(th0 + sweep)
        d1 = Dot(p1, radius=0.1, color=P.GREEN)
        d3 = Dot(p3, radius=0.1, color=P.GREEN)
        l1 = layout.label("passes 1 and 2", font_size=16, color=P.GREEN).next_to(d1, LEFT, buff=0.15)
        l3 = layout.label("pass 3", font_size=16, color=P.GREEN).next_to(d3, RIGHT, buff=0.15)
        hop = Arrow(p1, p3, buff=0.12, color=P.ORANGE, stroke_width=3)
        hop_lab = layout.label(f"{c['motion_arcmin']:.0f} arcmin in {c['motion_weeks']:.0f} weeks",
                               font_size=17, color=P.ORANGE)
        hop_lab.next_to(hop.get_center(), UP + RIGHT, buff=0.15)
        moon = Circle(radius=0.5 * 31.0 * per_arcmin, color=P.FG, stroke_width=1.4)
        moon.set_fill(P.FG, opacity=0.08).move_to([3.4, 0.45, 0])
        moon_lab = layout.label("the full Moon, same scale", font_size=16, color=P.FG)
        moon_lab.next_to(moon, DOWN, buff=0.45)
        cap = layout.caption("As the Earth circles the Sun, a distant body traces a small loop",
                             font_size=22)
        self.play(Create(ell), FadeIn(ell_lab), FadeIn(moon), FadeIn(moon_lab), FadeIn(cap),
                  run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.2)
        cap2 = layout.caption("IRAS would log it in one place, then as a lone source nearby",
                              font_size=22)
        self.play(FadeIn(d1), FadeIn(l1), FadeOut(cap), FadeIn(cap2))
        self.play(Create(hop), FadeIn(d3), FadeIn(l3), FadeIn(hop_lab), run_time=1.2)
        timing.hold_to_read(self, cap2, settle=0.6)
        self.play(FadeOut(VGroup(ell, ell_lab, d1, d3, l1, l3, hop, hop_lab, moon, moon_lab,
                                 cap2)))

        # 2. hundreds of such pairings, one survivor
        n = int(d["published_associations"])
        cols = 38
        dots = VGroup(*[Dot(radius=0.045, color=P.GREEN) for _ in range(n)])
        dots.arrange_in_grid(cols=cols, buff=0.11).move_to(UP * 0.55)
        count = layout.label(
            f"{n} pairings examined by eye ({d['published_pairs']:.0f} pairs, "
            f"{d['published_triplets']:.0f} triplets, {d['published_close_pairs']:.0f} close pairs)",
            font_size=18, color=P.FG)
        count.next_to(dots, UP, buff=0.3)
        cap3 = layout.caption(
            f"Chance alone predicts {d['chance_associations']:,.0f} coincidences: "
            "most pairs are unrelated sources", font_size=22)
        self.play(FadeIn(dots, lag_ratio=0.002), FadeIn(count), FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, settle=0.3)
        keep = dots[n // 2 + cols // 2]
        others = VGroup(*[x for x in dots if x is not keep])
        cap4 = layout.caption(
            "Checking every pass in the raw scans leaves a single candidate", font_size=22)
        self.play(others.animate.set_color(P.RED).set_opacity(0.25),
                  keep.animate.scale(2.2).set_color(P.GREEN), FadeOut(cap3), FadeIn(cap4),
                  run_time=1.6)
        timing.hold_to_read(self, cap4, settle=0.6)
        self.play(FadeOut(VGroup(dots, count, cap4)))

        # 3. how far away a 0.57 Jy source would be
        dist = np.array(d["distance_au"])
        plot = Plot([100, 450], [0.05, 3.0], [100, 150, 200, 250, 300, 350, 400, 450],
                    [0.1, 0.3, 1, 3], "distance from the Sun (AU)", "60 µm flux density (Jy)",
                    y_log=True, centre=(-0.4, 0.5), width=8.8, height=4.3)
        curves = VGroup()
        for k, cv in enumerate(d["curves"]):
            curves.add(plot.curve(dist, cv["flux_60um_jy"], P.BLUE,
                                  stroke_width=[1.8, 3.0, 1.8][k]))
        cl = layout.label(f"{d['published_mass_lo']:.0f}-{d['published_mass_hi']:.0f} Earth "
                          f"masses at {d['model_temp_k']:.0f} K", font_size=16, color=P.BLUE)
        cl.next_to(plot.p(330, 0.2), UP + RIGHT, buff=0.1)
        flux = plot.hline(c["flux_60um_jy"], P.GREEN)
        flux_lab = layout.label(f"candidate {c['flux_60um_jy']:.2f} Jy", font_size=16,
                                color=P.GREEN)
        flux_lab.next_to(plot.p(450, c["flux_60um_jy"]), UP + LEFT, buff=0.08)
        pub = plot.band(d["published_distance_au"] - d["published_distance_err_au"],
                        d["published_distance_au"] + d["published_distance_err_au"], 0.05, 3.0,
                        P.FG, opacity=0.1)
        pub_lab = layout.label(f"paper: {d['published_distance_au']:.0f} ± "
                               f"{d['published_distance_err_au']:.0f} AU", font_size=15,
                               color=P.FG)
        pub_lab.next_to(plot.p(d["published_distance_au"], 3.0), DOWN, buff=0.12)
        self.play(FadeIn(plot))
        cap5 = layout.caption("A body's heat fades as distance squared", font_size=22)
        self.play(Create(curves), FadeIn(cl), FadeIn(cap5), run_time=1.4)
        timing.hold_to_read(self, cap5, settle=0.2)
        hits = VGroup(*[Dot(plot.p(cv["implied_distance_au"], c["flux_60um_jy"]), radius=0.07,
                            color=P.GREEN) for cv in d["curves"]])
        lo = min(cv["implied_distance_au"] for cv in d["curves"])
        hi = max(cv["implied_distance_au"] for cv in d["curves"])
        box = _readout("implied distance", f"{lo:.0f}-{hi:.0f} AU",
                       f"paper: {d['published_distance_au']:.0f} ± "
                       f"{d['published_distance_err_au']:.0f} AU", P.GREEN)
        box.to_edge(RIGHT, buff=0.35).shift(UP * 1.2)
        temps = [cv["temperature_at_published_distance_k"] for cv in d["curves"]]
        cap6 = layout.caption(
            f"At {d['model_temp_k']:.0f} K it lies nearer than the paper's fit; "
            f"{min(temps):.0f}-{max(temps):.0f} K would put it at {d['published_distance_au']:.0f} AU",
            font_size=22)
        self.play(Create(flux), FadeIn(flux_lab), FadeIn(hits), FadeIn(pub), FadeIn(pub_lab),
                  FadeIn(box), FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, box, settle=0.6)
        self.play(FadeOut(VGroup(plot, curves, cl, flux, flux_lab, hits, pub, pub_lab, box,
                                 cap6)))

        # 4. where it is, against where Planet Nine is predicted to be
        m = sky.SkyMap(width=11.0, dec_range=(-80, 80), centre=(0.0, 0.5, 0.0))
        ecl, gal = m.reference_curves()
        pop = m.dots(d["population"], color=P.TEAL, radius=0.028)
        star = Star(n=5, outer_radius=0.17, color=P.GREEN).set_fill(P.GREEN, 1.0)
        star.move_to(m.p(c["ra_deg"], c["dec_deg"]))
        star_lab = layout.label(
            f"the IRAS candidate: {c['ecliptic_lat_deg']:.0f}° from the ecliptic", font_size=16,
            color=P.GREEN)
        star_lab.next_to(star, RIGHT, buff=0.15)
        self.play(FadeIn(m), Create(ecl), Create(gal), run_time=1.0)
        cap7 = layout.caption(
            f"Predicted Planet Nines stay within {d['population_max_ecliptic_lat_deg']:.0f}° "
            "of the ecliptic", font_size=22)
        self.play(FadeIn(pop, lag_ratio=0.02), FadeIn(cap7), run_time=1.2)
        timing.hold_to_read(self, cap7, settle=0.2)
        cap8 = layout.caption("The candidate is far off every predicted orbit, and unconfirmed",
                              font_size=22)
        self.play(FadeIn(star, scale=2.0), FadeIn(star_lab), FadeOut(cap7), FadeIn(cap8))
        timing.hold_to_read(self, cap8, settle=1.0)
        self.play(FadeOut(cap8))

        layout.show_takeaway(
            self, "One faint IRAS mover survives, but not where Planet Nine should be.")
