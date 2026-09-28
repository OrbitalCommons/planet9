"""Holman & Payne (2016) -- Cassini range observations as a finder chart.

The paper fits the Cassini range residuals with a tidal perturber that is free
in direction and strength, and turns the ranging into a preferred patch of sky
and a preferred tide (mass over distance cubed). The reproduction crate lays
the Batygin & Brown orbit on the sky and asks, point by point, whether a planet
there would push the Earth-Saturn range past the ranging precision.

Everything drawn comes from anim.json -> papers -> p9-2016-holman-payne: the
orbit's sky track with its range signal, the reproduced favoured position, the
heaviest planet that stays under the precision at each distance, and the
paper's region and tidal range as labelled comparisons.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2016-holman-payne"


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


def runs_by_flag(track):
    """Split the closed sky track into runs of constant verdict; neighbouring
    runs share their boundary point so the drawn curve has no gaps."""
    n = len(track)
    start = next(k for k in range(n) if track[k]["excluded"] != track[k - 1]["excluded"])
    out, cur = [], [track[start]]
    for j in range(1, n + 1):
        pt = track[(start + j) % n]
        if pt["excluded"] != cur[-1]["excluded"]:
            out.append((cur[-1]["excluded"], cur + [pt]))
            cur = [pt]
        else:
            cur.append(pt)
    return out


class HolmanPayne2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        orbit, pub, fav = d["orbit"], d["published"], d["favored"]
        track = d["track"]

        self.add(paper.scene_header(CRATE))

        # 1. the reference orbit on the sky, judged by the ranging
        m = sky.SkyMap(width=10.0, dec_range=(-75, 75), centre=(-1.55, 0.35, 0.0))
        ecl, gal = m.reference_curves()
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), run_time=1.0)

        lines = VGroup()
        for excluded, pts in runs_by_flag(track):
            lines.add(m.polyline([(p["ra_deg"], p["dec_deg"]) for p in pts],
                                 color=P.RED if excluded else P.TEAL, stroke_width=4.0))

        def swatch(color, text, dashed=False):
            mark = Line(LEFT * 0.2, RIGHT * 0.2, color=color, stroke_width=4)
            return VGroup(mark, layout.label(text, font_size=14, color=color, line_spacing=0.9)
                          ).arrange(RIGHT, buff=0.15)

        key = VGroup(
            swatch(P.RED, f"range signal\nabove {d['precision_m']:.0f} m"),
            swatch(P.TEAL, "range signal\nbelow it"),
            swatch(P.ORANGE, "ecliptic"),
            swatch(P.PURPLE, "galactic plane"),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        key.move_to([5.1, 1.8, 0])
        cap = layout.caption(
            f"Reproduced: a {orbit['mass_earth']:.0f} M⊕ planet at each point of the "
            f"{orbit['a_au']:.0f} AU orbit, judged by Cassini's range", font_size=22)
        self.play(Create(lines), FadeIn(key), FadeIn(cap), run_time=2.0)
        timing.hold_to_read(self, cap, key, settle=0.6)

        # 1b. walk a planet along the track and read what Cassini would see
        k = ValueTracker(0.0)

        def here():
            return track[int(round(k.get_value())) % len(track)]

        walker = always_redraw(lambda: Dot(m.p(here()["ra_deg"], here()["dec_deg"]), radius=0.1,
                                           color=P.BLUE).set_z_index(6))

        def meter():
            p = here()
            col = P.RED if p["excluded"] else P.TEAL
            g = VGroup(
                layout.label("a planet here:", font_size=16, color=P.MUTED),
                layout.label(f"{p['r_au']:.0f} AU from the Sun", font_size=17),
                layout.label(f"range signal {p['residual_m']:.0f} m", font_size=20, color=col,
                             weight="BOLD"),
                layout.label("Cassini would see it" if p["excluded"] else "hidden in the noise",
                             font_size=16, color=col),
            ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
            return g.move_to([5.35, -1.6, 0])

        readout = always_redraw(meter)
        capw = layout.caption("Near perihelion the tug is loud; far out it drops under the noise",
                              font_size=22)
        self.play(FadeIn(walker), FadeIn(readout), FadeOut(cap), FadeIn(capw), run_time=0.6)
        self.play(k.animate.set_value(len(track) - 1), run_time=8.0, rate_func=lambda t: t)
        timing.hold_to_read(self, capw, settle=0.4)
        self.play(FadeOut(walker), FadeOut(readout), run_time=0.5)
        cap = capw

        # 2. the paper's preferred region and the reproduced favoured spot
        h = pub["half_extent_deg"]
        box = m.box(pub["ra_deg"] - h, pub["ra_deg"] + h, pub["dec_deg"] - h,
                    pub["dec_deg"] + h, color=P.GREEN, opacity=0.22, stroke_width=2.0)
        spot = Dot(m.p(fav["ra_deg"], fav["dec_deg"]), radius=0.09, color=P.ORANGE).set_z_index(4)
        found = VGroup(
            VGroup(Polygon([-0.2, -0.12, 0], [0.2, -0.12, 0], [0.2, 0.12, 0], [-0.2, 0.12, 0],
                           color=P.GREEN, stroke_width=2).set_fill(P.GREEN, opacity=0.22),
                   layout.label(
                       f"paper: RA {pub['ra_deg']:.0f}°,\nDec {pub['dec_deg']:.0f}°, ±{h:.0f}°",
                       font_size=14, color=P.GREEN, line_spacing=0.9)).arrange(RIGHT, buff=0.15),
            VGroup(Dot(radius=0.09, color=P.ORANGE).shift(RIGHT * 0.11),
                   layout.label(
                       f"reproduced: RA {fav['ra_deg']:.0f}°,\nDec {fav['dec_deg']:.0f}°",
                       font_size=14, color=P.ORANGE, line_spacing=0.9)).arrange(RIGHT, buff=0.26),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        found.next_to(key, DOWN, buff=0.4, aligned_edge=LEFT)
        cap2 = layout.caption(
            f"The favoured spot lands {d['offset_deg']:.0f}° from the centre of the "
            "paper's preferred region", font_size=22)
        self.play(FadeIn(box), FadeIn(found[0]), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        self.play(FadeIn(spot), FadeIn(found[1]), run_time=0.7)
        timing.hold_to_read(self, cap2, found, settle=1.2)
        self.play(FadeOut(VGroup(m, ecl, gal, lines, key, box, spot, found, cap2)))

        # 3. how heavy, how far
        out = [p for p in track if p["nu_deg"] <= 180.0]
        r = np.array([p["r_au"] for p in out])
        ceiling = np.array([p["max_mass_earth"] for p in out])
        x0, x1, y0, y1 = 250.0, 1150.0, 0.3, 100.0
        ax = Plot(
            (x0, x1), (y0, y1),
            {v: f"{v}" for v in range(300, 1101, 200)},
            {0.3: "0.3", 1: "1", 3: "3", 10: "10", 30: "30", 100: "100"},
            "distance from the Sun along the orbit  (AU)", "planet mass  (M⊕)",
            size=(10.0, 4.1), centre=(0.2, 0.55), ylog=True)

        rr = np.linspace(x0, x1, 60)
        unit_m, unit_r = pub["tidal_unit_mass_earth"], pub["tidal_unit_distance_au"]
        t_lo, t_hi = pub["tidal_range"]
        edge_lo = np.clip(unit_m * t_lo * (rr / unit_r) ** 3, y0, y1)
        edge_hi = np.clip(unit_m * t_hi * (rr / unit_r) ** 3, y0, y1)
        wedge = Polygon(*[ax.c2p(x, y) for x, y in zip(rr, edge_lo)],
                        *[ax.c2p(x, y) for x, y in zip(rr[::-1], edge_hi[::-1])],
                        stroke_width=0).set_fill(P.GREEN, opacity=0.3)
        wedge_lab = layout.label(
            f"paper: favoured tide, {t_lo:g} to {t_hi:g} times\n"
            f"that of {unit_m:.0f} M⊕ at {unit_r:.0f} AU",
            font_size=16, color=P.GREEN, line_spacing=0.9)
        wedge_lab.move_to(ax.c2p(930, 0.75))

        keep = ceiling >= y0
        curve = widgets.curve(ax, r[keep], ceiling[keep], color=P.RED, stroke_width=3.5)
        curve_lab = layout.label(
            f"reproduced: heaviest planet\nthat stays under {d['precision_m']:.0f} m",
            font_size=16, color=P.RED, line_spacing=0.9)
        curve_lab.move_to(ax.c2p(930, 3.4))
        nominal = Line(ax.c2p(x0, orbit["mass_earth"]), ax.c2p(x1, orbit["mass_earth"]),
                       color=P.BLUE, stroke_width=1.5).set_stroke(opacity=0.7)
        nominal_lab = layout.label(f"{orbit['mass_earth']:.0f} M⊕ Planet Nine", font_size=14,
                                   color=P.BLUE)
        nominal_lab.next_to(ax.c2p(x1, orbit["mass_earth"]), DOWN, buff=0.08, aligned_edge=RIGHT)

        cap3 = layout.caption(
            "A tide is mass over distance cubed: a closer planet must be lighter",
            font_size=22)
        self.play(FadeIn(ax), run_time=0.9)
        self.play(Create(nominal), FadeIn(nominal_lab), Create(curve), FadeIn(curve_lab),
                  FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, curve_lab, settle=0.8)
        cap4 = layout.caption(
            "The paper's fit prefers a tide at or above this ceiling: a hint, not a limit",
            font_size=22)
        self.play(FadeIn(wedge), FadeIn(wedge_lab), FadeOut(cap3), FadeIn(cap4), run_time=1.0)
        timing.hold_to_read(self, cap4, wedge_lab, settle=1.2)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, f"Cassini ranging becomes a finder chart: a {2 * h:.0f}°-wide patch of sky "
                  "to search.")
