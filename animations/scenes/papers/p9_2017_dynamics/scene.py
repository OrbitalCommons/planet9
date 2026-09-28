"""Batygin & Morbidelli (2017) -- Dynamical evolution induced by Planet Nine.

Orbit-averaged (secular) theory puts the anti-aligned objects on an island of
apsidal libration -- but every orbit on that island crosses Planet Nine's path,
so pure secular theory predicts the whole cluster should be cleared out. The
paper resolves this with mean-motion resonances: an object locked in, say, the
3:1 resonance only crosses Planet Nine's orbit when the planet is far away.
Resonances take over wherever the objects' aphelia reach Planet Nine's
perihelion, beyond a period ratio P/P9 ~ 0.1.

The portrait is Planet Nine's exactly ring-averaged potential plus the giant
planets' J2 field; the orbit-crossing mask, the resonance chain and the
critical period ratio are computed in p9-2017-dynamics
(anim.json -> papers -> p9-2017-dynamics).
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
    Line,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2017-dynamics"

X0 = -90.0  # left edge of the Δϖ axis: both islands appear whole


def _wrap(dv):
    """Δϖ in degrees, wrapped into [X0, X0 + 360)."""
    return (np.asarray(dv) - X0) % 360.0 + X0


def _contour_segments(xs, ys, z, level):
    """Marching squares: line segments of z(y, x) = level in data coordinates."""
    segs = []
    for i in range(len(ys) - 1):
        for j in range(len(xs) - 1):
            c = [(xs[j], ys[i], z[i, j]), (xs[j + 1], ys[i], z[i, j + 1]),
                 (xs[j + 1], ys[i + 1], z[i + 1, j + 1]), (xs[j], ys[i + 1], z[i + 1, j])]
            pts = []
            for k in range(4):
                (x1, y1, z1), (x2, y2, z2) = c[k], c[(k + 1) % 4]
                if (z1 - level) * (z2 - level) < 0:
                    f = (level - z1) / (z2 - z1)
                    pts.append((x1 + f * (x2 - x1), y1 + f * (y2 - y1)))
            if len(pts) == 2:
                segs.append((pts[0], pts[1]))
            elif len(pts) == 4:
                segs.append((pts[0], pts[1]))
                segs.append((pts[2], pts[3]))
    return segs


def _chain(segs):
    """Join segments sharing endpoints into polylines."""
    key = lambda p: (round(p[0], 6), round(p[1], 6))  # noqa: E731
    ends = {}
    for n, (a, b) in enumerate(segs):
        ends.setdefault(key(a), []).append(n)
        ends.setdefault(key(b), []).append(n)
    used = [False] * len(segs)
    lines = []
    for n in range(len(segs)):
        if used[n]:
            continue
        used[n] = True
        line = [segs[n][0], segs[n][1]]
        for forward in (True, False):
            while True:
                tip = line[-1] if forward else line[0]
                nxt = [m for m in ends.get(key(tip), []) if not used[m]]
                if not nxt:
                    break
                m = nxt[0]
                used[m] = True
                a, b = segs[m]
                p = b if key(a) == key(tip) else a
                if forward:
                    line.append(p)
                else:
                    line.insert(0, p)
        lines.append(line)
    return lines


def _polylines(ax, lines, color, width=1.6, opacity=0.8):
    vm = VMobject(stroke_color=color, stroke_width=width, stroke_opacity=opacity)
    for line in lines:
        pts = [ax.c2p(x, y) for x, y in line]
        vm.start_new_path(pts[0])
        vm.add_points_as_corners(pts[1:])
    return vm


class Dynamics2017(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        por = d["portrait"]
        a9, e9, q9 = d["a9"], d["e9"], d["q9"]

        self.add(paper.scene_header(CRATE))

        # ---- 1. the secular portrait and the orbit-crossing zone -------------
        e = np.asarray(por["e"])
        dv = _wrap(por["dvarpi_deg"])
        order = np.argsort(dv)
        xs = np.append(dv[order], dv[order][0] + 360.0)
        h = np.asarray(por["h"])[:, order]
        h = np.hstack([h, h[:, :1]])
        cross = np.asarray(por["crossing"])[:, order]

        ax, labels = widgets.labeled_axes(
            [X0, X0 + 360, 90], [0, 1.0, 0.2], x_label="Δϖ = ϖ − ϖ₉  (deg, perihelion relative to Planet Nine)",
            y_label="eccentricity e", y_rotate=True, numbers=True,
            x_length=9.4, y_length=4.5, shift_down=-0.05)
        ax.shift(RIGHT * 0.9)
        labels.shift(RIGHT * 0.9)
        self.play(Create(ax), FadeIn(labels), run_time=1.0)

        de = e[1] - e[0]
        zone = VGroup()
        for i, ev in enumerate(e):
            j = 0
            while j < len(cross[i]):
                if not cross[i][j]:
                    j += 1
                    continue
                k = j
                while k + 1 < len(cross[i]) and cross[i][k + 1]:
                    k += 1
                x_lo = xs[j] - 2.5
                x_hi = min(xs[k] + 2.5, X0 + 360)
                p0, p1 = ax.c2p(max(x_lo, X0), ev - de / 2), ax.c2p(x_hi, ev + de / 2)
                r = Rectangle(width=p1[0] - p0[0], height=p1[1] - p0[1], stroke_width=0)
                zone.add(r.set_fill(P.RED, opacity=0.16).move_to((p0 + p1) / 2))
                j = k + 1

        band = (e > 0.05) & (e < 0.93)
        lo, hi = np.percentile(h[band], 4), np.max(h[band])
        levels = np.linspace(lo, hi, 16)[1:-1]
        contours = VGroup(*[
            _polylines(ax, _chain(_contour_segments(xs, e, h, lv)), P.FG, width=1.4, opacity=0.55)
            for lv in levels])
        cap = layout.caption("Orbit-averaged theory: each curve is a path an orbit's (Δϖ, e) follows",
                             font_size=21)
        self.play(Create(contours, lag_ratio=0.05), FadeIn(cap), run_time=2.2)
        timing.hold_to_read(self, cap, settle=0.6)

        # the anti-aligned island: the maximum of H on the Δϖ = 180° column
        j180 = int(np.argmin(np.abs(xs - 180.0)))
        i_top = int(np.argmax(np.where(band, h[:, j180], -np.inf)))
        centre = (180.0, float(e[i_top]))
        lv = h[i_top, j180] - 0.28 * (h[i_top, j180] - levels[-4])
        loops = _chain(_contour_segments(xs, e, h, lv))
        loop = min(loops, key=lambda ln: np.hypot(np.mean([p[0] for p in ln]) - centre[0],
                                                   60 * (np.mean([p[1] for p in ln]) - centre[1])))
        path = [ax.c2p(x, y) for x, y in loop]
        seg = np.cumsum([0.0] + [np.linalg.norm(np.subtract(path[k + 1], path[k]))
                                 for k in range(len(path) - 1)])

        def along(s):
            s = (s % 1.0) * seg[-1]
            k = int(np.clip(np.searchsorted(seg, s) - 1, 0, len(path) - 2))
            f = (s - seg[k]) / max(seg[k + 1] - seg[k], 1e-9)
            return np.asarray(path[k]) + f * (np.asarray(path[k + 1]) - np.asarray(path[k]))

        t = ValueTracker(0.0)
        body = always_redraw(lambda: Dot(along(t.get_value()), radius=0.08, color=P.GREEN))
        tag = layout.label("anti-aligned island", font_size=16, color=P.GREEN)
        tag.next_to(ax.c2p(180, 1.0), UP, buff=0.12)
        j0 = int(np.argmin(np.abs(xs - 0.0)))
        i_al = int(np.argmax(np.where(band, h[:, j0], -np.inf)))
        q_al = por["a"] * (1.0 - e[i_al])
        tag_al = layout.label(f"aligned island: q ≈ {q_al:.0f} AU", font_size=16, color=P.FG)
        tag_al.next_to(ax.c2p(0, 1.0), UP, buff=0.12)
        cap2 = layout.caption("Anti-aligned orbits circle the island around Δϖ = 180°: the clustering",
                              font_size=21)
        self.play(FadeOut(cap), FadeIn(cap2), FadeIn(tag), FadeIn(tag_al), FadeIn(body))
        self.play(t.animate.set_value(1.0), run_time=3.0, rate_func=rate_functions.linear)

        key = VGroup(*[layout.label(t, font_size=16, color=P.RED)
                       for t in ("red:", "the orbit", "crosses", "Planet Nine's", "orbit")]
                     ).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        key.next_to(labels[1], LEFT, buff=0.3)
        cap3 = layout.caption("But every orbit on that island crosses Planet Nine: secular theory "
                              "says the cluster should be destroyed", font_size=20)
        self.play(FadeIn(zone), FadeIn(key), FadeOut(cap2), FadeIn(cap3), run_time=1.2)
        self.play(t.animate.set_value(2.0), run_time=3.0, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap3, settle=0.4)
        self.remove(body)
        self.play(FadeOut(VGroup(ax, labels, zone, contours, tag, tag_al, key, cap3)), run_time=0.8)

        # ---- 2. resonance: crossing orbits that never meet ------------------
        chain = d["chain"]
        res = next(c for c in chain if c["j"] == 3)
        s = 3.2 / a9
        sun_at = np.array([2.35, 0.1, 0.0])
        sun = orbits.sun(0.1).move_to(sun_at)
        p9_orbit = orbits.ellipse_orbit(a9 * s, e9, color=P.BLUE, varpi=0.0).shift(sun_at)
        kbo_orbit = orbits.ellipse_orbit(res["a"] * s, res["e"], color=P.GREEN, varpi=np.pi,
                                         stroke_width=2.4).shift(sun_at)
        nep = orbits.ellipse_orbit(30.1 * s, 0.0, color=P.MUTED, stroke_width=1.2).shift(sun_at)
        l9 = VGroup(layout.label("Planet Nine", font_size=18, color=P.BLUE, weight="BOLD"),
                    layout.label(f"a = {a9:.0f} AU,  e = {e9:.1f}", font_size=16, color=P.BLUE)
                    ).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        lk = VGroup(layout.label("3:1 resonant object", font_size=18, color=P.GREEN, weight="BOLD"),
                    layout.label(f"a = {res['a']:.0f} AU,  q = {d['q_peri']:.0f} AU", font_size=16,
                                 color=P.GREEN)).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        VGroup(l9, lk).arrange(DOWN, buff=0.4, aligned_edge=LEFT).to_edge(LEFT, buff=0.6).shift(UP * 1.2)
        self.play(Create(p9_orbit), FadeIn(sun), FadeIn(nep), FadeIn(l9), run_time=1.0)
        self.play(Create(kbo_orbit), FadeIn(lk), run_time=0.9)

        tau = ValueTracker(0.0)

        def pos(a_au, ecc, varpi, mean):
            nu = orbits.nu_from_mean_anomaly(ecc, mean)
            return sun_at + orbits.orbit_point(a_au * s, ecc, nu, varpi)

        p9_dot = always_redraw(lambda: Dot(pos(a9, e9, 0.0, np.pi + 2 * np.pi * tau.get_value()),
                                           radius=0.11, color=P.BLUE))
        kbo_dot = always_redraw(lambda: Dot(
            pos(res["a"], res["e"], np.pi, 6 * np.pi * tau.get_value()), radius=0.07, color=P.GREEN))
        cap4 = layout.caption("Three laps per Planet Nine lap: it crosses the planet's path "
                              "only while the planet is far away", font_size=20)
        self.add(p9_dot, kbo_dot)
        self.play(FadeIn(cap4))
        self.play(tau.animate.set_value(1.0), run_time=6.0, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap4, settle=0.2)
        self.remove(p9_dot, kbo_dot)
        self.play(FadeOut(VGroup(p9_orbit, kbo_orbit, nep, sun, l9, lk, cap4)), run_time=0.7)

        # ---- 3. where resonances take over ---------------------------------
        crit = d["critical_ratio"]
        ax2, labels2 = widgets.labeled_axes(
            [0, 0.55, 0.1], [0, 900, 200], x_label="orbital period ratio P / P₉",
            y_label="aphelion distance (AU)", y_rotate=True, numbers=True,
            x_length=8.0, y_length=4.1, shift_down=-0.4)
        ax2.shift(LEFT * 1.9)
        labels2.shift(LEFT * 1.9)
        self.play(Create(ax2), FadeIn(labels2), run_time=0.9)
        q9_line = DashedLine(ax2.c2p(0, q9), ax2.c2p(0.55, q9), color=P.BLUE, stroke_width=2)
        q9_lab = VGroup(layout.label("Planet Nine's", font_size=15, color=P.BLUE),
                        layout.label(f"perihelion, {q9:.0f} AU", font_size=15, color=P.BLUE)
                        ).arrange(DOWN, buff=0.06)
        q9_lab.next_to(ax2.c2p(0.415, q9), UP, buff=0.1)

        stems = VGroup()
        for c in chain:
            x = c["period_ratio"]
            if x > 0.55:
                continue
            col = P.ORANGE if c["reaches_p9"] else P.FG
            stems.add(VGroup(Line(ax2.c2p(x, 0), ax2.c2p(x, c["aphelion"]), color=col,
                                  stroke_width=3, stroke_opacity=0.8),
                             Dot(ax2.c2p(x, c["aphelion"]), radius=0.055, color=col)))
        names = VGroup(*[
            layout.label(f"{c['j']}:1", font_size=13, color=P.ORANGE)
            .next_to(ax2.c2p(c["period_ratio"], c["aphelion"]), UP, buff=0.08)
            for c in chain if c["j"] in (2, 3, 4, 5)])
        cap5 = layout.caption(f"The N:1 resonances for objects with q = {d['q_peri']:.0f} AU: "
                              "how far out each one's orbit reaches", font_size=20)
        self.play(Create(stems, lag_ratio=0.08), FadeIn(names), FadeIn(cap5), run_time=1.6)
        self.play(Create(q9_line), FadeIn(q9_lab))
        timing.hold_to_read(self, cap5, settle=0.4)

        mark = DashedLine(ax2.c2p(crit, 0), ax2.c2p(crit, 900), color=P.ORANGE, stroke_width=2)
        legend = VGroup(
            layout.label("reaches Planet Nine: resonance-protected", font_size=15, color=P.ORANGE),
            layout.label("stays inside: secular, isolated", font_size=15, color=P.FG),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        legend.next_to(ax2.c2p(0.2, 880), RIGHT, buff=0.0)
        readout = paper.result_readout("secular → resonant at", f"P/P₉ = {crit:.2f}",
                                       color=P.ORANGE).scale(0.8)
        pub = layout.label("paper: P/P₉ ≳ 0.1", font_size=17, color=P.FG)
        alt = VGroup(layout.label("for a₉ = 600 AU, e₉ = 0.5:", font_size=15, color=P.FG),
                     layout.label(f"{d['critical_ratio_600']:.2f}   (paper ≳ 0.15)", font_size=15,
                                  color=P.FG)).arrange(DOWN, buff=0.08)
        side = VGroup(readout, pub, alt).arrange(DOWN, buff=0.18)
        side.next_to(ax2, RIGHT, buff=0.5)
        cap6 = layout.caption("Beyond P/P₉ ≈ 0.1 every object reaches Planet Nine's orbit "
                              "and survives only inside a resonance", font_size=20)
        self.play(Create(mark), FadeIn(legend), FadeOut(cap5), FadeIn(cap6), run_time=1.0)
        self.play(FadeIn(side, shift=LEFT * 0.1))
        timing.hold_to_read(self, cap6, side, settle=1.0)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "Secular theory draws the island; resonances keep it alive past P/P₉ ≈ 0.1.")
