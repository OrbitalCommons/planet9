"""Khain, Batygin & Brown (2018) -- the distant Kuiper belt from a broad
primordial perihelion distribution.

At a = 345 AU (the paper's Figure 4) the secular paths of Planet Nine's
potential are drawn in the (Δϖ, q) plane, with the region where an orbit
crosses Planet Nine's shaded. Anti-aligned orbits live inside that region and
need resonances to survive; aligned orbits survive only on the high-perihelion
paths that never cross it. The crate's own N-body integration (Neptune, Planet
Nine and the J2-averaged giants) then evolves a narrow (q = 30-36 AU) and a
broad (q = 30-300 AU) primordial disk for the first 60 Myr of the paper's
4 Gyr. The 4 Gyr survivor numbers are the paper's, shown as such. Data from
p9-2018-kuiper-belt (anim.json -> papers -> p9-2018-kuiper-belt).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    Cross,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2018-kuiper-belt"

X0 = -90.0  # left edge of the Δϖ axis: both islands appear whole
Q_TOP = 345.0


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


def _contours(ax, xs, ys, z, level, color, width=1.4, opacity=0.55):
    vm = VMobject(stroke_color=color, stroke_width=width, stroke_opacity=opacity)
    for (x1, y1), (x2, y2) in _contour_segments(xs, ys, z, level):
        vm.start_new_path(ax.c2p(x1, y1))
        vm.add_line_to(ax.c2p(x2, y2))
    return vm


def _track_pos(p, t):
    """(Δϖ, q, a) of a particle at time t (Myr), interpolated between
    snapshots; Δϖ is unwrapped before interpolating."""
    ts = np.asarray(p["t_myr"])
    dv = np.degrees(np.unwrap(np.radians(p["dvarpi_deg"])))
    x = float(np.interp(t, ts, dv))
    return _wrap(x), float(np.interp(t, ts, p["q"])), float(np.interp(t, ts, p["a"]))


class KuiperBelt2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        por = d["portrait"]
        a_p = por["a"]
        t_end = d["t_myr"]

        self.add(paper.scene_header(CRATE))

        # ---- 1. the phase plane in perihelion distance ----------------------
        e = np.asarray(por["e"])
        q = a_p * (1.0 - e)
        dv = _wrap(por["dvarpi_deg"])
        order = np.argsort(dv)
        xs = np.append(dv[order], dv[order][0] + 360.0)
        h = np.asarray(por["h"])[:, order]
        h = np.hstack([h, h[:, :1]])
        cross = np.asarray(por["crossing"])[:, order]
        rows = np.argsort(q)
        qs, h, cross = q[rows], h[rows], cross[rows]

        ax, labels = widgets.labeled_axes(
            [X0, X0 + 360, 90], [0, 350, 50], x_label="Δϖ  (deg, perihelion relative to Planet Nine)",
            y_label="perihelion q (AU)", y_rotate=True, numbers=True,
            x_length=8.0, y_length=4.3, shift_down=-0.4)
        ax.shift(LEFT * 1.9)
        labels.shift(LEFT * 1.9)
        nep = DashedLine(ax.c2p(X0, 30), ax.c2p(X0 + 360, 30), color=P.MUTED, stroke_width=2)
        nep_lab = layout.label("Neptune, 30 AU", font_size=14, color=P.MUTED)
        nep_lab.next_to(ax.c2p(X0 + 360, 30), RIGHT, buff=0.12)
        self.play(Create(ax), FadeIn(labels), Create(nep), FadeIn(nep_lab), run_time=1.0)

        dq = qs[1] - qs[0]
        zone = VGroup()
        for i, qv in enumerate(qs):
            j = 0
            while j < len(cross[i]):
                if not cross[i][j]:
                    j += 1
                    continue
                k = j
                while k + 1 < len(cross[i]) and cross[i][k + 1]:
                    k += 1
                p0 = ax.c2p(max(xs[j] - 2.5, X0), max(qv - dq / 2, 0))
                p1 = ax.c2p(min(xs[k] + 2.5, X0 + 360), min(qv + dq / 2, 350))
                r = Rectangle(width=p1[0] - p0[0], height=p1[1] - p0[1], stroke_width=0)
                zone.add(r.set_fill(P.RED, opacity=0.16).move_to((p0 + p1) / 2))
                j = k + 1

        band = (qs > 15) & (qs < 330)
        lo, hi = np.percentile(h[band], 4), np.max(h[band])
        levels = np.linspace(lo, hi, 16)[1:-1]
        contours = VGroup(*[_contours(ax, xs, qs, h, lv, P.FG) for lv in levels])
        cap = layout.caption(f"Secular paths at a = {a_p:.0f} AU, in perihelion distance q",
                             font_size=21)
        self.play(Create(contours, lag_ratio=0.05), FadeIn(cap), run_time=2.0)
        timing.hold_to_read(self, cap, settle=0.3)

        key = VGroup(layout.label("red: the orbit crosses", font_size=16, color=P.RED),
                     layout.label("Planet Nine's orbit", font_size=16, color=P.RED)
                     ).arrange(DOWN, buff=0.06, aligned_edge=LEFT)
        key.next_to(ax, RIGHT, buff=0.45).align_to(ax, UP)
        cap2 = layout.caption("Red: orbits that cross Planet Nine's path and are eventually "
                              "scattered away", font_size=21)
        self.play(FadeIn(zone), FadeIn(key), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        timing.hold_to_read(self, cap2, settle=0.3)

        # the aligned paths that never touch the red: all of them sit high up
        level, q_safe = d["aligned_safe_level"], d["aligned_safe_q_min"]
        near = np.abs(xs) < 90
        safe_levels = np.linspace(level, h[:, near].max(), 6)[:-1]
        safe_paths = VGroup()
        for lv in safe_levels:
            vm = VMobject(stroke_color=P.GREEN, stroke_width=2.4)
            for (x1, y1), (x2, y2) in _contour_segments(xs, qs, h, lv):
                if abs(x1) < 90 and abs(x2) < 90 and max(y1, y2) < 290:
                    vm.start_new_path(ax.c2p(x1, y1))
                    vm.add_line_to(ax.c2p(x2, y2))
            safe_paths.add(vm)
        floor = DashedLine(ax.c2p(X0, q_safe), ax.c2p(90, q_safe), color=P.GREEN, stroke_width=2)
        tag = layout.label("aligned, never crossing", font_size=16, color=P.GREEN, weight="BOLD")
        tag.next_to(ax.c2p(0, 350), UP, buff=0.1)
        safe_read = paper.result_readout("safe aligned orbits reach down to",
                                         f"q = {q_safe:.0f} AU", color=P.GREEN).scale(0.75)
        safe_read.next_to(key, DOWN, buff=0.45).align_to(key, LEFT)
        cap2b = layout.caption("Aligned orbits that never meet Planet Nine all keep a high perihelion",
                               font_size=21)
        self.play(Create(safe_paths), Create(floor), FadeIn(tag), FadeOut(cap2), FadeIn(cap2b),
                  run_time=1.4)
        self.play(FadeIn(safe_read))
        timing.hold_to_read(self, cap2b, safe_read, settle=0.6)
        self.play(FadeOut(safe_read), FadeOut(key), FadeOut(cap2b))

        # ---- 2. the crate's N-body: narrow vs broad primordial disks -------
        clock = ValueTracker(0.0)

        def particle(p, color):
            ts = p["t_myr"]
            t_gone = ts[-1] if not p["survived"] else None

            def draw():
                t = clock.get_value()
                if t_gone is not None and t > t_gone:
                    x, y, _ = _track_pos(p, t_gone)
                    return Cross(scale_factor=0.09, stroke_color=P.RED, stroke_width=3).move_to(
                        ax.c2p(x, min(y, 345)))
                x, y, a = _track_pos(p, t)
                off = abs(a - a_p) > 60 or y > 345
                dot = Dot(ax.c2p(x, min(y, 345)), radius=0.065, color=color)
                return dot.set_opacity(0.3 if off else 1.0)

            return always_redraw(draw)

        narrow = VGroup(*[particle(p, P.ORANGE) for p in d["narrow"]])
        broad = VGroup(*[particle(p, P.TEAL) for p in d["broad"]])
        nq, bq = d["narrow_q"], d["broad_q"]
        legend = VGroup(
            layout.label(f"narrow disk: q = {nq[0]:.0f}–{nq[1]:.0f} AU", font_size=16, color=P.ORANGE),
            layout.label(f"broad disk: q = {bq[0]:.0f}–{bq[1]:.0f} AU", font_size=16, color=P.TEAL),
            layout.label("red ×: removed from the system", font_size=16, color=P.RED),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        legend.next_to(ax, RIGHT, buff=0.45).align_to(ax, UP)
        t_lab = always_redraw(lambda: layout.label(
            f"t = {clock.get_value():4.0f} Myr", font_size=20, color=P.FG, weight="BOLD")
            .next_to(legend, DOWN, buff=0.35).align_to(legend, LEFT))
        cap3 = layout.caption(f"The crate's N-body: two primordial disks with a ≈ {a_p:.0f} AU",
                              font_size=21)
        self.play(FadeIn(cap3), FadeIn(legend), FadeIn(t_lab),
                  FadeIn(narrow), FadeIn(broad), run_time=1.0)
        timing.hold_to_read(self, cap3, settle=0.2)
        cap4 = layout.caption("Neptune scatters the narrow disk; the broad disk glides "
                              "along the secular paths", font_size=21)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.5)
        self.play(clock.animate.set_value(t_end), run_time=6.0, rate_func=rate_functions.linear)

        lost = paper.result_readout(f"lost in {t_end:.0f} Myr",
                                    f"narrow {100 * d['narrow_lost']:.0f}%   "
                                    f"broad {100 * d['broad_lost']:.0f}%", color=P.ORANGE).scale(0.8)
        lost.next_to(t_lab, DOWN, buff=0.35).align_to(legend, LEFT)
        dim = layout.label("faded: scattered away from a ≈ 345 AU", font_size=15, color=P.FG)
        dim.next_to(lost, DOWN, buff=0.2).align_to(legend, LEFT)
        self.play(FadeIn(lost), FadeIn(dim))
        timing.hold_to_read(self, cap4, lost, settle=0.6)

        # ---- 3. the paper's 4 Gyr outcome ------------------------------
        self.play(FadeOut(VGroup(lost, dim, t_lab, legend)), FadeOut(cap4))
        head = layout.label("paper, after 4 Gyr", font_size=18, color=P.FG, weight="BOLD")
        rows_ = VGroup(
            layout.label("aligned survivors", font_size=16, color=P.GREEN),
            layout.label("narrow disk: 0.14%", font_size=16, color=P.FG),
            layout.label("broad disk: 3.70%", font_size=16, color=P.FG),
            layout.label("median q, detectable", font_size=16, color=P.GREEN),
            layout.label("anti-aligned: 44 AU", font_size=16, color=P.FG),
            layout.label("aligned: 100 AU", font_size=16, color=P.FG),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        rows_[3].shift(DOWN * 0.15)
        rows_[4:].shift(DOWN * 0.15)
        panel = VGroup(head, rows_).arrange(DOWN, buff=0.25, aligned_edge=LEFT)
        panel.next_to(ax, RIGHT, buff=0.45).align_to(ax, UP)
        cap5 = layout.caption("Only a broad disk fills the high-q aligned paths: two perihelion peaks",
                              font_size=21)
        self.play(FadeIn(panel), FadeIn(cap5))
        timing.hold_to_read(self, cap5, panel, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Aligned survivors need high q, so the early disk reached far past 36 AU.")
