"""Li, Hadden, Payne & Holman (2018) -- the secular dynamics of TNOs and Planet Nine.

Averaged over their orbits, Planet Nine and a distant TNO exchange only slow
torques, so the TNO's state is a point in the eccentricity-vector plane
(k, h) = e (cos Δϖ, sin Δϖ) that drifts along a level curve of the conserved
secular Hamiltonian. Two TNOs with the same 60 AU perihelion are followed:
at 150 AU the giant planets' precession wins and the apse circulates; at 300 AU
Planet Nine's torque wins and the apse librates about anti-alignment. The
switch-over sits where the giants' free precession period grows longer than the
age of the Solar System. Reproduced in p9-2018-secular-dynamics; the Hamiltonian
maps, traced trajectories, precession periods and a_crit shown here are the
crate's own (anim.json -> papers -> p9-2018-secular-dynamics).
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
    DashedVMobject,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2018-secular-dynamics"

AGE_GYR = 4.5              # age of the Solar System, a scale reference
PLANE_CENTRE = np.array([-3.55, 0.05, 0.0])
PLANE_SIDE = 4.7
SUN_POS = np.array([4.05, 0.25, 0.0])
AU_SCALE = 0.0044          # scene units per AU in the top view


def contour_path(axis, grid, level):
    """Marching squares on the exported Hamiltonian grid: every segment where
    H crosses ``level``, as (k, h) point pairs. grid[i][j] is at
    (k, h) = (axis[j], axis[i]); cells touching a null corner are skipped."""
    segs = []
    n = len(axis)
    for i in range(n - 1):
        for j in range(n - 1):
            c = (grid[i][j], grid[i][j + 1], grid[i + 1][j + 1], grid[i + 1][j])
            if any(v is None for v in c):
                continue
            xy = ((axis[j], axis[i]), (axis[j + 1], axis[i]),
                  (axis[j + 1], axis[i + 1]), (axis[j], axis[i + 1]))
            pts = []
            for m in range(4):
                v0, v1 = c[m] - level, c[(m + 1) % 4] - level
                if (v0 < 0) != (v1 < 0):
                    t = v0 / (v0 - v1)
                    p0, p1 = xy[m], xy[(m + 1) % 4]
                    pts.append((p0[0] + t * (p1[0] - p0[0]), p0[1] + t * (p1[1] - p0[1])))
            for m in range(0, len(pts) - 1, 2):
                segs.append((pts[m], pts[m + 1]))
    return segs


def contour_map(ax, portrait, n_levels=22):
    """Level curves of H at evenly spaced quantiles, as one faint VMobject."""
    axis = portrait["axis"]
    grid = portrait["h_grid"]
    vals = np.array([v for row in grid for v in row if v is not None])
    levels = np.quantile(vals, np.linspace(0.02, 0.98, n_levels))
    vm = VMobject(stroke_color=P.MUTED, stroke_width=1.3, stroke_opacity=0.75)
    for lev in levels:
        for p0, p1 in contour_path(axis, grid, lev):
            vm.start_new_path(ax.c2p(*p0))
            vm.add_line_to(ax.c2p(*p1))
    return vm


def plane():
    """The (k, h) eccentricity-vector plane, P9's perihelion along +k."""
    lim = 0.92
    ax = widgets.axes([-lim, lim, 0.5], [-lim, lim, 0.5], x_length=PLANE_SIDE,
                      y_length=PLANE_SIDE, shift_down=0)
    ax.move_to(PLANE_CENTRE)
    ax.set_stroke(opacity=0.5)
    rim = DashedVMobject(Circle(radius=ax.x_axis.unit_size * 0.9, color=P.MUTED,
                                stroke_width=1.5), num_dashes=60).move_to(ax.c2p(0, 0))
    rim_lab = layout.label("e = 0.9", font_size=13, color=P.MUTED)
    rim_lab.move_to(ax.c2p(0.5, 0.83))
    xl = layout.label("k = e cos Δϖ", font_size=15, color=P.FG)
    xl.next_to(ax, DOWN, buff=0.08)
    yl = layout.label("h = e sin Δϖ", font_size=15, color=P.FG).rotate(np.pi / 2)
    yl.next_to(ax, LEFT, buff=0.08)
    toward = Arrow(ax.c2p(0.45, -0.62), ax.c2p(0.88, -0.62), buff=0, color=P.BLUE,
                   stroke_width=3, max_tip_length_to_length_ratio=0.2)
    toward_lab = layout.label("toward P9's perihelion", font_size=13, color=P.BLUE)
    toward_lab.next_to(toward, DOWN, buff=0.06).align_to(toward, RIGHT)
    return ax, VGroup(ax, rim, rim_lab, xl, yl, toward, toward_lab)


def top_view(p9):
    """Sun, Neptune and Planet Nine's orbit seen from above, perihelion to +x."""
    s = AU_SCALE
    p9_orbit = orbits.ellipse_orbit(p9["a_au"] * s, p9["e"], color=P.BLUE, varpi=0.0,
                                    stroke_width=2.5).shift(SUN_POS)
    peri = SUN_POS + np.array([p9["a_au"] * (1 - p9["e"]) * s, 0, 0])
    p9_lab = layout.label(f"Planet Nine  a = {p9['a_au']:.0f} AU, e = {p9['e']:.1f}",
                          font_size=13, color=P.BLUE)
    p9_lab.move_to(SUN_POS + np.array([-1.75, 2.35, 0]))
    nep = Circle(radius=30 * s, color=P.MUTED, stroke_width=1).move_to(SUN_POS)
    sun = orbits.sun(radius=0.08).move_to(SUN_POS)
    title = layout.label("seen from above", font_size=15, color=P.MUTED)
    title.move_to(SUN_POS + np.array([-1.75, 2.75, 0]))
    return VGroup(title, p9_orbit, p9_lab, nep, sun, Dot(peri, radius=0.05, color=P.BLUE))


def tracer(ax, portrait):
    """One traced secular trajectory driven by a clock in Myr: the state point
    in the plane, its eccentricity vector, its trail, and the orbit it stands
    for in the top view. Returns (clock, period, trail, live mobjects)."""
    path = portrait["path"]
    k = np.array(path["k"])
    h = np.array(path["h"])
    t = np.array(path["t_myr"])
    a = portrait["a_au"]
    clock = ValueTracker(0.0)

    def state():
        tt = clock.get_value()
        return float(np.interp(tt, t, k)), float(np.interp(tt, t, h))

    trail = VMobject(stroke_color=P.GREEN, stroke_width=3)
    trail.set_points_as_corners([ax.c2p(k[0], h[0]), ax.c2p(k[0], h[0])])

    def grow_trail(m):
        n = int(np.searchsorted(t, clock.get_value()))
        pts = [ax.c2p(k[i], h[i]) for i in range(max(n, 1))]
        pts.append(ax.c2p(*state()))
        m.set_points_as_corners(pts)

    trail.add_updater(grow_trail)
    vec = always_redraw(lambda: Arrow(ax.c2p(0, 0), ax.c2p(*state()), buff=0,
                                      color=P.GREEN, stroke_width=3,
                                      max_tip_length_to_length_ratio=0.08))
    dot = always_redraw(lambda: Dot(ax.c2p(*state()), radius=0.07, color=P.GREEN))

    def etno():
        kk, hh = state()
        return orbits.ellipse_orbit(a * AU_SCALE, float(np.hypot(kk, hh)), color=P.GREEN,
                                    varpi=float(np.arctan2(hh, kk)),
                                    stroke_width=2.5).shift(SUN_POS)

    orbit = always_redraw(etno)
    stamp = always_redraw(lambda: layout.label(
        f"t = {clock.get_value():3.0f} Myr", font_size=16, color=P.FG).move_to(
            ax.c2p(-0.66, 0.84)))
    return clock, float(t[-1]), trail, VGroup(vec, dot, orbit, stamp)


def run(scene, clock, period, trail, run_time):
    scene.add(trail)
    scene.play(clock.animate.set_value(period), run_time=run_time, rate_func=linear)
    trail.clear_updaters()


class SecularDynamics2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        p9 = d["p9"]
        inner, outer = d["portraits"]
        self.add(paper.scene_header(CRATE))

        # 1. the plane: one point is one orbit
        ax, frame = plane()
        view = top_view(p9)
        self.play(FadeIn(frame), FadeIn(view), run_time=1.2)
        lab_a = layout.label(f"TNO at a = {inner['a_au']:.0f} AU, perihelion "
                             f"{inner['q_au']:.0f} AU", font_size=15, color=P.GREEN)
        lab_a.move_to(SUN_POS + np.array([0.3, -2.45, 0]))
        clock, T_in, trail, live = tracer(ax, inner)
        cap = layout.caption("One point is one orbit: its distance from the centre is e, "
                             "its angle is where the perihelion points", font_size=20)
        self.play(FadeIn(live), FadeIn(lab_a), FadeIn(cap), run_time=0.8)
        timing.hold_to_read(self, cap, settle=0.6)

        # 2. at 150 AU the apse circulates
        cmap = contour_map(ax, inner)
        cap2 = layout.caption("Grey curves: paths of constant secular energy that an orbit must follow",
                              font_size=20)
        self.play(Create(cmap), FadeOut(cap), FadeIn(cap2), run_time=1.6)
        timing.hold_to_read(self, cap2, settle=0.3)
        run(self, clock, T_in, trail, run_time=6.0)
        cap3 = layout.caption(f"At {inner['a_au']:.0f} AU the giants win: the apse swings all the "
                              f"way round every {T_in:.0f} Myr", font_size=20)
        self.play(FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.8)

        # 3. at 300 AU the apse librates about anti-alignment
        cmap2 = contour_map(ax, outer)
        lab_b = layout.label(f"TNO at a = {outer['a_au']:.0f} AU, perihelion "
                             f"{outer['q_au']:.0f} AU", font_size=15, color=P.GREEN)
        lab_b.move_to(lab_a)
        self.play(FadeOut(VGroup(live, trail, cmap, lab_a, cap3)), run_time=0.8)
        clock2, T_out, trail2, live2 = tracer(ax, outer)
        self.play(Create(cmap2), FadeIn(lab_b), FadeIn(live2), run_time=1.4)
        run(self, clock2, T_out, trail2, run_time=6.0)
        lo, hi = outer["dvarpi_min_deg"], outer["dvarpi_max_deg"]
        rays = VGroup(*[DashedLine(ax.c2p(0, 0), ax.c2p(0.9 * np.cos(np.radians(g)),
                                                        0.9 * np.sin(np.radians(g))),
                                   color=P.ORANGE, stroke_width=2) for g in (lo, hi)])
        half = 0.5 * (hi - lo)
        cap4 = layout.caption(f"At {outer['a_au']:.0f} AU P9 wins: the apse only rocks "
                              f"±{half:.0f}° about anti-alignment, every {T_out:.0f} Myr",
                              font_size=20)
        self.play(Create(rays), FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=1.0)
        self.play(FadeOut(VGroup(frame, view, cmap2, live2, trail2, rays, lab_b, cap4)))

        # 4. why the switch: the giants' precession clock runs out
        pre = d["precession"]
        a_arr = np.array(pre["a_au"])
        lg = np.log10(np.array(pre["period_myr"]) / 1000.0) + 2.0  # 0 = 10 Myr
        ax2, labels = widgets.labeled_axes(
            [75, 475, 50], [0, 4, 1], x_label="TNO semi-major axis a (AU)",
            y_label="giants' precession period", y_rotate=True, numbers=False,
            x_length=9.6, y_length=4.2, shift_down=-0.45, font_size=20)
        ax2.x_axis.add_numbers(range(100, 475, 50))
        labels[0].next_to(ax2, DOWN, buff=0.12)
        yticks = VGroup(*[layout.label(t, font_size=14, color=P.FG)
                          .next_to(ax2.c2p(75, v), LEFT, buff=0.12)
                          for v, t in zip(range(0, 5), ["10 Myr", "100 Myr", "1 Gyr",
                                                           "10 Gyr", "100 Gyr"])])
        labels[1].next_to(yticks, LEFT, buff=0.12)
        curve = widgets.curve(ax2, a_arr, lg, color=P.ORANGE)
        age = DashedLine(ax2.c2p(75, np.log10(AGE_GYR) + 2), ax2.c2p(475, np.log10(AGE_GYR) + 2),
                         color=P.MUTED, stroke_width=2)
        age_lab = layout.label(f"age of the Solar System, {AGE_GYR} Gyr", font_size=14,
                               color=P.MUTED).next_to(age, DOWN, buff=0.08).align_to(age, RIGHT)
        self.play(Create(ax2), FadeIn(labels), FadeIn(yticks))
        self.play(Create(curve), run_time=1.6)
        self.play(Create(age), FadeIn(age_lab))
        cap5 = layout.caption("The giant planets turn a wide orbit ever more slowly: period grows as a^3.5",
                              font_size=20)
        self.play(FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=0.6)

        a_crit = d["a_crit_au"]
        a_pub = d["published_a_threshold_au"]
        m_crit = widgets.marker_line(ax2, a_crit, (0, 4), f"reproduced a = {a_crit:.0f} AU",
                                     color=P.GREEN, side=LEFT)
        m_pub = widgets.marker_line(ax2, a_pub, (0, 3.55), f"paper a ≳ {a_pub:.0f} AU",
                                    color=P.FG, side=RIGHT)
        zone = layout.label("apse librates:\nP9 in control", font_size=16, color=P.GREEN)
        zone.move_to(ax2.c2p(410, 1.0))
        spin = layout.label("apse circulates:\ngiants in control", font_size=16, color=P.MUTED)
        spin.move_to(ax2.c2p(140, 3.3))
        cap6 = layout.caption("Near where one giant-planet turn outlasts the Solar System, "
                              "P9 takes over", font_size=20)
        self.play(Create(m_crit), Create(m_pub), FadeIn(zone), FadeIn(spin), FadeOut(cap5),
                  FadeIn(cap6), run_time=1.4)
        timing.hold_to_read(self, cap6, settle=1.2)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "Secular torques alone confine the apses once the giants' precession is too slow.")
