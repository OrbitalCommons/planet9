"""Anderson & Kaib (2021) -- a distant planet in the detached belt's inclinations.

The paper's N-body formation models find that, with eight planets, distant
detached orbits are made only by Kozai cycles inside Neptune's resonances,
which never leave them below i ~ 20 degrees; a ninth planet adds a
low-inclination detached group. The reproduction crate
(p9-2021-detached-inclinations) does not run those formation models; it
computes the secular tug of war over the orbit planes -- the giant planets
against an inclined Planet Nine -- and the forced tilt of the belt's midplane
that results (anim.json -> papers -> p9-2021-detached-inclinations). The
first beat shows the paper's published finding, labelled as such.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Axes,
    Create,
    DashedLine,
    FadeIn,
    FadeOut,
    Line,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2021-detached-inclinations"


def strip(x0, x1, y, h, color, opacity):
    r = Rectangle(width=x1 - x0, height=h, stroke_width=0)
    return r.set_fill(color, opacity=opacity).move_to([(x0 + x1) / 2, y, 0])


class DetachedInclinations2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        cv = d["curve"]
        a = np.array([c["a_au"] for c in cv])
        gp = np.array([c["giant_period_myr"] for c in cv])
        p9p = np.array([c["p9_period_myr"] for c in cv])
        forced = np.array([c["forced_deg"] for c in cv])
        floor = d["floor_deg"]
        i9 = d["p9"]["i_deg"]
        self.add(paper.scene_header(CRATE))

        # 1. the paper's finding: which inclinations each model can make
        x0, x1 = -5.0, 5.0
        imax = 60.0

        def xi(i):
            return x0 + (x1 - x0) * i / imax

        scale = VGroup(Line([x0, 0, 0], [x1, 0, 0], color=P.MUTED, stroke_width=2))
        for i in range(0, 61, 10):
            scale.add(Line([xi(i), -0.08, 0], [xi(i), 0.08, 0], color=P.MUTED, stroke_width=2))
            scale.add(layout.label(f"{i}°", font_size=16, color=P.MUTED)
                      .move_to([xi(i), -0.35, 0]))
        scale.add(layout.label("inclination of distant detached orbits", font_size=18)
                  .move_to([0, -0.8, 0]))
        eight_t = layout.label("8 planets: only Kozai cycles inside Neptune's resonances",
                               font_size=19).move_to([0, 2.3, 0])
        eight = strip(xi(floor), x1, 1.55, 0.7, P.FG, 0.25)
        gap8 = strip(x0, xi(floor), 1.55, 0.7, P.RED, 0.18)
        gap8_l = layout.label(f"none below ~{floor:.0f}°", font_size=18, color=P.RED)
        gap8_l.move_to(gap8)
        nine_t = layout.label("9 planets: a low-inclination detached group appears",
                              font_size=19, color=P.BLUE).move_to([0, -1.55, 0])
        nine = strip(x0, xi(floor), -2.25, 0.7, P.BLUE, 0.3)
        nine2 = strip(xi(floor), x1, -2.25, 0.7, P.FG, 0.25)
        tag = layout.label("paper's N-body formation models", font_size=15, color=P.MUTED)
        tag.to_edge(RIGHT, buff=0.5).shift(UP * 2.95)
        self.play(FadeIn(scale), FadeIn(tag))
        self.play(FadeIn(eight_t), FadeIn(eight), FadeIn(gap8), FadeIn(gap8_l))
        cap = layout.caption("Kozai cycles lift perihelia only by tilting orbits steeply",
                             font_size=22)
        self.play(FadeIn(cap))
        timing.hold_to_read(self, cap, eight_t, settle=0.4)
        self.play(FadeIn(nine_t), FadeIn(nine), FadeIn(nine2))
        cap2 = layout.caption("So cold, distant, detached orbits would point to a ninth planet",
                              font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        timing.hold_to_read(self, cap2, nine_t, settle=0.8)
        self.play(FadeOut(VGroup(scale, tag, eight_t, eight, gap8, gap8_l, nine_t, nine, nine2,
                                 cap2)))

        # 2. who sets the orbit planes: the giant planets or Planet Nine
        ax = Axes(x_range=[150, 700, 50], y_range=[2, 3.7, 0.5], x_length=8.4, y_length=4.3,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.move_to([-1.2, 0.3, 0])
        ax.get_x_axis().add_numbers(font_size=16)
        yt = VGroup(*[layout.label(str(v), font_size=15).next_to(ax.c2p(150, np.log10(v)), LEFT,
                                                                  buff=0.1)
                      for v in (100, 300, 1000, 3000)])
        xl = layout.label("semi-major axis a (AU)", font_size=18).next_to(ax, DOWN, buff=0.45)
        yl = layout.label("time for the orbit plane to turn once (Myr)", font_size=15)
        yl.rotate(np.pi / 2).next_to(yt, LEFT, buff=0.12)
        g_curve = widgets.curve(ax, a, np.log10(gp), color=P.FG)
        p_curve = widgets.curve(ax, a, np.log10(p9p), color=P.BLUE)
        g_l = layout.label("giant planets", font_size=17).next_to(
            ax.c2p(a[-1], np.log10(gp[-1])), RIGHT, buff=0.12)
        p_l = layout.label("Planet Nine", font_size=17, color=P.BLUE).next_to(
            ax.c2p(a[-1], np.log10(p9p[-1])), RIGHT, buff=0.12)
        self.play(Create(ax), FadeIn(yt), FadeIn(xl), FadeIn(yl))
        cap3 = layout.caption(f"Two pulls on each plane (orbit with q = {d['curve_q_au']:.0f} AU; "
                              f"Planet Nine tilted {i9:.0f}°): the faster one wins",
                              font_size=22)
        self.play(FadeIn(cap3), Create(g_curve), FadeIn(g_l), run_time=1.4)
        self.play(Create(p_curve), FadeIn(p_l), run_time=1.4)
        timing.hold_to_read(self, cap3, settle=0.4)
        ae = d["a_equal_au"]
        mark = widgets.marker_line(ax, ae, (2, 3.7), f"{ae:.0f} AU", color=P.ORANGE, side=UP)
        cap4 = layout.caption(f"Beyond ~{ae:.0f} AU Planet Nine turns the planes faster "
                              "than the giant planets do", font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        self.play(Create(mark))
        timing.hold_to_read(self, cap4, settle=0.8)
        self.play(FadeOut(VGroup(ax, yt, xl, yl, g_curve, p_curve, g_l, p_l, mark, cap4)))

        # 3. the belt's midplane bends toward Planet Nine
        sx = 10.0 / (a[-1] - 0)
        origin = np.array([-5.4, -1.2, 0])
        sun = orbits.sun(radius=0.12).move_to(origin)
        ecl = DashedLine(origin, origin + RIGHT * 11.0, color=P.ORANGE, stroke_width=2)
        ecl_l = layout.label("planets' plane", font_size=16, color=P.ORANGE).next_to(
            origin + RIGHT * 11.0, UP, buff=0.12).shift(LEFT * 0.7)
        th9 = np.deg2rad(i9)
        p9_line = DashedLine(origin, origin + 10.4 * np.array([np.cos(th9), np.sin(th9), 0]),
                             color=P.BLUE, stroke_width=2)
        p9_l = layout.label(f"Planet Nine's plane ({i9:.0f}°)", font_size=16, color=P.BLUE)
        p9_l.rotate(th9).next_to(p9_line.get_end(), LEFT, buff=0.2).shift(UP * 0.35 + RIGHT * 0.1)
        dist = VGroup()
        for av in (150, ae, 700):
            x = origin + sx * av * RIGHT
            dist.add(Line(x + DOWN * 0.08, x + UP * 0.08, color=P.MUTED, stroke_width=2))
            dist.add(layout.label(f"{av:.0f} AU", font_size=15, color=P.MUTED)
                     .next_to(x, DOWN, buff=0.15))
        reach = ValueTracker(a[0])

        def warp():
            g = VGroup()
            for ak, fk in zip(a[::4], forced[::4]):
                if ak > reach.get_value():
                    break
                th = np.deg2rad(fk)
                c = origin + sx * ak * np.array([np.cos(th), np.sin(th), 0])
                u = 0.2 * np.array([np.cos(th), np.sin(th), 0])
                g.add(Line(c - u, c + u, color=P.GREEN, stroke_width=5))
            return g

        belt = always_redraw(warp)
        readout = always_redraw(lambda: layout.label(
            f"a = {reach.get_value():.0f} AU:  midplane tilt "
            f"{np.interp(reach.get_value(), a, forced):.1f}°", font_size=20,
            color=P.GREEN).move_to([-2.2, 2.5, 0]))
        cap5 = layout.caption("Edge-on: each green tick is the belt's local midplane, "
                              "as the secular model forces it", font_size=22)
        self.play(FadeIn(sun), Create(ecl), FadeIn(ecl_l), Create(p9_line), FadeIn(p9_l),
                  FadeIn(dist), FadeIn(cap5))
        self.add(belt, readout)
        self.play(reach.animate.set_value(a[-1]), run_time=6.0, rate_func=linear)
        belt.clear_updaters()
        readout.clear_updaters()
        timing.hold_to_read(self, cap5, settle=0.4)
        cap6 = layout.caption("The far belt leans toward Planet Nine's plane: its inclinations "
                              "carry the planet's mark", font_size=22)
        self.play(FadeOut(cap5), FadeIn(cap6), run_time=0.6)
        timing.hold_to_read(self, cap6, settle=0.8)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, f"Beyond ~{ae:.0f} AU Planet Nine, not Neptune, sets the orbit planes.")
