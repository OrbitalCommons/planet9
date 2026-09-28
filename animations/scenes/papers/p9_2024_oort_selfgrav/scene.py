"""Batygin & Nesvorný (2024) -- self-gravity of the inner Oort cloud.

If the inner Oort cloud holds a few Earth masses, its own smooth gravity
(a flattened Miyamoto-Nagai halo) drives slow von Zeipel-Lidov-Kozai-like
cycles in which a distant orbit trades perihelion for inclination. The cycles
reach observationally interesting perihelia, but one takes longer than the age
of the Sun. Reproduced in p9-2024-oort-selfgrav: the secular trajectories and
the cycle timescales are the crate's own (anim.json -> papers ->
p9-2024-oort-selfgrav).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Circle,
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
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2024-oort-selfgrav"


def shade(ax, x0, x1, y0, y1, color, opacity):
    p0, p1 = np.array(ax.c2p(x0, y0)), np.array(ax.c2p(x1, y1))
    r = Rectangle(width=p1[0] - p0[0], height=p1[1] - p0[1], stroke_width=0)
    return r.set_fill(color, opacity=opacity).move_to((p0 + p1) / 2)


class OortSelfgrav2024(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        age = d["age_gyr"]
        cloud = d["cloud"]
        self.add(paper.scene_header(CRATE))

        # 1. the cloud to scale
        k = 2.6 / cloud["b_tilde_au"]
        c = np.array([-3.2, 0.1, 0])
        halo = VGroup(*[Circle(radius=k * r, color=P.PURPLE, stroke_width=0)
                        .set_fill(P.PURPLE, opacity=0.07).move_to(c)
                        for r in (1000, 750, 500, 300)])
        sun = orbits.sun(radius=0.03).move_to(c)
        nep = Circle(radius=k * 30.0, color=P.ORANGE, stroke_width=3).move_to(c)
        orbit = orbits.ellipse_orbit(k * d["a_au"], 1 - 300.0 / d["a_au"], color=P.GREEN,
                                     stroke_width=2.5, varpi=np.deg2rad(200)).shift(c)
        lab = VGroup(
            layout.label(f"inner Oort cloud: ~{cloud['m_ioc_earth']:.0f} Earth masses", font_size=19,
                         color=P.PURPLE),
            layout.label(f"spread over ~{cloud['b_tilde_au']:.0f} AU", font_size=19,
                         color=P.PURPLE),
            layout.label("Neptune's orbit (30 AU): the orange speck", font_size=17,
                         color=P.ORANGE),
            layout.label(f"green: a distant orbit, a = {d['a_au']:.0f} AU", font_size=17,
                         color=P.GREEN),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([3.6, 0.6, 0])
        self.play(FadeIn(halo), FadeIn(sun), Create(nep), run_time=1.2)
        self.play(FadeIn(lab[:3]))
        self.play(Create(orbit), FadeIn(lab[3]))
        cap = layout.caption("Could the cloud's own gravity reshape the orbits "
                             "inside it, as the planets do?", font_size=22)
        self.play(FadeIn(cap))
        timing.hold_to_read(self, cap, lab, settle=0.6)
        self.play(FadeOut(VGroup(halo, sun, nep, orbit, lab, cap)))

        # 2. one orbit, followed: perihelion traded for tilt
        tr = d["headline_track"]
        t = np.array(tr["t_gyr"])
        q = np.array(tr["q_au"])
        inc = np.array(tr["i_deg"])
        t_end = float(t[-1])
        step = 10 if t_end <= 60 else 20
        axq = widgets.axes([0, t_end, step], [240, 305, 20], x_length=8.8, y_length=2.1,
                           shift_down=0).move_to([-0.6, 1.55, 0])
        axi = widgets.axes([0, t_end, step], [62, 66, 1], x_length=8.8, y_length=2.1,
                           shift_down=0).move_to([-0.6, -1.4, 0])
        axq.add_coordinates()
        axi.add_coordinates()
        lq = layout.label("perihelion q (AU)", font_size=15).rotate(np.pi / 2).next_to(
            axq, LEFT, buff=0.12)
        li = layout.label("inclination (deg)", font_size=15).rotate(np.pi / 2).next_to(
            axi, LEFT, buff=0.12)
        lt = layout.label("time (Gyr)", font_size=15).next_to(axi, DOWN, buff=0.08)
        sun_q = shade(axq, 0, age, 240, 305, P.TEAL, 0.14)
        sun_i = shade(axi, 0, age, 62, 66, P.TEAL, 0.14)
        sun_l = layout.label(f"age of the Sun, {age:.1f} Gyr", font_size=15, color=P.TEAL)
        sun_l.next_to(axq.c2p(age, 305), UP, buff=0.08).shift(RIGHT * 0.6)
        prog = ValueTracker(0.0)

        def upto():
            return max(2, int(prog.get_value() * (len(t) - 1)) + 1)

        cq = always_redraw(lambda: widgets.curve(axq, t[:upto()], q[:upto()], color=P.GREEN))
        ci = always_redraw(lambda: widgets.curve(axi, t[:upto()], inc[:upto()], color=P.ORANGE))
        self.play(Create(axq), Create(axi), FadeIn(lq), FadeIn(li), FadeIn(lt))
        cap2 = layout.caption(f"Start at q = {tr['q0_au']:.0f} AU: the cloud's pull slowly "
                              "trades perihelion for tilt", font_size=22)
        self.play(FadeIn(cap2))
        self.add(cq, ci)
        self.play(prog.animate.set_value(1.0), run_time=7.0, rate_func=linear)
        cq.clear_updaters()
        ci.clear_updaters()
        timing.hold_to_read(self, cap2, settle=0.2)
        cap3 = layout.caption(f"One cycle: {tr['timescale_gyr']:.0f} Gyr. In the Sun's lifetime "
                              f"q only drifts {tr['q0_au']:.0f} → {tr['q_at_age_au']:.0f} AU",
                              font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), FadeIn(sun_q), FadeIn(sun_i), FadeIn(sun_l))
        timing.hold_to_read(self, cap3, settle=0.8)
        self.play(FadeOut(VGroup(axq, axi, lq, li, lt, cq, ci, sun_q, sun_i, sun_l, cap3)))

        # 3. the clock across the cloud
        tc = d["timescale_curves"]
        a = np.array(tc["a_au"])
        ax, labels = widgets.labeled_axes(
            [500, 3000, 500], [0, 90, 15], x_label="semi-major axis a (AU)",
            y_label="one cycle (Gyr)", y_rotate=True, numbers=True,
            x_length=8.4, y_length=4.2, shift_down=0)
        VGroup(ax, labels).move_to([-1.2, 0.25, 0])
        age_line = DashedLine(ax.c2p(500, age), ax.c2p(3000, age), color=P.TEAL, stroke_width=2.5)
        age_l = layout.label(f"age of the Sun, {age:.1f} Gyr", font_size=16, color=P.TEAL)
        age_l.next_to(ax.c2p(3000, age), RIGHT, buff=0.12)
        curves = VGroup()
        tags = VGroup()
        shown = [cv for cv in tc["curves"] if cv["m_ioc_earth"] <= cloud["m_ioc_earth"]]
        for cv, col in zip(shown, (P.PURPLE, P.FG)):
            ys = np.array(cv["timescale_gyr"])
            curves.add(widgets.curve(ax, a, ys, color=col))
            tags.add(layout.label(f"cloud of {cv['m_ioc_earth']:.0f} M⊕", font_size=16,
                                  color=col).next_to(ax.c2p(a[-1], ys[-1]), RIGHT, buff=0.12))
        cap4 = layout.caption(f"Time for one cycle (orbits with q = {tc['q_au']:.0f} AU, "
                              f"i = {tc['i_deg']:.0f}°) across the cloud", font_size=22)
        self.play(Create(ax), FadeIn(labels), Create(age_line), FadeIn(age_l), FadeIn(cap4))
        self.play(Create(curves), FadeIn(tags), run_time=1.6)
        timing.hold_to_read(self, cap4, settle=0.4)
        ts = d["timescale_gyr"]
        mark = Dot(ax.c2p(1000, ts), radius=0.09, color=P.RED).set_z_index(5)
        mark_l = layout.label(f"{ts:.0f} Gyr at 1000 AU: {d['timescale_over_age']:.1f}× the "
                              "Sun's age", font_size=17, color=P.RED).move_to(ax.c2p(1000, 40))
        lead = Line(mark_l.get_bottom() + DOWN * 0.05, mark.get_top(), color=P.RED,
                    stroke_width=1.5)
        cap5 = layout.caption("Only the innermost orbits get through a cycle; "
                              "the rest have barely begun one", font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), run_time=0.6)
        self.play(FadeIn(mark), FadeIn(mark_l), Create(lead))
        timing.hold_to_read(self, cap5, settle=0.8)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "The cloud's own gravity is real but too slow to shape the distant orbits.")
