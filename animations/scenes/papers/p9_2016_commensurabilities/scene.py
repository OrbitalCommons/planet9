"""de la Fuente Marcos, de la Fuente Marcos & Aarseth (2016) -- N-body experiments.

The six objects first linked to Planet Nine are integrated under the nominal
planet (10 Earth masses, a = 700 AU, e = 0.6) to see whether it keeps their
orbits confined. The scene shows the starting configuration, then sweeps time
forward through the computed perihelion distances: the orbits whose perihelia
are dragged down to Neptune's distance get kicked and come apart, the others
hold. The integration is a reduced-scale one (one clone per object) run with
the p9-core hybrid integrator; every curve and number is from
anim.json -> papers -> p9-2016-commensurabilities.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Axes,
    Create,
    Cross,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2016-commensurabilities"

Q_LO, Q_HI = 20.0, 90.0
NEPTUNE_AU = 30.0


def mini_axes(t_end, centre):
    """One small panel: time (Myr) against perihelion distance (AU)."""
    ax = Axes(x_range=[0, t_end, 10], y_range=[Q_LO, Q_HI, 10], x_length=3.5, y_length=1.55,
              axis_config={"color": P.MUTED, "include_tip": False, "stroke_width": 1.5})
    ax.move_to(centre)
    return ax


def verdict_text(o):
    if o["lost_myr"] is not None:
        return f"lost at {o['lost_myr']:.0f} Myr"
    if o["unstable_myr"] is not None:
        return f"a wanders {100 * o['a_excursion']:.0f}%"
    return f"holds: a within {100 * o['a_excursion']:.0f}%"


class Commensurabilities2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        planet = d["planet"]
        t_end = d["t_myr"]

        self.add(paper.scene_header(CRATE))

        # 1. the starting configuration, seen from above
        # The view is turned so the planet's long axis lies along the frame.
        scale = 1.0 / 290.0
        turn = -np.radians(planet["varpi_deg"])
        centre = np.array([-1.9, 0.2, 0.0])
        sun = orbits.sun(radius=0.08)
        p9 = orbits.ellipse_orbit(planet["a_au"] * scale, planet["e"], color=P.BLUE,
                                  varpi=np.radians(planet["varpi_deg"]) + turn,
                                  stroke_width=3.0)
        swarm = VGroup(*[
            orbits.ellipse_orbit(o["a_au"][0] * scale, o["e0"], color=P.GREEN,
                                 varpi=np.radians(o["varpi0_deg"]) + turn, stroke_width=2.0,
                                 opacity=0.9)
            for o in objs])
        view = VGroup(sun, p9, swarm).shift(centre)
        names = VGroup(
            layout.label("the six objects", font_size=20, color=P.GREEN, weight="BOLD"),
            *[layout.label(f"{o['name']}   a = {o['a_au'][0]:.0f} AU", font_size=17, color=P.FG)
              for o in objs],
            layout.label("nominal Planet Nine", font_size=20, color=P.BLUE, weight="BOLD"),
            layout.label(
                f"{planet['mass_earth']:.0f} M⊕   a = {planet['a_au']:.0f} AU   "
                f"e = {planet['e']:.1f}   i = {planet['i_deg']:.0f}°",
                font_size=17, color=P.FG),
        ).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
        names[7].shift(DOWN * 0.25)
        names[8].shift(DOWN * 0.25)
        names.move_to([4.9, 0.3, 0])
        cap = layout.caption("Six real orbits and the planet proposed to hold them",
                             font_size=22)
        self.play(FadeIn(sun),
                  LaggedStart(*[Create(o) for o in swarm], lag_ratio=0.15),
                  FadeIn(names[:7]), run_time=1.8)
        self.play(Create(p9), FadeIn(names[7:]), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, names[8], settle=1.0)
        self.play(FadeOut(VGroup(view, names, cap)))

        # 2. sweep time forward: perihelion distance of each object
        centres = [np.array([x, y, 0]) for y in (1.3, -1.2) for x in (-4.3, 0.0, 4.3)]
        panels = VGroup()
        furniture = VGroup()
        for k, (o, c) in enumerate(zip(objs, centres)):
            ax = mini_axes(t_end, c)
            nline = DashedLine(ax.c2p(0, NEPTUNE_AU), ax.c2p(t_end, NEPTUNE_AU),
                               color=P.ORANGE, stroke_width=1.6, dash_length=0.08)
            title = layout.label(o["name"], font_size=18, color=P.FG, weight="BOLD")
            title.next_to(ax, UP, buff=0.12).align_to(ax, LEFT).shift(RIGHT * 0.1)
            ticks = VGroup()
            for q in (30, 60, 90):
                ticks.add(layout.label(f"{q}", font_size=12, color=P.MUTED)
                          .next_to(ax.c2p(0, q), LEFT, buff=0.08))
            for t in (0, t_end / 2, t_end):
                ticks.add(layout.label(f"{t:.0f}", font_size=12, color=P.MUTED)
                          .next_to(ax.c2p(t, Q_LO), DOWN, buff=0.08))
            panels.add(ax)
            furniture.add(VGroup(nline, title, ticks))
        nep_key = layout.label("Neptune's orbit (30 AU)", font_size=13, color=P.ORANGE)
        nep_key.next_to(panels[2].c2p(t_end, NEPTUNE_AU), UP, buff=0.05).align_to(
            panels[2], RIGHT)
        y_name = layout.label("perihelion distance (AU)", font_size=15, color=P.FG)
        y_name.rotate(np.pi / 2).next_to(VGroup(panels[0], panels[3]), LEFT, buff=0.42)
        x_name = layout.label("time from now (Myr)", font_size=15, color=P.FG)
        x_name.next_to(panels[4], DOWN, buff=0.3)
        cap2 = layout.caption(
            "Dashed line: Neptune's distance. A perihelion that sinks to it gets kicked",
            font_size=22)
        self.play(FadeIn(panels), FadeIn(furniture), FadeIn(y_name), FadeIn(x_name),
                  FadeIn(nep_key), FadeIn(cap2), run_time=1.2)

        clock = ValueTracker(0.0)

        def colour_at(o, t):
            u = o["unstable_myr"]
            return P.RED if (u is not None and t >= u) else P.GREEN

        def trace(o, ax):
            def build():
                t = clock.get_value()
                pts = [ax.c2p(tt, min(max(q, Q_LO), Q_HI))
                       for tt, q in zip(o["t_myr"], o["q_au"]) if tt <= t]
                m = VMobject(stroke_width=2.4, color=colour_at(o, t))
                if len(pts) >= 2:
                    m.set_points_as_corners(pts)
                else:
                    m.set_points_as_corners([ax.c2p(0, o["q_au"][0])] * 2)
                return m
            return always_redraw(build)

        traces = VGroup(*[trace(o, ax) for o, ax in zip(objs, panels)])
        time_lab = always_redraw(lambda: layout.label(
            f"t = {clock.get_value():4.1f} Myr", font_size=17, color=P.FG)
            .move_to([5.4, 2.95, 0]))
        timing.hold_to_read(self, cap2, settle=0.4)
        cap2b = layout.caption(
            f"Each orbit integrated {t_end:.0f} Myr forward under the nominal planet",
            font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap2b))
        self.add(traces, time_lab)
        self.play(clock.animate.set_value(t_end), run_time=7.0, rate_func=rate_functions.linear)
        for tr in traces:
            tr.clear_updaters()
        time_lab.clear_updaters()

        marks = VGroup()
        verdicts = VGroup()
        for o, ax in zip(objs, panels):
            col = P.RED if o["unstable_myr"] is not None else P.GREEN
            v = layout.label(verdict_text(o), font_size=14, color=col, weight="BOLD")
            v.next_to(ax, UP, buff=0.14).align_to(ax, RIGHT)
            verdicts.add(v)
            if o["lost_myr"] is not None:
                q_last = min(max(o["q_au"][-1], Q_LO), Q_HI)
                marks.add(Cross(Dot(ax.c2p(o["t_myr"][-1], q_last), radius=0.08),
                                color=P.RED, stroke_width=3))
        n_bad = d["n_unstable"]
        cap3 = layout.caption(
            f"{n_bad} of {d['n_objects']} perihelia sink to Neptune; those orbits are kicked apart",
            font_size=22)
        self.play(FadeIn(verdicts), FadeIn(marks), FadeOut(cap2b), FadeIn(cap3))
        timing.hold_to_read(self, cap3, verdicts, settle=1.2)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, "The nominal planet holds three of its six orbits and scatters the rest.")
