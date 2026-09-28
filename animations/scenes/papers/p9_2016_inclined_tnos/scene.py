"""Batygin & Brown (2016) -- generation of highly inclined TNOs by Planet Nine.

The same Planet Nine that clusters the distant orbits also tips scattered-disk
objects out of the plane of the planets, some past 90 degrees onto retrograde
orbits; Neptune then pulls them inward, making objects like Drac and Niku. The
scene follows one test orbit as the planet tilts it, then shows the whole
reduced-scale run: the inclination distribution before and after, and where the
tilted orbits sit against the real high-inclination objects. Reproduced in
p9-2016-inclined-tnos; every series is the crate's own run
(anim.json -> papers -> p9-2016-inclined-tnos).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arrow,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2016-inclined-tnos"

DIAL_AT = np.array([4.1, 0.55, 0.0])
DIAL_W = 3.0


def tilted_orbit(i_deg):
    """An orbit seen edge-on, tipped by ``i_deg`` from the planets' plane, with
    its spin axis (orbital angular momentum): up for prograde, below the
    plane once the orbit runs backwards."""
    orbit = Line(DIAL_AT + LEFT * DIAL_W / 2, DIAL_AT + RIGHT * DIAL_W / 2,
                 color=P.ORANGE, stroke_width=3.2)
    spin = Arrow(DIAL_AT, DIAL_AT + UP * 1.35, buff=0, color=P.ORANGE, stroke_width=5,
                 max_tip_length_to_length_ratio=0.2)
    g = VGroup(orbit, spin)
    g.rotate(np.radians(i_deg), about_point=DIAL_AT)
    return g


def log_axes(x_lo, x_hi, y_hi, x_length, y_length):
    """Axes with log10(a) along x; returns (axes, tick labels)."""
    ax = widgets.axes([np.log10(x_lo), np.log10(x_hi), 1], [0, y_hi, 30],
                      x_length=x_length, y_length=y_length, shift_down=0)
    ticks = VGroup()
    for a in (30, 100, 300):
        ticks.add(layout.label(f"{a}", font_size=13, color=P.MUTED)
                  .next_to(ax.c2p(np.log10(a), 0), DOWN, buff=0.1))
        ticks.add(Line(ax.c2p(np.log10(a), 0), ax.c2p(np.log10(a), 0) + UP * 0.08,
                       color=P.MUTED, stroke_width=1.5))
    for i in range(30, int(y_hi) + 1, 30):
        ticks.add(layout.label(f"{i}°", font_size=13, color=P.MUTED)
                  .next_to(ax.c2p(np.log10(x_lo), i), LEFT, buff=0.1))
    return ax, ticks


class InclinedTnos2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        t_end = d["t_myr"]
        track = d["tracks"][0]
        high = d["high_i_deg"]

        self.add(paper.scene_header(CRATE))

        # 1. one test orbit, tipped by the planet
        ax, labels = widgets.labeled_axes(
            [0, t_end, 10], [0, 120, 30], x_label="time (Myr)",
            y_label="inclination (deg)", y_rotate=True, numbers=True,
            x_length=6.4, y_length=4.2, shift_down=-0.25)
        VGroup(ax, labels).shift(LEFT * 2.9)
        perp = DashedLine(ax.c2p(0, 90), ax.c2p(t_end, 90), color=P.MUTED, stroke_width=1.6)
        perp_lab = layout.label("90°: perpendicular to the planets", font_size=14,
                                color=P.MUTED).next_to(ax.c2p(0, 90), UP + RIGHT, buff=0.08)

        plane = Line(DIAL_AT + LEFT * (DIAL_W / 2 + 0.4), DIAL_AT + RIGHT * (DIAL_W / 2 + 0.4),
                     color=P.MUTED, stroke_width=2.2)
        plane_arrow = Arrow(DIAL_AT, DIAL_AT + UP * 1.35, buff=0, color=P.MUTED,
                            stroke_width=3, max_tip_length_to_length_ratio=0.2)
        plane_lab = VGroup(
            layout.label("edge-on view", font_size=15, color=P.FG, weight="BOLD"),
            layout.label("grey: planets' plane and spin axis", font_size=14, color=P.MUTED),
            layout.label("orange: the test orbit and its spin axis", font_size=14,
                         color=P.ORANGE),
        ).arrange(DOWN, buff=0.08, aligned_edge=LEFT)
        plane_lab.next_to(DIAL_AT + DOWN * 1.75, DOWN, buff=0)
        sun = orbits.sun(radius=0.07).move_to(DIAL_AT)

        clock = ValueTracker(0.0)
        ts = np.array(track["t_myr"])
        incl = np.array(track["i_deg"])

        def i_now():
            return float(np.interp(clock.get_value(), ts, incl))

        def grow():
            t = clock.get_value()
            pts = [ax.c2p(x, y) for x, y in zip(ts, incl) if x <= t]
            pts.append(ax.c2p(t, i_now()))
            m = VMobject(color=P.ORANGE, stroke_width=3)
            m.set_points_as_corners(pts if len(pts) > 1 else pts * 2)
            return m

        curve = always_redraw(grow)
        head = always_redraw(lambda: Dot(ax.c2p(clock.get_value(), i_now()), radius=0.07,
                                         color=P.ORANGE))
        dial = always_redraw(lambda: tilted_orbit(i_now()))
        readout = always_redraw(lambda: layout.label(
            f"t = {clock.get_value():4.1f} Myr    i = {i_now():5.1f}°", font_size=18,
            color=P.ORANGE).move_to(DIAL_AT + np.array([0, 2.05, 0])))
        cap = layout.caption(
            f"One test object starting at a = {track['a0_au']:.0f} AU, "
            f"Planet Nine at {d['planet']['a_au']:.0f} AU", font_size=22)
        self.play(Create(ax), FadeIn(labels), Create(perp), FadeIn(perp_lab),
                  FadeIn(plane), FadeIn(plane_arrow), FadeIn(plane_lab), FadeIn(sun),
                  FadeIn(cap))
        self.add(curve, head, dial, readout)
        timing.hold_to_read(self, cap, settle=0.3)
        cap2 = layout.caption(
            "Planet Nine's slow pull tips the orbit; once its spin axis dips below the plane "
            "(i > 90°) it runs backwards", font_size=20)
        self.play(FadeOut(cap), FadeIn(cap2))
        self.play(clock.animate.set_value(t_end), run_time=8.0, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap2, settle=0.8)
        for m in (curve, head, dial, readout):
            m.clear_updaters()
        self.play(FadeOut(VGroup(ax, labels, perp, perp_lab, plane, plane_arrow, plane_lab,
                                 sun, curve, head, dial, readout, cap2)))

        # 2. the whole run: inclinations before and after
        h = d["histogram"]
        centres = np.array(h["i_deg"])
        width = centres[1] - centres[0]
        edges = np.append(centres - width / 2, centres[-1] + width / 2)
        top = max(max(h["start"]), max(h["end"])) + 2
        ax2, labels2 = widgets.labeled_axes(
            [0, 120, 30], [0, top, 5], x_label="inclination (deg)",
            y_label="test objects", y_rotate=True, numbers=True,
            x_length=9.0, y_length=3.9, shift_down=-0.25)
        start = widgets.histogram(ax2, edges[:13], h["start"][:12], color=P.MUTED, opacity=0.5)
        end = widgets.histogram(ax2, edges[:13], h["end"][:12], color=P.ORANGE, opacity=0.65)
        mark = widgets.marker_line(ax2, high, (0, top), f"{high:.0f}°", color=P.FG, side=UP)
        key = VGroup(
            layout.label(f"start: all within {edges[np.nonzero(h['start'])[0][-1] + 1]:.0f}°",
                         font_size=16, color=P.MUTED),
            layout.label(f"after {t_end:.0f} Myr", font_size=16, color=P.ORANGE),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        key.move_to(ax2.c2p(92, 0.78 * top))
        cap3 = layout.caption(
            f"{d['n_particles']} scattered-disk objects, integrated {t_end:.0f} Myr", font_size=22)
        self.play(Create(ax2), FadeIn(labels2), FadeIn(cap3))
        self.play(FadeIn(start, lag_ratio=0.1), FadeIn(key[0]), run_time=1.0)
        self.play(FadeIn(end, lag_ratio=0.1), FadeIn(key[1]), start.animate.set_opacity(0.25),
                  run_time=1.4)
        self.play(Create(mark))
        tally = paper.result_readout(
            f"tipped past {high:.0f}°", f"{d['n_high']} of {d['n_particles']}",
            color=P.ORANGE).scale(0.62)
        tally.move_to(ax2.c2p(92, 0.4 * top))
        cap4 = layout.caption(
            f"{d['n_high']} climb past {high:.0f}°, {d['n_retrograde']} beyond 90° "
            f"(up to {d['max_i_deg']:.0f}°)", font_size=22)
        self.play(FadeIn(tally), FadeOut(cap3), FadeIn(cap4))
        timing.hold_to_read(self, cap4, tally, settle=1.0)
        self.play(FadeOut(VGroup(ax2, labels2, start, end, mark, key, tally, cap4)))

        # 3. where they sit against the real steep objects
        ax3, ticks3 = log_axes(20, 700, 150, x_length=9.0, y_length=4.0)
        VGroup(ax3, ticks3).shift(np.array([-0.4, 0.35, 0]) - ax3.get_center())
        xl = layout.label("semi-major axis (AU, log scale)", font_size=17, color=P.FG)
        xl.next_to(ax3, DOWN, buff=0.4)
        yl = layout.label("inclination", font_size=15, color=P.FG).rotate(np.pi / 2)
        yl.next_to(ax3, LEFT, buff=0.55)
        a_lo, a_hi = ax3.c2p(np.log10(20), high), ax3.c2p(np.log10(100), 150)
        zone = Polygon(a_lo, [a_hi[0], a_lo[1], 0], a_hi, [a_lo[0], a_hi[1], 0],
                       stroke_width=1.2, color=P.GREEN).set_fill(P.GREEN, opacity=0.07)
        zone_lab = layout.label("steep objects inside 100 AU", font_size=14, color=P.GREEN)
        zone_lab.next_to(zone, DOWN, buff=0.08)

        start_dots = VGroup(*[
            Dot(ax3.c2p(np.log10(c["a0_au"]), c["i0_deg"]), radius=0.045, color=P.MUTED)
            for c in d["cloud"]])
        peak_dots = VGroup(*[
            Dot(ax3.c2p(np.log10(c["a_peak_au"]), c["i_peak_deg"]), radius=0.055,
                color=P.ORANGE if c["i_peak_deg"] > high else P.MUTED)
            for c in d["cloud"]])
        known = VGroup()
        for o in d["known"]:
            p = ax3.c2p(np.log10(o["a_au"]), o["i_deg"])
            known.add(VGroup(Dot(p, radius=0.07, color=P.GREEN),
                             layout.label(o["name"], font_size=14, color=P.GREEN)
                             .next_to(p, RIGHT, buff=0.08)))
        cap5 = layout.caption("Each test object: where it started, then its steepest state",
                              font_size=22)
        self.play(Create(ax3), FadeIn(ticks3), FadeIn(xl), FadeIn(yl), FadeIn(start_dots),
                  FadeIn(cap5))
        self.play(*[d0.animate.move_to(d1.get_center()).set_color(d1.get_color())
                    for d0, d1 in zip(start_dots, peak_dots)], run_time=2.0)
        timing.hold_to_read(self, cap5, settle=0.4)
        cap6 = layout.caption(
            "Real steep objects (green) sit inside 100 AU; ours still need Neptune to pull them in",
            font_size=22)
        self.play(FadeIn(zone), FadeIn(zone_lab), FadeIn(known), FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=1.0)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "Planet Nine tips distant orbits steep and even retrograde.")
