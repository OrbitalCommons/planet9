"""Lu & Laughlin (2022) -- tilting Uranus with a migrating Planet Nine.

As Planet Nine migrates outward, the precession of Uranus' orbit slows until it
matches the precession of Uranus' spin axis. Caught in that secular spin-orbit
resonance, the spin axis is dragged over. The paper tests this in N-body runs
and finds it works only if Uranus once precessed about a hundred times faster
than it does today.

Everything drawn comes from anim.json -> papers -> p9-2022-uranus-tilt: the
resonant equilibrium, the integrated sweep through resonance, the peak tilt
against precession rate, and Uranus' present precession constant, all from the
reproduction crate's single-resonance spin model; the paper's N-body numbers
appear as labelled comparisons.
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
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2022-uranus-tilt"


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


TILT_TICKS = {v: f"{v}°" for v in range(0, 121, 30)}


class UranusTilt2022(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub, show, band = d["published"], d["showcase"], d["band"]
        uranus = d["uranus_obliquity_deg"]
        alpha_show = d["showcase_alpha_arcsec_yr"]

        self.add(paper.scene_header(CRATE))

        # 1. the resonant equilibrium the spin axis follows
        state = d["cassini_state_2"]
        ax = Plot(
            (0.1, 25), (0, 120),
            {0.1: "0.1", 0.3: "0.3", 1: "1", 3: "3", 10: "10"}, TILT_TICKS,
            "spin precession rate ÷ orbit precession rate", "tilt of the spin axis",
            size=(9.6, 4.1), centre=(-0.2, 0.45), xlog=True)
        branch = widgets.curve(ax, state["ratio"], state["obliquity_deg"], color=P.TEAL,
                               stroke_width=3.5)
        today = DashedLine(ax.c2p(0.1, uranus), ax.c2p(25, uranus), color=P.GREEN,
                           stroke_width=2.5)
        today_lab = layout.label(f"Uranus today: {uranus:.0f}°", font_size=16, color=P.GREEN)
        today_lab.next_to(ax.c2p(0.1, uranus), UP, buff=0.08, aligned_edge=LEFT).shift(RIGHT * 0.1)
        match = DashedLine(ax.c2p(1, 0), ax.c2p(1, 90), color=P.MUTED, stroke_width=2)
        match_lab = layout.label("rates match", font_size=16, color=P.MUTED)
        match_lab.next_to(ax.c2p(1, 8), RIGHT, buff=0.12)
        branch_lab = layout.label("where a spin axis caught\nin the resonance sits",
                                  font_size=16, color=P.TEAL, line_spacing=0.9)
        branch_lab.move_to(ax.c2p(7, 45))

        cap = layout.caption(
            "As Planet Nine moves out, Uranus' orbit precesses more slowly: the ratio climbs",
            font_size=22)
        self.play(FadeIn(ax), Create(today), FadeIn(today_lab), run_time=1.0)
        self.play(Create(match), FadeIn(match_lab), Create(branch), FadeIn(branch_lab),
                  FadeIn(cap), run_time=1.8)
        timing.hold_to_read(self, cap, branch_lab, settle=0.8)
        self.play(FadeOut(VGroup(ax, branch, today, today_lab, match, match_lab, branch_lab,
                                 cap)))

        # 2. the integrated sweep
        t_end = show["t_myr"][-1]
        ax2 = Plot(
            (0, t_end), (0, 120),
            {v: f"{v}" for v in range(0, int(t_end) + 1, 50)}, TILT_TICKS,
            "time  (million years)", "tilt of Uranus' spin axis",
            size=(8.6, 4.1), centre=(-1.7, 0.45))
        sweep = widgets.curve(ax2, show["t_myr"], show["obliquity_deg"], color=P.ORANGE,
                              stroke_width=3.5)
        today2 = DashedLine(ax2.c2p(0, uranus), ax2.c2p(t_end, uranus), color=P.GREEN,
                            stroke_width=2.5)
        today2_lab = layout.label(f"Uranus today: {uranus:.0f}°", font_size=16, color=P.GREEN)
        today2_lab.next_to(ax2.c2p(0, uranus), UP, buff=0.08, aligned_edge=LEFT)
        today2_lab.shift(RIGHT * 0.1)
        notes = VGroup(
            layout.label(f"spin precession constant\nα = {alpha_show:.1f}″ per year",
                         font_size=15, color=P.ORANGE, weight="BOLD", line_spacing=0.9),
            layout.label(f"reproduced peak:  {show['peak_deg']:.0f}°", font_size=15,
                         color=P.ORANGE),
            layout.label(f"paper, N-body:  {pub['showcase_peak_deg']:.0f}°", font_size=15),
            layout.label(
                f"paper, all stable runs:\n{pub['peak_range_deg'][0]:.0f}° to "
                f"{pub['peak_range_deg'][1]:.0f}°", font_size=15, line_spacing=0.9),
            layout.label("one resonance alone\ncannot pass 90°", font_size=16, color=P.MUTED,
                         line_spacing=0.9),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        notes.move_to([5.0, 1.0, 0])

        # Uranus itself, its spin axis tipping as the sweep is integrated
        clock = ValueTracker(0.0)

        def tilt_now():
            return float(np.interp(clock.get_value(), show["t_myr"], show["obliquity_deg"]))

        hub = np.array([5.0, -1.85, 0.0])
        globe = Circle(radius=0.38, color=P.GREEN, stroke_width=2).set_fill(P.GREEN, opacity=0.25)
        globe.move_to(hub)
        orbit_plane = Line(hub + LEFT * 1.0, hub + RIGHT * 1.0, color=P.MUTED, stroke_width=1.6)
        axis = always_redraw(lambda: Line(hub + DOWN * 0.7, hub + UP * 0.7, color=P.ORANGE,
                                          stroke_width=3.5)
                             .rotate(-np.radians(tilt_now()), about_point=hub))
        tilt_lab = always_redraw(lambda: layout.label(
            f"spin axis tilted {tilt_now():.0f}°", font_size=16, color=P.ORANGE)
            .next_to(hub + DOWN * 0.72, DOWN, buff=0.05))
        tracer = always_redraw(lambda: Dot(ax2.c2p(clock.get_value(), tilt_now()), radius=0.07,
                                           color=P.FG).set_z_index(5))
        cap2 = layout.caption(
            "Reproduced: the slow sweep through resonance drags the spin axis onto its side",
            font_size=22)
        self.play(FadeIn(ax2), Create(today2), FadeIn(today2_lab), FadeIn(globe),
                  FadeIn(orbit_plane), FadeIn(axis), FadeIn(tilt_lab), run_time=1.0)
        self.add(tracer)
        self.play(Create(sweep), clock.animate.set_value(t_end), FadeIn(notes[0]), FadeIn(cap2),
                  rate_func=lambda s: s, run_time=5.0)
        self.play(FadeIn(notes[1:]), run_time=0.8)
        timing.hold_to_read(self, cap2, notes, settle=1.0)
        self.play(FadeOut(VGroup(ax2, sweep, today2, today2_lab, notes, cap2, globe, orbit_plane,
                                 axis, tilt_lab, tracer)))

        # 3. the catch: how fast Uranus must precess
        alpha = np.array(band["alpha_arcsec_yr"])
        peak = np.array(band["peak_deg"])
        a_now = d["alpha_today_arcsec_yr"]
        ax3 = Plot(
            (0.01, 7), (0, 120),
            {0.01: "0.01", 0.1: "0.1", 1: "1", 5: "5"}, TILT_TICKS,
            "spin precession constant α  (arcseconds per year)", "peak tilt reached",
            size=(9.6, 4.0), centre=(-0.2, 0.6), xlog=True)
        reach = widgets.curve(ax3, alpha, peak, color=P.ORANGE, stroke_width=3.5)
        today3 = DashedLine(ax3.c2p(0.01, uranus), ax3.c2p(7, uranus), color=P.GREEN,
                            stroke_width=2.5)
        today3_lab = layout.label(f"Uranus today: {uranus:.0f}°", font_size=16, color=P.GREEN)
        today3_lab.next_to(ax3.c2p(0.01, uranus), UP, buff=0.08, aligned_edge=LEFT)
        today3_lab.shift(RIGHT * 0.1)
        now_dot = Dot(ax3.c2p(a_now, band["peak_today_deg"]), radius=0.09,
                      color=P.GREEN).set_z_index(4)
        now_lab = layout.label(
            f"today's α: {a_now:.3f}″ per year\n(paper: {pub['alpha_today_arcsec_yr']:.3f})\n"
            f"reaches {band['peak_today_deg']:.0f}°",
            font_size=16, color=P.GREEN, line_spacing=0.9)
        now_lab.next_to(now_dot, DOWN, buff=0.3, aligned_edge=LEFT).shift(RIGHT * 0.35)
        need_dot = Dot(ax3.c2p(alpha_show, band["peak_showcase_deg"]), radius=0.09,
                       color=P.ORANGE).set_z_index(4)
        need_lab = layout.label(
            f"the paper's run: {alpha_show:.1f}″ per year,\n"
            f"{d['alpha_enhancement']:.0f} times today's",
            font_size=16, color=P.ORANGE, line_spacing=0.9)
        need_lab.next_to(need_dot, DOWN, buff=0.15, aligned_edge=RIGHT)
        cap3 = layout.caption(
            f"Same sweep (orbit precession {band['g_initial_arcsec_yr']:g}″ to "
            f"{band['g_final_arcsec_yr']:g}″ per year), different spin precession rates",
            font_size=22)
        self.play(FadeIn(ax3), Create(today3), FadeIn(today3_lab), run_time=1.0)
        self.play(Create(reach), FadeIn(cap3), run_time=1.6)
        self.play(FadeIn(now_dot), FadeIn(now_lab), FadeIn(need_dot), FadeIn(need_lab),
                  run_time=0.8)
        timing.hold_to_read(self, cap3, now_lab, need_lab, settle=0.8)
        cap4 = layout.caption(
            f"Today's rate stalls at {band['peak_today_deg']:.0f}°; one resonance never passes "
            "90°, so 98° needs more", font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        timing.hold_to_read(self, cap4, settle=1.0)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, "Planet Nine can tip Uranus, but 98° took a Uranus precessing ~100× faster.")
