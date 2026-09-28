"""Gomes, Deienno & Morbidelli (2016) -- the planetary plane tilts, the Sun stays.

In this paper the Sun's spin axis is the fixed reference. The giant planets'
invariant plane precesses about the total angular momentum of the Solar System,
Planet Nine included, and so walks away from the solar equator. The paper maps
the tilt reached after 4.5 Gyr onto the planet's mass, semi-major axis,
eccentricity and inclination, and reads off which orbits give the observed 5.9
degrees.

Everything drawn comes from anim.json -> papers -> p9-2016-obliquity-gomes: the
track of the planetary-plane pole, the tilt against eccentricity, and the
eccentricity at which the tilt reaches the target; the paper's eccentricities
appear as labelled comparisons.
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
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2016-obliquity-gomes"

DEG = 0.36  # scene units per degree on the pole plot
ORIGIN = np.array([-2.6, 0.2, 0.0])


def pole(x_deg, y_deg):
    return ORIGIN + DEG * np.array([x_deg, y_deg, 0.0])


def pole_track(case, color):
    """The cone the planetary-plane pole walks around, the part walked in
    4.5 Gyr, and where the pole is today -- all relative to the Sun's spin axis,
    which sits where the pole started."""
    r = case["forced_inclination_deg"]
    pts = [pole(x - r, y) for x, y in case["pole_track_deg"]]
    cone = DashedVMobject(
        Circle(radius=r * DEG, color=color, stroke_width=1.6).move_to(pole(-r, 0)),
        num_dashes=48)
    cone.set_stroke(opacity=0.55)
    walked = widgets.curve(_Identity(), [p[0] for p in pts], [p[1] for p in pts],
                           color=color, stroke_width=4.5)
    axis = Dot(pole(-r, 0), radius=0.07, color=color)
    return cone, walked, axis, pts


def walker(pts, t_gyr, tilt_deg, color):
    """A pole that walks the track as ``clock`` runs over the age of the Solar
    System, the chord from the Sun's spin axis to it (the tilt, to scale) and a
    live readout of that tilt."""
    clock = ValueTracker(0.0)

    def k():
        return min(int(round(clock.get_value())), len(pts) - 1)

    tip = always_redraw(lambda: Dot(pts[k()], radius=0.1, color=color).set_z_index(5))
    chord = always_redraw(lambda: Line(pole(0, 0), pts[k()], color=P.FG, stroke_width=2.4))

    def readout():
        lab = layout.label(f"{t_gyr[k()]:.1f} Gyr:  tilt {tilt_deg[k()]:.1f}°", font_size=17,
                           color=color, weight="BOLD")
        return lab.move_to([0.0, -2.35, 0], aligned_edge=LEFT)

    return clock, tip, chord, always_redraw(readout)


class _Identity:
    """Stands in for an Axes when the points are already scene coordinates."""

    @staticmethod
    def c2p(x, y):
        return np.array([x, y, 0.0])


class ObliquityGomes2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        nominal, solved = d["nominal"], d["solved"]
        target = d["target_tilt_deg"]

        self.add(paper.scene_header(CRATE))

        # 1. the pole of the planets' plane walks away from the Sun's spin axis
        spin = Dot(pole(0, 0), radius=0.1, color=P.SUN).set_z_index(5)
        spin_lab = layout.label("Sun's spin axis\n(stays put)", font_size=16, color=P.SUN,
                                line_spacing=0.9)
        spin_lab.next_to(spin, RIGHT, buff=0.15)
        ring = DashedVMobject(Circle(radius=target * DEG, color=P.GREEN, stroke_width=2.2)
                              .move_to(pole(0, 0)), num_dashes=60)
        ring_lab = layout.label(f"observed tilt  {target:.1f}°", font_size=16, color=P.GREEN)
        ring_lab.next_to(pole(0, target), UP, buff=0.08)

        cap = layout.caption(
            "Looking down the Sun's spin axis: where does the planets' plane point?",
            font_size=22)
        self.play(FadeIn(spin), FadeIn(spin_lab), Create(ring), FadeIn(ring_lab), FadeIn(cap),
                  run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.6)

        def notes(case, color, title):
            g = VGroup(
                layout.label(title, font_size=18, color=color, weight="BOLD"),
                layout.label(
                    f"{case['mass_earth']:.0f} M⊕, {case['a_au']:.0f} AU, "
                    f"e = {case['e']:.2f}, inclined {case['i9_deg']:.0f}°",
                    font_size=16, color=color),
                layout.label(f"circle radius  {case['forced_inclination_deg']:.1f}°",
                             font_size=16),
                layout.label(f"one lap takes  {case['precession_period_gyr']:.0f} Gyr",
                             font_size=16),
                layout.label(f"tilt after {d['age_gyr']:.1f} Gyr  {case['tilt_today_deg']:.1f}°",
                             font_size=16, weight="BOLD"),
            ).arrange(DOWN, buff=0.13, aligned_edge=LEFT)
            return g

        cone, walked, axis, pts = pole_track(nominal, P.BLUE)
        clock, tip, chord, live = walker(pts, nominal["t_gyr"], nominal["tilt_deg"], P.BLUE)
        n1 = notes(nominal, P.BLUE, "Batygin & Brown's planet")
        n1.move_to([4.2, 1.55, 0])
        pivot_txt = layout.label("pivot: total angular momentum,\nPlanet Nine included",
                                 font_size=15, color=P.FG, line_spacing=0.9)
        pivot_txt.move_to([axis.get_x(), -2.45, 0])
        pivot_lab = VGroup(pivot_txt, Arrow(pivot_txt.get_top(), axis.get_center(), buff=0.1,
                                            color=P.MUTED, stroke_width=2,
                                            max_tip_length_to_length_ratio=0.08))
        cap2 = layout.caption(
            "Planet Nine makes the planets' pole circle a pivot, slowly leaving the Sun behind",
            font_size=22)
        self.play(Create(cone), FadeIn(axis), FadeIn(pivot_lab), FadeIn(n1[:2]), FadeOut(cap),
                  FadeIn(cap2), run_time=1.2)
        self.add(chord, tip, live)
        self.play(Create(walked), clock.animate.set_value(len(pts) - 1), rate_func=lambda s: s,
                  run_time=4.0)
        self.play(FadeIn(n1[2:]), run_time=0.6)
        timing.hold_to_read(self, cap2, n1, settle=0.8)

        cone2, walked2, axis2, pts2 = pole_track(solved, P.ORANGE)
        clock2, tip2, chord2, live2 = walker(pts2, solved["t_gyr"], solved["tilt_deg"], P.ORANGE)
        n2 = notes(solved, P.ORANGE, "The orbit that reaches the tilt")
        n2.next_to(n1, DOWN, buff=0.4, aligned_edge=LEFT)
        cap3 = layout.caption(
            "A more eccentric orbit turns the plane faster and reaches the observed tilt",
            font_size=22)
        self.play(Create(cone2), FadeIn(axis2), FadeIn(n2[:2]), FadeOut(cap2), FadeIn(cap3),
                  FadeOut(chord), FadeOut(live), FadeOut(pivot_lab), run_time=1.2)
        self.add(chord2, tip2, live2)
        self.play(Create(walked2), clock2.animate.set_value(len(pts2) - 1),
                  rate_func=lambda s: s, run_time=4.0)
        self.play(FadeIn(n2[2:]), run_time=0.6)
        timing.hold_to_read(self, cap3, n2, settle=1.0)
        self.play(FadeOut(VGroup(spin, spin_lab, ring, ring_lab, cone, walked, axis, tip,
                                 cone2, walked2, axis2, tip2, chord2, live2, n1, n2, cap3)))

        # 2. which eccentricity gives the tilt
        top = 9.0
        ax, labels = widgets.labeled_axes(
            [0, 0.9, 0.1], [0, top, 1], x_label="eccentricity of Planet Nine",
            y_label=f"tilt after {d['age_gyr']:.1f} Gyr  (deg)", y_rotate=True, numbers=True,
            x_length=9.6, y_length=4.2, shift_down=-0.35, font_size=22)
        VGroup(ax, labels).shift(LEFT * 0.9)
        goal = DashedLine(ax.c2p(0, target), ax.c2p(0.9, target), color=P.GREEN,
                          stroke_width=2.5)
        goal_lab = layout.label(f"observed  {target:.1f}°", font_size=16, color=P.GREEN)
        goal_lab.next_to(ax.c2p(0, target), UP, buff=0.08).shift(RIGHT * 0.95)

        curves, marks, legend = VGroup(), VGroup(), VGroup()
        for scan, color in zip(d["eccentricity_scans"], (P.TEAL, P.ORANGE)):
            e = np.array(scan["e"])
            tilt = np.minimum(np.array(scan["tilt_deg"]), top)
            curves.add(widgets.curve(ax, e, tilt, color=color, stroke_width=3.2))
            marks.add(Dot(ax.c2p(scan["required_e"], target), radius=0.08, color=color)
                      .set_z_index(4))
            marks.add(DashedLine(ax.c2p(scan["published_required_e"], 0),
                                 ax.c2p(scan["published_required_e"], target), color=P.FG,
                                 stroke_width=1.6))
            legend.add(layout.label(
                f"{scan['a_au']:.0f} AU:  reproduced e = {scan['required_e']:.2f},  "
                f"paper {scan['published_required_e']:.2f}", font_size=16, color=color))
        legend.add(layout.label(
            f"all for {nominal['mass_earth']:.0f} M⊕ inclined {nominal['i9_deg']:.0f}°",
            font_size=16, color=P.MUTED))
        legend.arrange(DOWN, buff=0.13, aligned_edge=LEFT)
        legend.move_to(ax.c2p(0.25, 7.9))
        bb = Dot(ax.c2p(nominal["e"], nominal["tilt_today_deg"]), radius=0.08,
                 color=P.BLUE).set_z_index(4)
        bb_lab = layout.label("Batygin & Brown", font_size=16, color=P.BLUE)
        bb_lab.next_to(bb, UP, buff=0.08).shift(LEFT * 0.6)

        cap4 = layout.caption(
            "The analytic map, read backwards: the tilt fixes the eccentricity", font_size=22)
        self.play(Create(ax), FadeIn(labels), Create(goal), FadeIn(goal_lab), run_time=1.0)
        self.play(Create(curves), FadeIn(cap4), run_time=1.8)
        self.play(FadeIn(marks), FadeIn(legend), FadeIn(bb), FadeIn(bb_lab), run_time=0.8)
        timing.hold_to_read(self, cap4, legend, settle=1.2)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, "Read this way, the Sun's tilt is a measurement of Planet Nine's orbit.")
