"""Lai (2016) -- the solar obliquity from Planet Nine as a closed form.

Bailey et al. and Gomes et al. integrate the secular equations numerically.
Lai solves them: the tilt of the Sun's equator is set by two frequencies, the
rate at which Planet Nine drives the tilt and the mismatch between the
precession of the Sun's spin and that of the planets' plane.

Everything drawn comes from anim.json -> papers -> p9-2016-obliquity-lai: the
two precession periods, the closed-form tilt history set against the numerical
integration of the Bailey et al. crate for the same planet, and the tilt
against effective semi-major axis; the paper's coefficients and range appear as
labelled comparisons.
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
    Polygon,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2016-obliquity-lai"


class ObliquityLai2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        cmp, pub, ref = d["compare"], d["published"], d["reference"]
        observed = d["observed_obliquity_deg"]
        age = d["age_gyr"]

        self.add(paper.scene_header(CRATE))

        # 1. the closed form, term by term
        eq = layout.explain_equation(
            self,
            [r"\theta", "=", r"\dfrac{2\,\Omega_y}{\Omega_z}",
             r"\sin\!\left(\dfrac{\Omega_z\,t}{2}\right)"],
            [
                (2, "Planet Nine's drive, divided by how far two precessions are out of step"),
                (3, "the tilt swings as that mismatch winds up over the time t"),
            ],
            scale=1.1, where=UP * 0.6)
        self.play(FadeOut(eq))

        # 2. the formula against the numerical integration
        top = 8.0
        ax, labels = widgets.labeled_axes(
            [0, age, 0.5], [0, top, 2], x_label="time since the Solar System formed  (Gyr)",
            y_label="tilt of the Sun's equator  (deg)", y_rotate=True, numbers=True,
            x_length=8.4, y_length=4.2, shift_down=-0.2, font_size=22)
        VGroup(ax, labels).shift(LEFT * 2.0)
        goal = DashedLine(ax.c2p(0, observed), ax.c2p(age, observed), color=P.GREEN,
                          stroke_width=2.5)
        goal_lab = layout.label(f"observed today: {observed:.0f}°", font_size=16, color=P.GREEN)
        goal_lab.next_to(ax.c2p(0, observed), UP, buff=0.08).shift(RIGHT * 1.15)
        formula = widgets.curve(ax, cmp["t_gyr"], cmp["analytic_deg"], color=P.ORANGE,
                                stroke_width=3.5)
        dots = VGroup(*[Dot(ax.c2p(t, y), radius=0.05, color=P.FG)
                        for t, y in zip(cmp["numerical_t_gyr"], cmp["numerical_deg"])])

        notes = VGroup(
            layout.label(
                f"{cmp['mass_earth']:.0f} M⊕ at {cmp['a_au']:.0f} AU, e = {cmp['e']:.1f},\n"
                f"inclined {cmp['inclination_deg']:.0f}°",
                font_size=15, color=P.BLUE, weight="BOLD", line_spacing=0.9),
            layout.label(f"closed form:  {cmp['analytic_final_deg']:.1f}°", font_size=15,
                         color=P.ORANGE),
            layout.label(f"numerical integration:  {cmp['numerical_final_deg']:.1f}°",
                         font_size=15),
            layout.label(
                f"planets' plane turns once in\n{d['plane_precession_period_gyr']:.1f} Gyr"
                f"   (paper: {pub['plane_precession_period_gyr']:.1f})",
                font_size=16, line_spacing=0.9),
            layout.label(
                f"Sun's spin axis turns once in\n{d['spin_precession_period_gyr']:.1f} Gyr"
                f"   (paper: {pub['spin_precession_period_gyr']:.1f})",
                font_size=16, line_spacing=0.9),
            layout.label(
                f"periods for {ref['mass_earth']:.0f} M⊕ at {ref['a_tilde_au']:.0f} AU\n"
                f"and a {ref['spin_period_days']:.0f} day solar rotation",
                font_size=15, color=P.MUTED, line_spacing=0.9),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        notes.move_to([4.75, 0.5, 0])

        cap = layout.caption(
            "The formula, against a 4.5 Gyr numerical integration of the same planet",
            font_size=22)
        self.play(Create(ax), FadeIn(labels), Create(goal), FadeIn(goal_lab), run_time=1.0)
        self.play(Create(formula), FadeIn(notes[0]), FadeIn(notes[1]), FadeIn(cap), run_time=1.6)
        self.play(FadeIn(dots, lag_ratio=0.1), FadeIn(notes[2]), run_time=1.0)
        timing.hold_to_read(self, cap, notes[:3], settle=0.4)
        cap2 = layout.caption(
            "Its two frequencies come straight from the planets' masses and the Sun's shape",
            font_size=22)
        self.play(FadeIn(notes[3:]), FadeOut(cap), FadeIn(cap2), run_time=0.8)
        timing.hold_to_read(self, cap2, notes[3:], settle=0.6)
        self.play(FadeOut(VGroup(ax, labels, goal, goal_lab, formula, dots, notes, cap2)))

        # 3. one formula maps every planet onto a tilt
        a_tilde = np.array(d["a_tilde_au"])
        top3 = 12.0
        ax3, labels3 = widgets.labeled_axes(
            [250, 750, 50], [0, top3, 2],
            x_label="effective semi-major axis of Planet Nine,  a√(1−e²)  (AU)",
            y_label=f"tilt after {age:.1f} Gyr  (deg)", y_rotate=True, numbers=True,
            x_length=10.0, y_length=4.2, shift_down=-0.35, font_size=22)
        lo_a, hi_a = pub["a_tilde_range_au"]
        p0, p1 = ax3.c2p(lo_a, 0), ax3.c2p(hi_a, top3)
        band = Polygon(p0, [p1[0], p0[1], 0], p1, [p0[0], p1[1], 0], stroke_width=0)
        band.set_fill(P.GREEN, opacity=0.16)
        band_lab = layout.label(f"paper: {lo_a:.0f} to {hi_a:.0f} AU", font_size=16,
                                color=P.GREEN)
        band_lab.next_to(ax3.c2p(0.5 * (lo_a + hi_a), top3), UP, buff=0.08)
        goal3 = DashedLine(ax3.c2p(250, observed), ax3.c2p(750, observed), color=P.GREEN,
                           stroke_width=2.5)
        goal3_lab = layout.label(f"observed: {observed:.0f}°", font_size=16, color=P.GREEN)
        goal3_lab.next_to(ax3.c2p(610, observed), DOWN, buff=0.08)

        mass = ref["mass_earth"]
        family = [c for c in d["tilt_curves"] if c["mass_earth"] == mass]
        shades = (P.TEAL, P.ORANGE, P.RED)
        curves, marks, legend = VGroup(), VGroup(), VGroup()
        for c, color in zip(family, shades):
            curves.add(widgets.curve(ax3, a_tilde, np.minimum(c["tilt_deg"], top3), color=color,
                                     stroke_width=3.2))
            need = c["a_tilde_for_observed_au"]
            text = f"inclined {c['inclination_deg']:.0f}°:  "
            if need is None:
                text += f"never reaches {observed:.0f}°"
            else:
                text += f"{observed:.0f}° at {need:.0f} AU"
                marks.add(Dot(ax3.c2p(need, observed), radius=0.08, color=color).set_z_index(4))
            legend.add(layout.label(text, font_size=16, color=color))
        legend.add(layout.label(
            f"{mass:.0f} M⊕, solar rotation {d['spin_period_days']:.0f} days", font_size=15,
            color=P.MUTED))
        legend.arrange(DOWN, buff=0.13, aligned_edge=LEFT)
        legend.next_to(ax3.c2p(750, top3), DOWN, buff=0.1, aligned_edge=RIGHT)

        cap3 = layout.caption(
            "One formula maps every planet onto a tilt, with a peak where the precessions match",
            font_size=22)
        self.play(Create(ax3), FadeIn(labels3), Create(goal3), FadeIn(goal3_lab), run_time=1.0)
        self.play(Create(curves), FadeIn(legend), FadeIn(cap3), run_time=2.0)

        # slide a planet outward along the 30° curve and read the tilt it leaves
        lead = family[1]
        lead_tilt = np.array(lead["tilt_deg"])
        slide = ValueTracker(a_tilde[0])

        def tilt_at():
            return float(np.interp(slide.get_value(), a_tilde, lead_tilt))

        bead = always_redraw(lambda: Dot(ax3.c2p(slide.get_value(), min(tilt_at(), top3)),
                                         radius=0.1, color=P.FG).set_z_index(6))

        def reading():
            g = VGroup(
                layout.label(f"a√(1−e²) = {slide.get_value():.0f} AU", font_size=17),
                layout.label(f"tilt {tilt_at():.1f}°", font_size=22, color=P.ORANGE,
                             weight="BOLD"),
            ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
            return g.move_to(ax3.c2p(680, 7.1))

        live = always_redraw(reading)
        cap4 = layout.caption(
            "Too close, the tilt swings back and forth; too far, the pull is too weak",
            font_size=22)
        self.play(FadeIn(bead), FadeIn(live), run_time=0.5)
        self.play(slide.animate.set_value(a_tilde[-1]), rate_func=lambda s: s, run_time=5.5)
        self.play(FadeOut(bead), FadeOut(live), FadeOut(cap3), FadeIn(cap4), run_time=0.5)
        self.play(FadeIn(band), FadeIn(band_lab), FadeIn(marks), run_time=0.8)
        timing.hold_to_read(self, cap4, band_lab, settle=1.0)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, f"A formula replaces the integrations: {observed:.0f}° needs a√(1−e²) of "
                  f"{lo_a:.0f} to {hi_a:.0f} AU.")
