"""Bailey, Batygin & Brown (2016) -- solar obliquity induced by Planet Nine.

The Sun's equator is tilted six degrees to the planets' plane. The paper
integrates the secular spin-orbit equations for 4.5 Gyr and asks, over a grid
of Planet Nine mass, semi-major axis and eccentricity, what inclination the
planet needs to produce that tilt from an aligned start.

Everything drawn comes from anim.json -> papers -> p9-2016-obliquity: the
integrated tilt histories and the solved inclination in every grid cell, both
from the reproduction crate's integrator.
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
    ManimColor,
    Rectangle,
    Scene,
    SurroundingRectangle,
    ValueTracker,
    VGroup,
    always_redraw,
    interpolate_color,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2016-obliquity"

CELL_W, CELL_H = 0.6, 0.56


def survey_panel(survey, a_grid, e_grid, lo, hi, show_e):
    """One mass of the survey: a (across) by e (up), each cell carrying the
    inclination that yields the observed tilt."""
    g = VGroup()
    g.boxes = {}
    cells = {(c["a_au"], round(c["e"], 2)): c for c in survey["cells"]}
    for ix, a in enumerate(a_grid):
        for iy, e in enumerate(e_grid):
            c = cells.get((a, round(e, 2)))
            if c is None:
                continue
            box = Rectangle(width=CELL_W * 0.94, height=CELL_H * 0.92, stroke_width=0)
            box.move_to([ix * CELL_W, iy * CELL_H, 0])
            g.boxes[(a, round(e, 2))] = box
            need = c["required_i_deg"]
            if need is None:
                box.set_fill(P.RED, opacity=0.3)
                g.add(box)
                continue
            shade = interpolate_color(ManimColor(P.TEAL), ManimColor(P.ORANGE),
                                      (need - lo) / (hi - lo))
            box.set_fill(shade, opacity=0.85)
            num = layout.label(f"{need:.0f}", font_size=13, color=P.BG, weight="BOLD")
            num.move_to(box)
            g.add(box, num)
    for ix, a in enumerate(a_grid):
        g.add(layout.label(f"{a:.0f}", font_size=13, color=P.FG)
              .move_to([ix * CELL_W, -0.62 * CELL_H - 0.1, 0]).set_opacity(0.75))
    if show_e:
        for iy, e in enumerate(e_grid):
            g.add(layout.label(f"{e:.1f}", font_size=13, color=P.FG)
                  .move_to([-0.62 * CELL_W - 0.14, iy * CELL_H, 0]).set_opacity(0.75))
    title = layout.label(f"{survey['mass_earth']:.0f} M⊕", font_size=17, color=P.BLUE,
                         weight="BOLD")
    title.move_to([0.5 * (len(a_grid) - 1) * CELL_W, (len(e_grid) - 0.5) * CELL_H + 0.18, 0])
    g.add(title)
    return g


class Obliquity2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        show, nominal, pub = d["showcase"], d["nominal"], d["published"]
        observed = d["observed_obliquity_deg"]
        age = d["age_gyr"]

        self.add(paper.scene_header(CRATE))

        # 1. the tilt grows over the age of the Solar System
        top = 10.0
        ax, labels = widgets.labeled_axes(
            [0, age, 0.5], [0, top, 2], x_label="time since the Solar System formed  (Gyr)",
            y_label="tilt of the Sun's equator to the planets' plane  (deg)", y_rotate=True,
            numbers=True, x_length=9.0, y_length=4.3, shift_down=-0.2, font_size=22)
        VGroup(ax, labels).shift(LEFT * 1.3)
        target = DashedLine(ax.c2p(0, observed), ax.c2p(age, observed), color=P.GREEN,
                            stroke_width=2.5)
        target_lab = layout.label(f"observed today: {observed:.0f}°", font_size=14,
                                  color=P.GREEN)
        target_lab.next_to(ax.c2p(0, observed), UP, buff=0.08).shift(RIGHT * 1.15)

        curves, tags = VGroup(), VGroup()
        for h in show["histories"]:
            curves.add(widgets.curve(ax, h["t_gyr"], h["obliquity_deg"], color=P.ORANGE,
                                     stroke_width=3.2))
            tags.add(layout.label(
                f"inclined {h['i9_deg']:.0f}°:  {h['final_obliquity_deg']:.1f}°",
                font_size=14, color=P.ORANGE)
                .next_to(ax.c2p(age, h["final_obliquity_deg"]), RIGHT, buff=0.12))
        head = layout.label(
            f"{show['mass_earth']:.0f} M⊕ at {show['a_au']:.0f} AU, e = {show['e']:.1f}",
            font_size=15, color=P.ORANGE, weight="BOLD")
        head.next_to(tags, UP, buff=0.3, aligned_edge=LEFT)

        # the tilt itself, drawn to scale: Sun's equator against the planets' plane,
        # driven by the 20° history as it is integrated
        lead = show["histories"][1]
        t_now = ValueTracker(0.0)

        def tilt_now():
            return float(np.interp(t_now.get_value(), lead["t_gyr"], lead["obliquity_deg"]))

        hub = np.array([4.9, -2.3, 0.0])
        sun = Dot(hub, radius=0.13, color=P.SUN).set_z_index(3)
        equator = Line(hub + LEFT * 1.3, hub + RIGHT * 1.3, color=P.SUN, stroke_width=2.5)
        eq_lab = layout.label("Sun's equator", font_size=14, color=P.SUN)
        eq_lab.next_to(hub + RIGHT * 1.3, DOWN, buff=0.1).align_to(hub + RIGHT * 1.3, RIGHT)
        plane = always_redraw(lambda: Line(hub + LEFT * 1.3, hub + RIGHT * 1.3, color=P.FG,
                                           stroke_width=2.5).rotate(np.radians(tilt_now()),
                                                                    about_point=hub))
        plane_lab = always_redraw(lambda: layout.label(
            f"planets' plane  {tilt_now():.1f}°", font_size=14, color=P.FG)
            .next_to(hub + RIGHT * 1.3, UP, buff=0.12 + 1.3 * np.tan(np.radians(tilt_now())))
            .align_to(hub + RIGHT * 1.3, RIGHT))
        tracer = always_redraw(lambda: Dot(ax.c2p(t_now.get_value(), tilt_now()), radius=0.07,
                                           color=P.FG).set_z_index(5))
        inset = VGroup(equator, eq_lab, sun)

        cap = layout.caption(
            "Integrated for 4.5 Gyr from an aligned start: the planet tips the planets' plane",
            font_size=22)
        self.play(Create(ax), FadeIn(labels), run_time=1.0)
        self.play(Create(target), FadeIn(target_lab), FadeIn(inset), FadeIn(plane),
                  FadeIn(plane_lab), FadeIn(cap), run_time=0.8)
        self.add(tracer)
        self.play(Create(curves[1]), t_now.animate.set_value(age), FadeIn(tags[1]),
                  rate_func=lambda s: s, run_time=4.0)
        self.play(Create(curves[0]), Create(curves[2]), FadeIn(head), FadeIn(tags[0]),
                  FadeIn(tags[2]), run_time=1.4)
        timing.hold_to_read(self, cap, tags, settle=0.8)

        lone = widgets.curve(ax, nominal["t_gyr"], nominal["obliquity_deg"], color=P.BLUE,
                             stroke_width=3.2)
        lo_t, hi_t = pub["nominal_tilt_deg"]
        lone_tag = layout.label(
            f"nominal {nominal['mass_earth']:.0f} M⊕, {nominal['a_au']:.0f} AU,\n"
            f"inclined {nominal['i9_deg']:.0f}°:  {nominal['final_obliquity_deg']:.1f}°\n"
            f"paper: {lo_t:.0f}° to {hi_t:.0f}°",
            font_size=14, color=P.BLUE, line_spacing=0.9)
        lone_tag.next_to(ax.c2p(age, nominal["final_obliquity_deg"]), RIGHT, buff=0.12)
        lone_tag.shift(DOWN * 0.45)
        cap2 = layout.caption(
            "The nominal 700 AU planet falls short: the tilt favours a closer or heavier one",
            font_size=22)
        self.play(Create(lone), FadeIn(lone_tag), FadeOut(cap), FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, lone_tag, settle=1.0)
        self.play(FadeOut(VGroup(ax, labels, target, target_lab, curves, tags, head, lone,
                                 lone_tag, cap2, inset, plane, plane_lab, tracer)))

        # 2. the survey: inclination needed in every grid cell
        lo, hi = d["required_i_min_deg"], d["required_i_max_deg"]
        a_grid = d["a_grid_au"]
        used = {round(c["e"], 2) for s in d["surveys"] for c in s["cells"]}
        e_grid = [e for e in d["e_grid"] if round(e, 2) in used]
        panels = VGroup(*[
            survey_panel(s, a_grid, e_grid, lo, hi, show_e=(k == 0))
            for k, s in enumerate(d["surveys"])
        ]).arrange(RIGHT, buff=0.55, aligned_edge=DOWN)
        panels.move_to([0.25, 0.55, 0])
        xlab = layout.label("semi-major axis  (AU)", font_size=15)
        xlab.next_to(panels, DOWN, buff=0.15)
        ylab = layout.label("eccentricity", font_size=15).rotate(np.pi / 2)
        ylab.next_to(panels, LEFT, buff=0.15)

        bar = VGroup(*[
            Rectangle(width=0.16, height=0.22, stroke_width=0).set_fill(
                interpolate_color(ManimColor(P.TEAL), ManimColor(P.ORANGE), k / 11.0),
                opacity=0.85)
            for k in range(12)
        ]).arrange(RIGHT, buff=0)
        key = VGroup(
            layout.label(f"inclination needed:  {lo:.0f}°", font_size=14),
            bar,
            layout.label(f"{hi:.0f}°", font_size=14),
            Rectangle(width=0.3, height=0.22, stroke_width=0).set_fill(P.RED, opacity=0.3),
            layout.label("cannot reach 6° in 4.5 Gyr", font_size=14, color=P.RED),
        ).arrange(RIGHT, buff=0.15)
        key[3].shift(RIGHT * 0.5)
        key[4].shift(RIGHT * 0.5)
        key.next_to(xlab, DOWN, buff=0.22)

        cap3 = layout.caption(
            f"{d['cells_total']} orbits with perihelia near 250 AU, each solved for the "
            "inclination that gives 6°", font_size=22)
        self.play(FadeIn(panels, lag_ratio=0.02), FadeIn(xlab), FadeIn(ylab), FadeIn(cap3),
                  run_time=2.0)
        self.play(FadeIn(key), run_time=0.6)
        timing.hold_to_read(self, cap3, key, settle=1.0)

        mass_of = [s["mass_earth"] for s in d["surveys"]]
        show_box = panels[mass_of.index(show["mass_earth"])].boxes[
            (show["a_au"], round(show["e"], 2))]
        nom_box = panels[mass_of.index(nominal["mass_earth"])].boxes[
            (nominal["a_au"], round(nominal["e"], 2))]
        ring_s = SurroundingRectangle(show_box, color=P.FG, buff=0.03, stroke_width=3)
        ring_n = SurroundingRectangle(nom_box, color=P.BLUE, buff=0.03, stroke_width=3)
        cap_s = layout.caption(
            f"White box: the orbit integrated above needs {show['required_i_deg']:.0f}°",
            font_size=22)
        self.play(Create(ring_s), FadeOut(cap3), FadeIn(cap_s), run_time=0.8)
        timing.hold_to_read(self, cap_s, settle=0.6)
        cap_n = layout.caption(
            "Blue box: Batygin & Brown's 700 AU planet never reaches 6° in 4.5 Gyr",
            font_size=22)
        self.play(Create(ring_n), FadeOut(cap_s), FadeIn(cap_n), run_time=0.8)
        timing.hold_to_read(self, cap_n, settle=0.8)
        cap3 = cap_n

        need_lo, need_hi = pub["required_i_deg"]
        cap4 = layout.caption(
            f"Reproduced: {lo:.0f}° to {hi:.0f}° where a solution exists   "
            f"(paper: {need_lo:.0f}° to {need_hi:.0f}°)", font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        timing.hold_to_read(self, cap4, settle=1.2)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, f"A 6° solar tilt needs Planet Nine inclined {lo:.0f}°-{hi:.0f}°: "
                  "the Sun constrains the orbit.")
