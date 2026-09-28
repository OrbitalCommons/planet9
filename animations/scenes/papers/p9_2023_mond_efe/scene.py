"""Brown & Mathur (2023) -- MOND as an alternative to Planet Nine.

In MOND the Galaxy's own field does not cancel out of the Solar System: it
leaves a quadrupole tide whose axis points at the Galactic centre. Averaged
over an eccentric orbit, that tide favours long axes lying along the
Galactic-centre line, and traps orbits whose axis starts close enough to it.
The scene draws the crate's tidal field, sweeps an orbit's apse through the
crate's orbit-averaged disturbing function, and then puts the ten Brown (2017)
objects against the predicted axis. Data: anim.json -> papers -> p9-2023-mond-efe.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    AnnularSector,
    Arrow,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    GrowArrow,
    LaggedStart,
    Polygon,
    Scene,
    Square,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2023-mond-efe"


def unit(deg):
    r = np.radians(deg)
    return np.array([np.cos(r), np.sin(r), 0.0])


class MondEfe2023(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        gc = d["gc_lon_deg"]
        lobe = d["lobe_half_width_deg"]

        self.add(paper.scene_header(CRATE))

        # 1. the tide: the crate's EFE quadrupole field in the ecliptic plane
        centre = np.array([-3.6, -0.1, 0.0])
        half = d["field_half_au"]
        s = 2.3 / half
        frame = Square(side_length=2 * half * s, color=P.MUTED, stroke_width=1.2)
        frame.move_to(centre).set_stroke(opacity=0.5)
        arrows = VGroup()
        for f in d["field"]:
            p = centre + s * np.array([f["x_au"], f["y_au"], 0.0])
            v = 0.52 * np.array([f["ax"], f["ay"], 0.0])
            arrows.add(Arrow(p - 0.5 * v, p + 0.5 * v, buff=0, color=P.ORANGE, stroke_width=3.2,
                             max_tip_length_to_length_ratio=0.45, tip_length=0.13))
        sun = orbits.sun(radius=0.06).move_to(centre)
        nep = Circle(radius=30.0 * s, color=P.FG, stroke_width=1.4).move_to(centre)
        nep_lab = layout.label("Neptune's orbit (30 AU)", font_size=14, color=P.FG)
        nep_lab.next_to(centre, RIGHT, buff=0.2).shift(UP * 0.02)
        gc_start = np.array([-0.55, -0.55, 0.0])
        gc_arrow = Arrow(gc_start, gc_start + unit(gc) * 1.3, buff=0, color=P.PURPLE,
                         stroke_width=3)
        gc_lab = layout.label("toward the Galactic centre", font_size=15, color=P.PURPLE)
        gc_lab.next_to(gc_arrow.get_end(), RIGHT, buff=0.15)
        scale_lab = layout.label(f"{2 * half:.0f} AU across, ecliptic plane seen from above",
                                 font_size=14, color=P.MUTED)
        scale_lab.next_to(frame, UP, buff=0.3)

        facts = VGroup(
            layout.label("MOND: gravity changes below", font_size=18),
            layout.label(f"a₀ = {d['a0_m_s2'] * 1e10:.1f}×10⁻¹⁰ m/s²", font_size=18, color=P.ORANGE),
            layout.label(f"The Sun's pull falls to a₀ at {d['mond_radius_au']:,.0f} AU", font_size=16,
                         color=P.MUTED),
            layout.label("The Galaxy's field does not cancel out:", font_size=16, color=P.MUTED),
            layout.label("it leaves a tide on the outer Solar System", font_size=16,
                         color=P.MUTED),
            layout.label("stretching along the Galactic-centre line,", font_size=16,
                         color=P.ORANGE),
            layout.label("squeezing across it", font_size=16, color=P.ORANGE),
        ).arrange(DOWN, buff=0.17, aligned_edge=LEFT).move_to([3.6, 0.75, 0])
        cap = layout.caption("The tide MOND predicts, computed on a grid around the Sun",
                             font_size=22)
        self.play(FadeIn(frame), FadeIn(sun), FadeIn(scale_lab), FadeIn(cap), run_time=0.8)
        self.play(Create(nep), FadeIn(nep_lab), FadeIn(facts[:3]), run_time=0.8)
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows], lag_ratio=0.02),
                  FadeIn(facts[3:]), run_time=2.2)
        self.play(GrowArrow(gc_arrow), FadeIn(gc_lab), run_time=0.8)
        timing.hold_to_read(self, cap, facts, settle=0.8)
        self.play(FadeOut(VGroup(frame, arrows, sun, nep, nep_lab, gc_arrow, gc_lab, scale_lab,
                                 facts, cap)), run_time=0.7)

        # 2. averaged over an orbit: which way should the long axis point?
        o = np.array([-4.3, 0.0, 0.0])
        a_s = 1.25
        axis = VGroup(
            DashedLine(o - unit(gc) * 2.5, o + unit(gc) * 2.5, color=P.PURPLE, stroke_width=2),
        )
        axis_lab = layout.label("Galactic-centre line", font_size=15, color=P.PURPLE)
        axis_lab.next_to(o + unit(gc) * 2.5, RIGHT, buff=0.1).shift(UP * 0.1)
        sun2 = orbits.sun(radius=0.07).move_to(o)
        varpi = ValueTracker(0.0)
        orbit = always_redraw(lambda: orbits.ellipse_orbit(
            a_s, d["e"], color=P.GREEN, varpi=np.radians(varpi.get_value()),
            stroke_width=2.6).shift(o))
        orbit_lab = layout.label(f"a test orbit: a = {d['a_au']:.0f} AU, e = {d['e']:.1f}",
                                 font_size=15, color=P.GREEN)
        orbit_lab.move_to(o + DOWN * 2.75)

        lon = np.array(d["lon_deg"])
        rbar = np.array(d["rbar"])
        top = float(np.ceil(rbar.max() * 10) / 10)
        bot = float(np.floor(rbar.min() * 10) / 10)
        ax, labs = widgets.labeled_axes(
            [0, 360, 90], [bot, top, 0.1], x_label="longitude of the orbit's long axis  (deg)",
            y_label="tidal energy, orbit-averaged", y_rotate=True, numbers=False,
            x_length=6.4, y_length=3.6, shift_down=0)
        ax.move_to([2.2, 0.35, 0])
        labs[0].next_to(ax, DOWN, buff=0.45)
        labs[1].next_to(ax, LEFT, buff=0.15)
        ticks = VGroup(*[layout.label(f"{x}", font_size=14, color=P.MUTED)
                         .next_to(ax.c2p(x, bot), DOWN, buff=0.12) for x in (0, 90, 180, 270, 360)])
        curve = widgets.curve(ax, lon, rbar, color=P.ORANGE)
        circ = DashedLine(ax.c2p(0, d["rbar_circular"]), ax.c2p(360, d["rbar_circular"]),
                          color=P.MUTED, stroke_width=1.5)
        circ_lab = layout.label("circular orbit", font_size=13, color=P.MUTED)
        circ_lab.next_to(ax.c2p(360, d["rbar_circular"]), RIGHT, buff=0.08)
        gc_marks = VGroup()
        for x in (gc - 180.0, gc):
            gc_marks.add(DashedLine(ax.c2p(x, bot), ax.c2p(x, top), color=P.PURPLE,
                                    stroke_width=1.6))
        gc_tag = layout.label(f"Galactic centre {gc:.0f}°", font_size=13, color=P.PURPLE)
        gc_tag.next_to(ax.c2p(gc, top), UP, buff=0.08)
        anti_tag = layout.label(f"{gc - 180:.0f}°", font_size=13, color=P.PURPLE)
        anti_tag.next_to(ax.c2p(gc - 180, top), UP, buff=0.08)
        rider = always_redraw(lambda: Dot(
            ax.c2p(varpi.get_value(), float(np.interp(varpi.get_value(), lon, rbar))),
            radius=0.07, color=P.GREEN))

        cap2 = layout.caption("Turn one orbit's long axis all the way round and add up the tide",
                              font_size=22)
        self.play(FadeIn(sun2), Create(axis), FadeIn(axis_lab), FadeIn(orbit), FadeIn(orbit_lab),
                  FadeIn(ax), FadeIn(labs), FadeIn(ticks), FadeIn(cap2), run_time=1.0)
        self.play(FadeIn(gc_marks), FadeIn(gc_tag), FadeIn(anti_tag), FadeIn(rider), run_time=0.6)
        trace = always_redraw(lambda: widgets.curve(
            ax, lon[lon <= varpi.get_value() + 1e-9], rbar[lon <= varpi.get_value() + 1e-9],
            color=P.ORANGE) if varpi.get_value() > 2 else VGroup())
        self.add(trace)
        self.play(varpi.animate.set_value(360.0), run_time=6.0, rate_func=linear)
        self.remove(trace)
        self.add(curve)
        timing.hold_to_read(self, cap2, settle=0.4)

        lobes = VGroup()
        for c in (gc - 180.0, gc):
            lo, hi = c - lobe, c + lobe
            a_, b_ = ax.c2p(lo, bot), ax.c2p(hi, top)
            band = Polygon(a_, [b_[0], a_[1], 0], b_, [a_[0], b_[1], 0], stroke_width=0)
            lobes.add(band.set_fill(P.ORANGE, opacity=0.12))
        cap3 = layout.caption(
            f"Axes within ±{lobe:.0f}° of the Galactic line are trapped: the tide swings them back",
            font_size=22)
        self.play(FadeIn(lobes), Create(circ), FadeIn(circ_lab), FadeOut(cap2), FadeIn(cap3),
                  varpi.animate.set_value(gc - 180.0 + 35.0), run_time=1.2)
        self.play(varpi.animate.set_value(gc - 180.0 - 35.0), run_time=1.6)
        self.play(varpi.animate.set_value(gc - 180.0 + 20.0), run_time=1.3)
        self.play(varpi.animate.set_value(gc - 180.0), run_time=0.8)
        timing.hold_to_read(self, cap3, settle=0.8)
        self.play(FadeOut(VGroup(sun2, axis, axis_lab, orbit, orbit_lab, ax, labs, ticks, curve,
                                 circ, circ_lab, gc_marks, gc_tag, anti_tag, rider, lobes, cap3)),
                  run_time=0.7)

        # 3. the test: the ten Brown (2017) objects against the predicted axis
        c3 = np.array([-3.0, 0.2, 0.0])
        R = 2.35
        dial = Circle(radius=R, color=P.MUTED, stroke_width=1.4).move_to(c3)
        lobes3 = VGroup(*[
            AnnularSector(inner_radius=0.0, outer_radius=R, angle=np.radians(2 * lobe),
                          start_angle=np.radians(c - lobe), arc_center=c3, color=P.ORANGE,
                          fill_opacity=0.13, stroke_width=0)
            for c in (gc - 180.0, gc)])
        gc_line = DashedLine(c3 - unit(gc) * R, c3 + unit(gc) * R, color=P.PURPLE,
                             stroke_width=2)
        lon_ticks = VGroup()
        for x in (0, 90, 180, 270):
            lon_ticks.add(layout.label(f"{x}°", font_size=13, color=P.MUTED)
                          .move_to(c3 + unit(x) * (R + 0.28)))
        gc_tag3 = layout.label("Galactic centre", font_size=14, color=P.PURPLE)
        gc_tag3.next_to(c3 + unit(gc) * R, RIGHT, buff=0.3).shift(UP * 0.12)
        objs = d["objects"]
        amax = max(ob["a_au"] for ob in objs)
        arrows3 = VGroup()
        for ob in objs:
            length = R * (0.55 + 0.4 * ob["a_au"] / amax)
            col = P.GREEN if ob["in_lobe"] else P.RED
            arrows3.add(Arrow(c3, c3 + unit(ob["varpi_deg"]) * length, buff=0, color=col,
                              stroke_width=2.6, max_tip_length_to_length_ratio=0.08))
        sun3 = orbits.sun(radius=0.08).move_to(c3)

        stats = VGroup(
            layout.label("ten distant objects (Brown 2017)", font_size=18, weight="BOLD"),
            layout.label("arrow = direction of perihelion", font_size=15, color=P.MUTED),
            layout.label(f"inside the trap: {d['n_in_lobe']} of {d['n_objects']}", font_size=18,
                         color=P.GREEN),
            layout.label(f"outside it: {d['n_objects'] - d['n_in_lobe']}", font_size=16,
                         color=P.RED),
            layout.label(f"by chance: {100 * d['lobe_chance_fraction']:.0f}%", font_size=16,
                         color=P.MUTED),
            layout.label(f"mean offset from the line: {d['mean_line_sep_deg']:.0f}°",
                         font_size=18, color=P.GREEN),
            layout.label(f"random directions: {d['random_line_sep_deg']:.0f}°", font_size=16,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT).move_to([3.3, 0.6, 0])
        cap4 = layout.caption("MOND's prediction: long axes along the Galactic-centre line",
                              font_size=22)
        self.play(FadeIn(dial), FadeIn(lon_ticks), Create(gc_line), FadeIn(gc_tag3),
                  FadeIn(lobes3), FadeIn(sun3), FadeIn(cap4), run_time=1.0)
        timing.hold_to_read(self, cap4, settle=0.3)
        cap5 = layout.caption("The observed perihelion directions, computed against that axis",
                              font_size=22)
        self.play(LaggedStart(*[GrowArrow(a) for a in arrows3], lag_ratio=0.12),
                  FadeIn(stats[:2]), FadeOut(cap4), FadeIn(cap5), run_time=2.0)
        self.play(FadeIn(stats[2:5]), run_time=0.6)
        self.play(FadeIn(stats[5:]), run_time=0.6)
        readout = paper.result_readout("mean apsidal offset", f"{d['mean_line_sep_deg']:.0f}°",
                                       color=P.GREEN).scale(0.62)
        readout.next_to(stats, DOWN, buff=0.4)
        self.play(FadeIn(readout), run_time=0.5)
        timing.hold_to_read(self, cap5, stats, settle=1.0)
        self.play(FadeOut(cap5), run_time=0.4)

        layout.show_takeaway(
            self, f"A testable axis: the objects lean toward it, {d['mean_line_sep_deg']:.0f}° vs "
                  f"{d['random_line_sep_deg']:.0f}° by chance.")
