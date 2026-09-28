"""Beust (2016) -- orbital clustering by Planet Nine: secular or resonant?

Beust redoes the secular model behind the clustering without expanding the
Hamiltonian in the ratio of semi-major axes, because for orbits a few hundred AU
out that series does not converge. The scene first sets the truncated apsidal
forcing against the full orbit average, then draws the phase portrait of the
full Hamiltonian and follows one orbit along a level curve of the family whose
apse never lines up with the planet's. Reproduced in p9-2016-secular-resonance;
every curve is the crate's own (anim.json -> papers -> p9-2016-secular-resonance).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Create,
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
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2016-secular-resonance"


def polyline(ax, points, color, stroke_width=2.0, opacity=1.0):
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners([ax.c2p(x, y) for x, y in points])
    return m.set_stroke(opacity=opacity)


def arc_fractions(points):
    """Cumulative arc length along a (Δϖ/360, e) curve, normalised to [0, 1]."""
    xy = np.array([[p[0] / 360.0, p[1]] for p in points])
    s = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(xy, axis=0), axis=1))])
    return s / s[-1]


class SecularResonance2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        h = d["harmonic"]
        port = d["portrait"]
        planet = d["planet"]

        self.add(paper.scene_header(CRATE))

        # 1. the truncated series against the full orbit average
        ax, labels = widgets.labeled_axes(
            [0, 1, 0.2], [0, 1, 0.25], x_label="eccentricity of the distant orbit",
            y_label="strength of the apsidal forcing", y_rotate=True, numbers=True,
            x_length=9.5, y_length=4.3, shift_down=-0.35)
        full = widgets.curve(ax, h["e"], [abs(v) for v in h["exact"]], color=P.ORANGE)
        cut = widgets.curve(ax, h["e"], [abs(v) for v in h["truncated"]], color=P.RED)
        key = VGroup(
            layout.label("full orbit average, no expansion", font_size=18, color=P.ORANGE),
            layout.label("series cut at the octupole term", font_size=18, color=P.RED),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        key.move_to(ax.c2p(0.3, 0.85))
        cap = layout.caption(
            f"Orbit at {h['a_au']:.0f} AU, planet at {planet['a_au']:.0f} AU: "
            "the series is not converging", font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(Create(cut), FadeIn(key[1]), run_time=1.4)
        self.play(Create(full), FadeIn(key[0]), FadeIn(cap), run_time=1.4)
        tally = paper.result_readout(
            "error of the truncated series", f"up to {100 * d['truncation_error']:.0f}%",
            color=P.RED).scale(0.75)
        tally.move_to(ax.c2p(0.8, 0.22))
        self.play(FadeIn(tally))
        timing.hold_to_read(self, cap, tally, settle=1.0)
        self.play(FadeOut(VGroup(ax, labels, full, cut, key, cap, tally)))

        # 2. the phase portrait of the full Hamiltonian
        ax2, labels2 = widgets.labeled_axes(
            [0, 360, 90], [0, 1, 0.25],
            x_label="apsidal angle from the planet, Δϖ (deg)",
            y_label="eccentricity", y_rotate=True, numbers=True,
            x_length=8.6, y_length=4.1, shift_down=-0.55)
        VGroup(ax2, labels2).shift(LEFT * 1.9)
        # "anti": spans 180° and never comes within 45° of alignment
        fam = [ln for ln in port["lines"] if ln["kind"] == "anti"]
        anti = VGroup(*[polyline(ax2, ln["points"], P.ORANGE, 2.4) for ln in fam])
        rest = VGroup(*[polyline(ax2, ln["points"], P.FG, 1.6, 0.55)
                        for ln in port["lines"] if ln["kind"] != "anti"])
        cap2 = layout.caption(
            f"Full Hamiltonian at {port['a_au']:.0f} AU: an orbit's (Δϖ, e) slides along one curve",
            font_size=22)
        self.play(Create(ax2), FadeIn(labels2))
        self.play(LaggedStart(*[Create(m) for m in [*rest, *anti]], lag_ratio=0.03),
                  FadeIn(cap2), run_time=3.0)
        timing.hold_to_read(self, cap2, settle=0.8)

        # 3. the anti-aligned family, and one orbit followed along its curve
        cap3 = layout.caption(
            "On the orange curves Δϖ never reaches 0°: the orbit stays opposite the planet",
            font_size=22)
        self.play(rest.animate.set_stroke(opacity=0.2), anti.animate.set_stroke(width=3.2),
                  FadeOut(cap2), FadeIn(cap3), run_time=1.2)
        timing.hold_to_read(self, cap3, settle=0.6)

        path = min(fam, key=lambda ln: abs(min(p[1] for p in ln["points"]) - 0.5))["points"]
        s_of = arc_fractions(path)

        def at(s):
            return (float(np.interp(s, s_of, [p[0] for p in path])),
                    float(np.interp(s, s_of, [p[1] for p in path])))

        # top view: planet apse along +x, particle apse at angle Δϖ
        a_p, e_p, a = planet["a_au"], planet["e"], port["a_au"]
        lo = -a_p * (1 + e_p)
        hi = a * 1.96
        scale = 3.6 / (hi - lo)
        sun_at = np.array([4.8 - 0.5 * (hi + lo) * scale, 0.75, 0.0])
        sun = orbits.sun(radius=0.06).move_to(sun_at)
        p9 = orbits.ellipse_orbit(a_p * scale, e_p, color=P.BLUE, stroke_width=2.2).shift(sun_at)
        p9_lab = layout.label("Planet Nine", font_size=14, color=P.BLUE)
        p9_lab.next_to(p9, UP, buff=0.08)
        s = ValueTracker(0.02)

        def body():
            w, e = at(s.get_value())
            return orbits.ellipse_orbit(a * scale, min(e, 0.97), color=P.ORANGE,
                                        varpi=np.radians(w), stroke_width=2.4).shift(sun_at)

        def pointer():
            w, e = at(s.get_value())
            return Dot(ax2.c2p(w, e), radius=0.08, color=P.ORANGE).set_z_index(6)

        def readout():
            w, e = at(s.get_value())
            return layout.label(f"Δϖ = {w:5.0f}°    e = {e:.2f}", font_size=16,
                                color=P.ORANGE).move_to([4.8, -0.9, 0])

        orbit_m = always_redraw(body)
        dot = always_redraw(pointer)
        read = always_redraw(readout)
        view_lab = layout.label(f"top view, orbit at {a:.0f} AU", font_size=14, color=P.MUTED)
        view_lab.move_to([4.8, -1.22, 0])
        cap4 = layout.caption(
            "Follow one: the apse swings across 180° as e dips and climbs, never lining up",
            font_size=22)
        self.play(FadeIn(sun), Create(p9), FadeIn(p9_lab), FadeIn(orbit_m), FadeIn(dot),
                  FadeIn(read), FadeIn(view_lab), FadeOut(cap3), FadeIn(cap4))
        self.play(s.animate.set_value(0.98), run_time=6.0, rate_func=rate_functions.linear)
        self.play(s.animate.set_value(0.5), run_time=3.0, rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap4, settle=0.4)

        tally = paper.result_readout(
            "widest anti-aligned swing", f"±{port['anti_halfwidth_deg']:.0f}° about 180°",
            color=P.ORANGE).scale(0.62)
        tally.move_to([4.8, -2.02, 0])
        self.play(FadeOut(cap4), FadeIn(tally))
        timing.hold_to_read(self, tally, settle=1.0)

        layout.show_takeaway(
            self, "Secular dynamics alone, no resonance, can keep orbits anti-aligned.")
