"""Madigan & McCourt (2015) -- the inclination instability.

A flat disc of eccentric orbits is unstable under its own gravity: inclinations
grow exponentially and the disc opens into a cone, with no planet involved. The
price is mass. The scene shows the mechanism (paced by the crate's e-folding
time), then the two computed results: how long an e-fold takes against disc
mass, and how much mass the disc needs before its self-gravity beats the
differential precession forced by the giant planets. Everything plotted is from
anim.json -> papers -> p9-2016-inclination-instability.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    Line,
    RIGHT,
    UP,
    Axes,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Polygon,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing

CRATE = "p9-2016-inclination-instability"


class LogAxes(VGroup):
    """Log-log axes between arbitrary limits, with labelled ticks at the given
    values. ``p(x, y)`` maps data to the scene; ``curve`` draws a polyline."""

    def __init__(self, x_lim, y_lim, x_ticks, y_ticks, x_label, y_label,
                 x_length=9.6, y_length=4.4, centre=(0.3, 0.35, 0)):
        super().__init__()
        self.x0, self.y0 = np.log10(x_lim[0]), np.log10(y_lim[0])
        self.ax = Axes(x_range=[0, np.log10(x_lim[1]) - self.x0, 10],
                       y_range=[0, np.log10(y_lim[1]) - self.y0, 10],
                       x_length=x_length, y_length=y_length,
                       axis_config={"color": P.MUTED, "include_tip": False,
                                    "include_ticks": False})
        self.ax.move_to(centre)
        self.add(self.ax)
        for v in x_ticks:
            at = self.p(v, y_lim[0])
            self.add(Line(at, at + DOWN * 0.08, color=P.MUTED, stroke_width=1.5))
            self.add(layout.label(f"{v:g}", font_size=15, color=P.MUTED)
                     .next_to(at, DOWN, buff=0.14))
        widest = 0.0
        for v in y_ticks:
            at = self.p(x_lim[0], v)
            self.add(Line(at, at + LEFT * 0.08, color=P.MUTED, stroke_width=1.5))
            lab = layout.label(f"{v:g}", font_size=15, color=P.MUTED).next_to(at, LEFT, buff=0.14)
            widest = max(widest, lab.width)
            self.add(lab)
        self.add(layout.label(x_label, font_size=17).next_to(self.ax, DOWN, buff=0.45))
        self.add(layout.label(y_label, font_size=17).rotate(np.pi / 2)
                 .next_to(self.ax, LEFT, buff=0.3 + widest))

    def p(self, x, y):
        return self.ax.c2p(np.log10(x) - self.x0, np.log10(y) - self.y0)

    def curve(self, xs, ys, color, stroke_width=3.0):
        m = VMobject(color=color, stroke_width=stroke_width)
        m.set_points_as_corners([self.p(x, y) for x, y in zip(xs, ys)])
        return m

    def band(self, x0, x1, y0, y1, color, opacity=0.13):
        a, b = self.p(x0, y0), self.p(x1, y1)
        r = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0, color=color)
        return r.set_fill(color, opacity=opacity)


def tilted_orbit(a, e, varpi, tilt, elevation, color):
    """An orbit pitched about its own latus rectum by ``tilt`` (aphelion up),
    turned to apsidal direction ``varpi``, seen from ``elevation`` above the
    disc plane."""
    nu = np.linspace(0, 2 * np.pi, 72)
    r = a * (1 - e * e) / (1 + e * np.cos(nu))
    x, y = r * np.cos(nu), r * np.sin(nu)
    # pitch about the y axis of the orbit frame: aphelion (x < 0) rises
    xp, zp = x * np.cos(tilt), -x * np.sin(tilt)
    X = xp * np.cos(varpi) - y * np.sin(varpi)
    Y = xp * np.sin(varpi) + y * np.cos(varpi)
    pts = np.column_stack([X, Y * np.sin(elevation) + zp * np.cos(elevation), np.zeros_like(X)])
    m = VMobject(color=color, stroke_width=1.8)
    m.set_points_as_corners(pts)
    return m.set_stroke(opacity=0.8)


class InclinationInstability2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        fid = d["fiducial"]
        tau = fid["tau_myr"]
        m_lo, m_hi = d["paper_mass_range_earth"]

        self.add(paper.scene_header(CRATE))

        # 1. mechanism: a flat eccentric disc opens into a cone, one e-fold per tau
        centre = np.array([-1.6, 0.0, 0.0])
        elev = np.radians(14)
        varpis = np.linspace(0, 2 * np.pi, 12, endpoint=False)
        i0 = 1.0
        efolds = ValueTracker(0.0)

        def disc():
            tilt = np.radians(i0 * np.exp(efolds.get_value()))
            g = VGroup(*[tilted_orbit(1.9, fid["e"], v, tilt, elev, P.GREEN) for v in varpis])
            return g.shift(centre)

        live = always_redraw(disc)
        sun = orbits.sun(radius=0.1).move_to(centre)
        clock = always_redraw(lambda: VGroup(
            layout.label(f"t = {efolds.get_value() * tau:4.0f} Myr", font_size=22),
            layout.label(f"inclination = {i0 * np.exp(efolds.get_value()):4.1f}°", font_size=22,
                         color=P.ORANGE),
        ).arrange(DOWN, buff=0.18, aligned_edge=LEFT).move_to([4.3, 1.6, 0], aligned_edge=LEFT)
            .shift(LEFT * 1.4))
        note = VGroup(
            layout.label(f"disc: {fid['mass_earth']:.0f} M⊕ at {fid['a_au']:.0f} AU", font_size=17,
                         color=P.MUTED),
            layout.label(f"orbital period {fid['period_yr']:,.0f} yr", font_size=17, color=P.MUTED),
            layout.label(f"e-folding time {tau:.0f} Myr", font_size=17, color=P.MUTED),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([2.9, -0.2, 0], aligned_edge=LEFT)
        cap = layout.caption("A flat disc of eccentric orbits, held together only by its own gravity",
                             font_size=22)
        self.play(FadeIn(sun), FadeIn(live), FadeIn(cap), run_time=1.0)
        timing.hold_to_read(self, cap, settle=0.5)
        cap2 = layout.caption("Every orbit pitches over the same way: the disc becomes a cone",
                              font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), FadeIn(clock), FadeIn(note), run_time=0.8)
        self.play(efolds.animate.set_value(3.5), run_time=5.0, rate_func=linear)
        timing.hold_to_read(self, cap2, note, settle=0.8)
        self.play(FadeOut(VGroup(live, sun, clock, note, cap2)), run_time=0.7)

        # 2. computed: the e-folding time against disc mass
        masses = np.array(d["masses_earth"])
        ax = LogAxes((0.1, 100), (1, 1e5), [0.1, 1, 10, 100], [1, 10, 100, 1000, 1e4, 1e5],
                     "disc mass  (Earth masses)", "e-folding time  (Myr)")
        self.play(FadeIn(ax), run_time=0.9)
        shade = ax.band(m_lo, m_hi, 1, 1e5, P.ORANGE)
        shade_lab = layout.label(f"paper: {m_lo:.0f}–{m_hi:.0f} M⊕", font_size=15, color=P.ORANGE)
        shade_lab.next_to(ax.p(np.sqrt(m_lo * m_hi), 1e5), DOWN, buff=0.1)
        age = DashedLine(ax.p(0.1, 4500), ax.p(100, 4500), color=P.RED, stroke_width=2)
        age_lab = layout.label("age of the Solar System", font_size=15, color=P.RED)
        age_lab.next_to(ax.p(100, 4500), UP, buff=0.08).align_to(ax.p(100, 1), RIGHT)
        cols = [P.TEAL, P.GREEN, P.PURPLE]
        curves, tags = VGroup(), VGroup()
        for c, col in zip(d["efolding"], cols):
            t = np.array(c["tau_myr"])
            curves.add(ax.curve(masses, t, col))
            tags.add(layout.label(f"{c['a_au']:.0f} AU", font_size=15, color=col)
                     .next_to(ax.p(100, t[-1]), RIGHT, buff=0.1))
        cap3 = layout.caption("Computed: one e-fold takes the orbital period × (M☉ / disc mass)",
                              font_size=22)
        self.play(FadeIn(shade), FadeIn(shade_lab), FadeIn(cap3), run_time=0.7)
        self.play(*[Create(c) for c in curves], FadeIn(tags), run_time=1.6)
        self.play(Create(age), FadeIn(age_lab), run_time=0.6)
        dot_lo = Dot(ax.p(m_lo, d["tau_low_mass_myr"]), radius=0.07, color=P.GREEN)
        dot_hi = Dot(ax.p(m_hi, d["tau_fiducial_myr"]), radius=0.07, color=P.GREEN)
        lab_lo = layout.label(f"{d['tau_low_mass_myr'] / 1000:.1f} Gyr", font_size=15,
                              color=P.GREEN)
        lab_lo.next_to(dot_lo, LEFT, buff=0.12).shift(DOWN * 0.18)
        lab_hi = layout.label(f"{d['tau_fiducial_myr']:.0f} Myr", font_size=15, color=P.GREEN)
        lab_hi.next_to(dot_hi, LEFT, buff=0.12).shift(DOWN * 0.18)
        self.play(FadeIn(dot_lo), FadeIn(dot_hi), FadeIn(lab_lo), FadeIn(lab_hi), run_time=0.6)
        timing.hold_to_read(self, cap3, lab_lo, lab_hi, settle=1.0)
        self.play(FadeOut(VGroup(ax, shade, shade_lab, age, age_lab, curves, tags, dot_lo,
                                 dot_hi, lab_lo, lab_hi, cap3)), run_time=0.7)

        # 3. computed: the mass needed to beat the giant planets
        radii = np.array(d["radii_au"])
        ax2 = LogAxes((50, 1000), (0.01, 1000), [50, 100, 200, 500, 1000],
                      [0.01, 0.1, 1, 10, 100, 1000], "disc semi-major axis  (AU)",
                      "critical disc mass  (M⊕)", y_length=4.2, centre=(0.3, 0.6, 0))
        shade2 = ax2.band(50, 1000, m_lo, m_hi, P.ORANGE)
        shade2_lab = layout.label(f"paper: {m_lo:.0f}–{m_hi:.0f} M⊕", font_size=15, color=P.ORANGE)
        shade2_lab.next_to(ax2.p(50, m_lo), UP + RIGHT, buff=0.1)
        crit, ctags = VGroup(), VGroup()
        for c, col in zip(d["critical_mass"], cols):
            m = np.array(c["m_crit_earth"])
            crit.add(ax2.curve(radii, m, col))
            ctags.add(layout.label(f"e = {c['e']:.1f}", font_size=15, color=col)
                      .next_to(ax2.p(1000, m[-1]), RIGHT, buff=0.1))
        regions = VGroup(
            layout.label("unstable: self-gravity wins", font_size=16, color=P.GREEN)
            .move_to(ax2.p(400, 250)),
            layout.label("stable: the giant planets shear the disc apart", font_size=16,
                         color=P.RED).move_to(ax2.p(180, 0.03)),
        )
        mark = Dot(ax2.p(fid["a_au"], d["m_crit_fiducial_earth"]), radius=0.08, color=P.GREEN)
        mark_lab = layout.label(
            f"{d['m_crit_fiducial_earth']:.1f} M⊕ at {fid['a_au']:.0f} AU, e = {fid['e']:.1f}",
            font_size=15, color=P.GREEN).move_to(ax2.p(140, 300))
        leader = Line(mark_lab.get_bottom() + DOWN * 0.06, mark.get_center(), color=P.GREEN,
                      stroke_width=1.2).set_stroke(opacity=0.6)
        mark_lab = VGroup(mark_lab, leader)
        cap4 = layout.caption(
            "Computed: the disc mass whose self-gravity outruns the planets' precession",
            font_size=22)
        self.play(FadeIn(ax2), FadeIn(cap4), run_time=0.9)
        self.play(*[Create(c) for c in crit], FadeIn(ctags), run_time=1.5)
        self.play(FadeIn(regions), FadeIn(shade2), FadeIn(shade2_lab), run_time=0.7)
        self.play(FadeIn(mark), FadeIn(mark_lab), run_time=0.5)
        timing.hold_to_read(self, cap4, regions, mark_lab, settle=1.2)
        self.play(FadeOut(cap4), run_time=0.4)

        layout.show_takeaway(
            self, f"No planet needed, but {m_lo:.0f}–{m_hi:.0f} Earth masses must sit at hundreds of AU.")
