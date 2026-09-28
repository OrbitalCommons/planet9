"""Preface F -- Indirect fingerprints, and where Planet Nine could still be.

Scenes: P12Ranging, P12Obliquity, P13Map.

Every plotted number comes from ``anim.json -> preface -> f_indirect_map``
(``crates/p9-anim-data/src/preface/f_indirect_map.rs``).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    PI,
    RIGHT,
    UP,
    Arc,
    Arrow,
    Circle,
    Create,
    DashedLine,
    DashedVMobject,
    Dot,
    Ellipse,
    FadeIn,
    FadeOut,
    Line,
    Rectangle,
    Scene,
    Text,
    ValueTracker,
    VGroup,
    VMobject,
    Write,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, timing


def _data():
    return dataio.section("preface")["f_indirect_map"]


def _header(scene, badge, title):
    scene.add(layout.concept_badge(badge))
    t = Text(title, color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.55)
    scene.play(Write(t), run_time=0.9)
    return t


def _say(scene, text, old=None, hold=True, settle=0.8):
    """Swap in a caption (fading out ``old``), optionally holding to read it."""
    cap = layout.caption(text)
    if cap.width > 13.4:
        cap.scale_to_fit_width(13.4)
    anims = [FadeIn(cap, shift=UP * 0.1)]
    if old is not None:
        anims.append(FadeOut(old))
    scene.play(*anims, run_time=0.5)
    if hold:
        timing.hold_to_read(scene, cap, settle=settle)
    return cap


class _Frame:
    """A plotting frame: a linear map from data (x, y) to the scene, with the
    axes drawn as an L at the range minima (never crossing mid-plot)."""

    def __init__(self, x0, x1, y0, y1, width, height, center):
        self.x0, self.x1, self.y0, self.y1 = x0, x1, y0, y1
        self.w, self.h = width, height
        self.c = np.array([center[0], center[1], 0.0])

    def p(self, x, y):
        u = (x - self.x0) / (self.x1 - self.x0) - 0.5
        v = (y - self.y0) / (self.y1 - self.y0) - 0.5
        return self.c + np.array([u * self.w, v * self.h, 0.0])

    def axes(self, xticks=(), yticks=(), xlabel=None, ylabel=None, font=17):
        g = VGroup()
        o = self.p(self.x0, self.y0)
        g.add(Line(o, self.p(self.x1, self.y0), color=P.MUTED, stroke_width=2))
        g.add(Line(o, self.p(self.x0, self.y1), color=P.MUTED, stroke_width=2))
        xl, yl = VGroup(), VGroup()
        for x, txt in xticks:
            a = self.p(x, self.y0)
            g.add(Line(a, a + DOWN * 0.08, color=P.MUTED, stroke_width=2))
            xl.add(layout.label(txt, font_size=font, color=P.MUTED).next_to(a, DOWN, buff=0.14))
        for y, txt in yticks:
            a = self.p(self.x0, y)
            g.add(Line(a, a + LEFT * 0.08, color=P.MUTED, stroke_width=2))
            yl.add(layout.label(txt, font_size=font, color=P.MUTED).next_to(a, LEFT, buff=0.14))
        g.add(xl, yl)
        if xlabel:
            lab = layout.label(xlabel, font_size=font + 2, color=P.FG)
            top = xl.get_bottom()[1] if len(xl) else o[1]
            lab.move_to([self.c[0], top - 0.1 - lab.height / 2, 0])
            g.add(lab)
        if ylabel:
            lab = layout.label(ylabel, font_size=font + 2, color=P.FG).rotate(PI / 2)
            left = yl.get_left()[0] if len(yl) else o[0]
            lab.move_to([left - 0.12 - lab.width / 2, self.c[1], 0])
            g.add(lab)
        return g

    def curve(self, xs, ys, color, stroke_width=3.0):
        m = VMobject(color=color, stroke_width=stroke_width)
        m.set_points_as_corners([self.p(x, y) for x, y in zip(xs, ys)])
        return m


def _sci(x, digits=1):
    """'2.9×10⁻¹³' style scientific notation."""
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10 ** exp
    sup = str(exp).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"{mant:.{digits}f}×10{sup}"


def _spread(ys, gap):
    """Push sorted label heights apart so neighbours sit at least ``gap`` apart."""
    order = np.argsort(ys)
    out = np.array(ys, dtype=float)
    for k in range(1, len(order)):
        a, b = order[k - 1], order[k]
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    return out


_STATUS_COLOR = {"excluded": P.RED, "allowed": P.TEAL, "undetectable": P.MUTED}


# ---------------------------------------------------------------------------
# P12 -- Cassini ranging: a tidal tug on Saturn
# ---------------------------------------------------------------------------

class P12Ranging(Scene):
    """Goal: an unseen planet still tugs on Saturn with a tidal acceleration
    ~ G M r / d^3; Cassini's ~75 m ranging would have felt it near P9's
    perihelion, so those positions along the orbit are ruled out."""

    def construct(self):
        d = _data()["ranging"]
        _header(self, "PREFACE 12", "Fingerprints I: a tug on Saturn")

        # ---- beat 1: the tidal tug (schematic) -------------------------------
        sun_pt = np.array([-4.2, 0.1, 0])
        p9_pt = np.array([2.7, 0.1, 0])
        r_sat = 1.5
        k = 1.3 * np.linalg.norm(p9_pt - sun_pt) ** 2

        def pull(pt):
            v = p9_pt - pt
            return k * v / np.linalg.norm(v) ** 3

        phi = ValueTracker(np.deg2rad(40))

        def sat():
            return sun_pt + r_sat * np.array([np.cos(phi.get_value()), np.sin(phi.get_value()), 0])

        sun = Dot(sun_pt, radius=0.16, color=P.SUN).set_z_index(3)
        s_orbit = Circle(radius=r_sat, color=P.MUTED, stroke_width=1.5).move_to(sun_pt)
        p9 = Dot(p9_pt, radius=0.14, color=P.BLUE)
        p9_lab = layout.label("Planet Nine", font_size=19, color=P.BLUE).next_to(p9, UP, buff=0.15)
        sat_dot = always_redraw(lambda: Dot(sat(), radius=0.09, color=P.FG).set_z_index(3))
        sun_arrow = Arrow(sun_pt, sun_pt + pull(sun_pt), buff=0, color=P.MUTED, stroke_width=5,
                          max_tip_length_to_length_ratio=0.2)
        sun_lab = layout.label("grey: pull toward Planet Nine", font_size=18, color=P.MUTED)
        sun_lab.move_to(sun_pt + np.array([0, -r_sat - 0.8, 0]))

        def sat_arrows():
            s = sat()
            a_sat = pull(s)
            tide = a_sat - pull(sun_pt)
            return VGroup(
                Arrow(s, s + a_sat, buff=0, color=P.MUTED, stroke_width=4,
                      max_tip_length_to_length_ratio=0.2),
                Arrow(s, s + 2 * tide, buff=0, color=P.ORANGE, stroke_width=7,
                      max_tip_length_to_length_ratio=0.4),
            )

        arrows = always_redraw(sat_arrows)
        tide_lab = layout.label("orange: Saturn's pull minus the Sun's = the tidal tug",
                                font_size=18, color=P.ORANGE)
        tide_lab.move_to(sun_pt + np.array([0.9, r_sat + 0.45, 0]))
        not_scale = layout.label("(not to scale; tug drawn ×2)", font_size=16, color=P.MUTED)
        not_scale.move_to([p9_pt[0], -0.45, 0])
        s_lab = always_redraw(lambda: layout.label("Saturn", font_size=18, color=P.FG).next_to(
            sat(), UP if np.sin(phi.get_value()) >= 0 else DOWN, buff=0.15))

        cap = _say(self, "Planet Nine pulls on the Sun and on Saturn almost equally.", hold=False)
        self.play(FadeIn(sun), Create(s_orbit), FadeIn(p9), FadeIn(p9_lab), FadeIn(not_scale),
                  FadeIn(sat_dot), FadeIn(s_lab), run_time=1.0)
        self.play(Create(sun_arrow), FadeIn(sun_lab))
        self.play(FadeIn(arrows), FadeIn(tide_lab))
        timing.hold_to_read(self, cap, settle=0.3)
        cap = _say(self, "Only the difference changes Saturn's orbit about the Sun: it stretches "
                         "along the line to Planet Nine.", cap, hold=False)
        self.play(phi.animate.set_value(np.deg2rad(40) + 2 * np.pi), run_time=6.0,
                  rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap, settle=0.2)

        fav_nu = d["favored_nu_deg"]
        tide_fav = float(np.interp(fav_nu, d["nu_deg"], d["tidal_m_s2"]))
        real = VGroup(
            layout.label(f"real sizes, Planet Nine {d['favored_distance_au']:.0f} AU away:",
                         font_size=18, color=P.FG),
            layout.label(f"pull on the Sun {_sci(d['sun_pull_m_s2'])} m/s²", font_size=18,
                         color=P.MUTED),
            layout.label(f"tidal tug on Saturn {_sci(tide_fav)} m/s²", font_size=18,
                         color=P.ORANGE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.1)
        real.move_to([2.6, -1.9, 0])
        sig_fav = 1e3 * float(np.interp(fav_nu, d["nu_deg"], d["prefit_km"]))
        cap = _say(self, f"Tiny, yet over Cassini's decade it shifts Saturn's range by ~{sig_fav:.0f} m:"
                         f" Cassini measured to ~{1e3 * d['floor_km']:.0f} m.", cap, hold=False)
        self.play(FadeIn(real))
        timing.hold_to_read(self, cap, settle=1.0)
        self.play(FadeOut(VGroup(sun, s_orbit, p9, p9_lab, sat_dot, sun_arrow, sun_lab, arrows,
                                 tide_lab, not_scale, s_lab, real)), FadeOut(cap))

        eq = layout.explain_equation(
            self,
            [r"a_{\rm tide}", r"\;\approx\;", r"G M_9", r"\cdot", r"r_{\rm Sat}", r"\cdot",
             r"d^{-3}"],
            [
                (2, "G M₉: Planet Nine's mass sets how hard it pulls"),
                (4, f"r_Sat: Saturn's {d['saturn_au']:.1f} AU lever arm from the Sun"),
                (6, "d⁻³: twice as far away, eight times weaker"),
            ],
            scale=1.3,
            where=UP * 0.8,
        )
        self.play(FadeOut(eq))

        # ---- beat 2: walk Planet Nine around its orbit ------------------------
        nu = np.asarray(d["nu_deg"])
        dist = np.asarray(d["distance_au"])
        sig_m = np.asarray(d["prefit_km"]) * 1e3
        status = d["status"]
        a, e = d["a"], d["e"]
        s = 4.9 / (a * (1 + e) + a * (1 - e))
        xc = np.array([-4.05, 0.35, 0])
        sun2 = xc + np.array([a * e * s, 0, 0])

        def orbit_pt(nu_deg):
            t = np.deg2rad(nu_deg)
            r = a * (1 - e * e) / (1 + e * np.cos(t))
            return sun2 + r * s * np.array([np.cos(t), np.sin(t), 0])

        segs = VGroup()
        for i in range(len(nu)):
            n0, n1 = nu[i], nu[i] + 2.0
            pts = [orbit_pt(v) for v in np.linspace(n0, n1, 4)]
            m = VMobject(color=_STATUS_COLOR[status[i]], stroke_width=6)
            m.set_points_as_corners(pts)
            segs.add(m)
        orbit_lab = VGroup(
            layout.label("Batygin & Brown (2016) orbit:", font_size=18, color=P.BLUE),
            layout.label(f"{d['mass_earth']:.0f} M⊕, a = {a:.0f} AU, e = {e:.1f}", font_size=18,
                         color=P.BLUE),
        ).arrange(DOWN, buff=0.06).move_to([xc[0], 2.62, 0])
        sun_d = Dot(sun2, radius=0.09, color=P.SUN)
        sun_l = layout.label("Sun", font_size=16, color=P.SUN).next_to(sun_d, DOWN, buff=0.1)
        peri_l = layout.label("perihelion", font_size=16, color=P.FG)
        peri_l.next_to(orbit_pt(0), RIGHT, buff=0.12)

        fr = _Frame(0, 360, 1.0, np.log10(3000), 5.4, 3.9, (3.95, 0.55))
        axes = fr.axes(
            xticks=[(v, f"{v}°") for v in (0, 90, 180, 270, 360)],
            yticks=[(np.log10(v), f"{v}") for v in (10, 30, 100, 300, 1000, 3000)],
            xlabel="where Planet Nine is on its orbit",
            ylabel="range signal (m)",
        )
        bands = VGroup()
        for i in range(len(nu)):
            a0, b0 = fr.p(nu[i], 1.0), fr.p(nu[i] + 2.0, np.log10(3000))
            bands.add(Rectangle(width=b0[0] - a0[0], height=b0[1] - a0[1], stroke_width=0)
                      .set_fill(_STATUS_COLOR[status[i]], opacity=0.13).move_to((a0 + b0) / 2))
        sig_curve = fr.curve(nu, np.log10(sig_m), P.ORANGE, stroke_width=3.5)
        floor_y = np.log10(d["floor_km"] * 1e3)
        floor = DashedLine(fr.p(0, floor_y), fr.p(360, floor_y), color=P.FG, stroke_width=2)
        floor_l = layout.label(f"Cassini precision ≈ {d['floor_km'] * 1e3:.0f} m", font_size=17,
                               color=P.FG)
        floor_l.next_to(fr.p(180, floor_y), DOWN, buff=0.1)

        keys = VGroup()
        for col, text in ((P.RED, "would have shown: ruled out"), (P.TEAL, "allowed"),
                          (P.MUTED, "too faint to tell")):
            sw = Rectangle(width=0.28, height=0.16, stroke_width=0).set_fill(col, opacity=0.9)
            keys.add(VGroup(sw, layout.label(text, font_size=16, color=P.FG)).arrange(RIGHT,
                                                                                    buff=0.12))
        keys.arrange(RIGHT, buff=0.35).move_to([0.0, -2.45, 0])

        cap = _say(self, "Now move Planet Nine around its orbit: would Cassini have noticed?",
                   hold=False)
        self.play(FadeIn(orbit_lab), FadeIn(sun_d), FadeIn(sun_l), FadeIn(axes), run_time=0.8)
        orbit_grey = VGroup(*[m.copy().set_color(P.MUTED).set_stroke(opacity=0.5) for m in segs])
        self.play(Create(orbit_grey), FadeIn(peri_l), Create(floor), FadeIn(floor_l), run_time=1.0)

        t = ValueTracker(0.0)

        def here():
            v = t.get_value()
            return v, float(np.interp(v, nu, dist)), float(np.interp(v, nu, np.log10(sig_m)))

        p9_dot = always_redraw(lambda: Dot(orbit_pt(t.get_value()), radius=0.11,
                                           color=P.BLUE).set_z_index(5))
        tracer = always_redraw(lambda: Dot(fr.p(here()[0], here()[2]), radius=0.08,
                                           color=P.ORANGE).set_z_index(5))

        def readout():
            v, dd, ls = here()
            return layout.label(f"{dd:,.0f} AU away:  signal {10 ** ls:,.0f} m", font_size=18,
                                color=P.FG).move_to([xc[0], -1.95, 0])

        curve_drawn = always_redraw(
            lambda: fr.curve(nu[nu <= t.get_value() + 1e-9],
                             np.log10(sig_m[nu <= t.get_value() + 1e-9]), P.ORANGE,
                             stroke_width=3.5) if t.get_value() > 2 else VGroup())
        read_m = always_redraw(readout)
        self.add(p9_dot, tracer, curve_drawn, read_m)

        def reveal(upto):
            idx = [i for i in range(len(nu)) if nu[i] < upto and i not in reveal.done]
            reveal.done.update(idx)
            return [FadeIn(segs[i]) for i in idx] + [FadeIn(bands[i]) for i in idx]

        reveal.done = set()
        for stop, text in ((110, "Near perihelion the tug would have spoiled Saturn's orbit fit: "
                                 "ruled out."),
                           (212, "Far out it is too weak to see either way; between the two, "
                                 "windows survive."),
                           (360, None)):
            self.play(t.animate.set_value(stop), *reveal(stop), run_time=(stop - t.get_value()) / 40,
                      rate_func=rate_functions.linear)
            if text:
                cap = _say(self, text, cap, settle=0.4)
        self.remove(curve_drawn)
        self.add(sig_curve)
        self.play(FadeIn(keys), run_time=0.6)

        fav_line = DashedLine(fr.p(fav_nu, 1.0), fr.p(fav_nu, np.log10(3000)), color=P.BLUE,
                              stroke_width=2.5)
        lo, hi = d["favored_interval_deg"]
        y_bar = np.log10(13)
        paper = VGroup(Line(fr.p(lo, y_bar), fr.p(hi, y_bar), color=P.FG, stroke_width=5),
                       Line(fr.p(lo, y_bar) + DOWN * 0.08, fr.p(lo, y_bar) + UP * 0.08,
                            color=P.FG, stroke_width=2),
                       Line(fr.p(hi, y_bar) + DOWN * 0.08, fr.p(hi, y_bar) + UP * 0.08,
                            color=P.FG, stroke_width=2))
        paper_l = layout.label(f"paper: {lo:.0f}°–{hi:.0f}°", font_size=16, color=P.FG)
        paper_l.next_to(paper, LEFT, buff=0.1)
        fav_l = layout.label(f"here: {fav_nu:.0f}°", font_size=16, color=P.BLUE)
        fav_l.next_to(fr.p(fav_nu, np.log10(700)), RIGHT, buff=0.08)
        cap = _say(self, f"The best fit sits in the window: this reproduction favours {fav_nu:.0f}°, "
                         f"{d['favored_distance_au']:.0f} AU from the Sun.", cap, hold=False)
        self.play(t.animate.set_value(fav_nu), run_time=1.5)
        self.play(Create(fav_line), FadeIn(fav_l), Create(paper), FadeIn(paper_l))
        timing.hold_to_read(self, cap, settle=1.0)
        self.play(FadeOut(cap), FadeOut(keys))
        layout.show_takeaway(self, "Cassini's ranging rules out Planet Nine near perihelion.")


# ---------------------------------------------------------------------------
# P12b -- the Sun's 6 degree tilt
# ---------------------------------------------------------------------------

class P12Obliquity(Scene):
    """Goal: the Sun's spin is tilted 6 degrees from the planets' plane; an
    inclined Planet Nine slowly twists the planets' plane over 4.5 Gyr while the
    Sun's spin stays put, which naturally produces that tilt."""

    def construct(self):
        d = _data()["obliquity"]
        _header(self, "PREFACE 12b", "Fingerprints II: the Sun's tilted spin")

        # ---- beat 1: what the obliquity is -----------------------------------
        c = np.array([0.0, -0.2, 0])
        disk = Ellipse(width=9.0, height=1.2, color=P.MUTED, stroke_width=2).move_to(c)
        planets = VGroup(*[Dot(c + np.array([4.5 * np.cos(t), 0.6 * np.sin(t), 0]), radius=0.06,
                               color=P.FG) for t in (0.4, 1.9, 3.3, 4.9)])
        sun = Circle(radius=0.55, color=P.SUN, stroke_width=0).set_fill(P.SUN, opacity=1.0)
        sun.move_to(c).set_z_index(3)
        tilt = np.deg2rad(d["observed_deg"])
        normal = Arrow(c + DOWN * 2.3, c + UP * 2.6, buff=0, color=P.FG, stroke_width=4)
        spin_dir = np.array([np.sin(tilt), np.cos(tilt), 0])
        spin = Arrow(c - 2.3 * spin_dir, c + 2.6 * spin_dir, buff=0, color=P.SUN,
                     stroke_width=4).set_z_index(4)
        arc = Arc(radius=2.2, start_angle=PI / 2 - tilt, angle=tilt, arc_center=c, color=P.ORANGE,
                  stroke_width=3)
        arc_l = layout.label(f"{d['observed_deg']:.0f}°", font_size=22, color=P.ORANGE)
        arc_l.next_to(arc, UP, buff=0.08).shift(RIGHT * 0.1)
        n_lab = layout.label("axis of the planets' orbits", font_size=18, color=P.FG)
        n_lab.next_to(normal.get_end(), LEFT, buff=0.15)
        s_lab = layout.label("the Sun's spin axis", font_size=18, color=P.SUN)
        s_lab.next_to(spin.get_end(), RIGHT, buff=0.15)
        cap = _say(self, "The planets orbit in nearly one plane, formed from the Sun's own "
                         "spinning disk.", hold=False)
        self.play(Create(disk), FadeIn(planets), FadeIn(sun), run_time=1.0)
        self.play(Create(normal), FadeIn(n_lab))
        timing.hold_to_read(self, cap, settle=0.3)
        cap = _say(self, "Yet the Sun spins about an axis tilted 6° from theirs. Why?", cap,
                   hold=False)
        self.play(Create(spin), FadeIn(s_lab))
        self.play(Create(arc), FadeIn(arc_l))
        timing.hold_to_read(self, cap, settle=0.8)
        self.play(FadeOut(VGroup(disk, planets, sun, normal, spin, arc, arc_l, n_lab, s_lab)),
                  FadeOut(cap))

        # ---- beat 2: 4.5 Gyr of secular spin-orbit evolution ---------------
        snaps = d["snapshots"]
        tg = np.array([s["t_gyr"] for s in snaps])
        obl = np.array([s["obliquity_deg"] for s in snaps])

        def deg_xy(v):
            v = np.asarray(v, dtype=float)
            r = np.linalg.norm(v)
            if r < 1e-12:
                return np.zeros(2)
            return np.degrees(np.arcsin(min(r, 1.0))) * v / r

        p9_xy = np.array([deg_xy(s["p9"]) for s in snaps])
        pl_xy = np.array([deg_xy(s["planets"]) for s in snaps])
        sun_xy = np.array([deg_xy(s["sun"]) for s in snaps])

        cc = np.array([-3.5, -0.05, 0])
        kdeg = 2.35 / 5.6

        def at(xy):
            return cc + kdeg * np.array([xy[0], xy[1], 0])

        rings = VGroup()
        for rdeg in (2, 4):
            ring = DashedVMobject(Circle(radius=kdeg * rdeg, color=P.MUTED, stroke_width=1.2),
                                  num_dashes=40).move_to(cc)
            lab = layout.label(f"{rdeg}°", font_size=15, color=P.MUTED)
            lab.move_to(cc + kdeg * rdeg * np.array([np.cos(2.2), np.sin(2.2), 0])
                        + np.array([0.12, 0.12, 0]))
            rings.add(ring, lab)
        rim = Circle(radius=kdeg * 5.6, color=P.MUTED, stroke_width=1.5).move_to(cc)
        centre_mark = VGroup(Line(cc + LEFT * 0.1, cc + RIGHT * 0.1, color=P.MUTED),
                             Line(cc + UP * 0.1, cc + DOWN * 0.1, color=P.MUTED))
        centre_lab = VGroup(layout.label("total spin of the", font_size=15, color=P.MUTED),
                            layout.label("solar system", font_size=15, color=P.MUTED)
                            ).arrange(DOWN, buff=0.04).next_to(cc, DOWN, buff=0.18)

        tt = ValueTracker(0.0)

        def interp(arr):
            x = tt.get_value()
            return np.array([np.interp(x, tg, arr[:, 0]), np.interp(x, tg, arr[:, 1])])

        def tips():
            pl, sn, p9 = at(interp(pl_xy)), at(interp(sun_xy)), interp(p9_xy)
            u = np.array([p9[0], p9[1], 0]) / max(np.linalg.norm(p9), 1e-9)
            g = VGroup(
                Line(sn, pl, color=P.ORANGE, stroke_width=3),
                Dot(pl, radius=0.09, color=P.FG).set_z_index(3),
                Dot(sn, radius=0.09, color=P.SUN).set_z_index(3),
                Arrow(cc + u * (kdeg * 5.6 + 0.05), cc + u * (kdeg * 5.6 + 0.6), buff=0,
                      color=P.BLUE, stroke_width=4),
            )
            return g

        def trail():
            x = tt.get_value()
            keep = tg <= x + 1e-9
            pts = [at(p) for p in pl_xy[keep]]
            if len(pts) < 2:
                return VGroup()
            m = VMobject(color=P.FG, stroke_width=2).set_points_as_corners(pts)
            return m.set_stroke(opacity=0.5)

        tip_m = always_redraw(tips)
        trail_m = always_redraw(trail)
        keys = VGroup(
            VGroup(Dot(radius=0.08, color=P.FG),
                   layout.label("planets' orbit axis", font_size=17, color=P.FG)).arrange(RIGHT,
                                                                                         buff=0.12),
            VGroup(Dot(radius=0.08, color=P.SUN),
                   layout.label("Sun's spin axis", font_size=17, color=P.SUN)).arrange(RIGHT,
                                                                                      buff=0.12),
            VGroup(Arrow(LEFT * 0.2, RIGHT * 0.2, buff=0, color=P.BLUE, stroke_width=4),
                   layout.label(f"toward Planet Nine's axis ({d['i9_required_deg']:.0f}° off)",
                                font_size=17, color=P.BLUE)).arrange(RIGHT, buff=0.12),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        keys.move_to([3.9, 2.35, 0])

        fr = _Frame(0.0, 4.5, 0.0, 7.0, 5.0, 2.7, (4.0, -0.5))
        axes = fr.axes(
            xticks=[(v, str(v)) for v in (0, 1, 2, 3, 4)],
            yticks=[(v, f"{v}°") for v in (0, 2, 4, 6)],
            xlabel="billions of years since the planets formed",
            ylabel="Sun's tilt",
        )
        obs = DashedLine(fr.p(0, d["observed_deg"]), fr.p(4.5, d["observed_deg"]), color=P.GREEN,
                         stroke_width=2)
        obs_l = layout.label(f"observed {d['observed_deg']:.0f}°", font_size=16, color=P.GREEN)
        obs_l.next_to(fr.p(0.1, d["observed_deg"]), UP, buff=0.08, aligned_edge=LEFT)
        obl_curve = always_redraw(
            lambda: fr.curve(tg[tg <= tt.get_value() + 1e-9], obl[tg <= tt.get_value() + 1e-9],
                             P.ORANGE, stroke_width=3.5) if tt.get_value() > 0.05 else VGroup())
        obl_dot = always_redraw(lambda: Dot(fr.p(tt.get_value(), np.interp(tt.get_value(), tg,
                                                                           obl)),
                                            radius=0.07, color=P.ORANGE))
        read = always_redraw(lambda: layout.label(
            f"{tt.get_value():.1f} billion years:  tilt {np.interp(tt.get_value(), tg, obl):.1f}°",
            font_size=20, color=P.ORANGE).move_to(fr.p(2.25, 7.0) + UP * 0.45))

        cap = _say(self, "Look straight down the solar system's total spin: each dot marks where "
                         "an axis points.", hold=False)
        self.play(Create(rim), FadeIn(rings), FadeIn(centre_mark), FadeIn(centre_lab),
                  FadeIn(keys), run_time=1.0)
        self.play(FadeIn(tip_m), FadeIn(axes), Create(obs), FadeIn(obs_l))
        self.add(trail_m, obl_curve, obl_dot, read)
        timing.hold_to_read(self, cap, settle=0.6)
        cap = _say(self, "An inclined Planet Nine slowly drags the planets' axis around; the Sun's "
                         "spin barely moves.", cap, hold=False)
        self.play(tt.animate.set_value(tg[-1]), run_time=9.0, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap, settle=0.4)
        cap = _say(self, f"{d['mass_earth']:.0f} M⊕ at {d['a']:.0f} AU, tilted "
                         f"{d['i9_required_deg']:.0f}°, opens a {obl[-1]:.1f}° gap in 4.5 billion "
                         "years (Bailey et al. 2016 model).", cap)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "The Sun's 6° tilt may itself be Planet Nine's fingerprint.")


# ---------------------------------------------------------------------------
# P13 -- where could it still be?
# ---------------------------------------------------------------------------

class P13Map(Scene):
    """Goal: the clustering predicts a cloud of (mass, distance) possibilities;
    each survey removes the part of the cloud it was deep enough (and pointed
    right) to see; the rest is within Rubin's reach."""

    def construct(self):
        d = _data()["map"]
        _header(self, "PREFACE 13", "Where could it still be?")

        y_hi = 1300.0
        fr = _Frame(2.0, 14.0, 150.0, y_hi, 8.4, 4.7, (-1.5, 0.5))
        axes = fr.axes(
            xticks=[(v, str(v)) for v in (2, 4, 6, 8, 10, 12, 14)],
            yticks=[(v, str(v)) for v in (200, 400, 600, 800, 1000, 1200)],
            xlabel="mass (Earth masses)",
            ylabel="distance from the Sun today (AU)",
        )
        draws = [o for o in d["draws"] if 2.0 <= o["mass"] <= 14.0 and o["r_au"] <= y_hi]
        dots = VGroup(*[Dot(fr.p(o["mass"], max(o["r_au"], 150.0)), radius=0.03, color=P.TEAL)
                        .set_opacity(0.85) for o in draws])
        cap = _say(self, "Each dot is one Planet Nine the clustering allows: its mass and where it "
                         "is today.", hold=False)
        self.play(FadeIn(axes), run_time=0.8)
        self.play(FadeIn(dots, lag_ratio=0.002), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.4)

        mass = np.asarray(d["mass_grid"])

        def reach_curve(key, color, dashed=False, dotted=False):
            r = np.array([np.nan if v is None else v for v in d["reach"][key]])
            keep = (mass >= 2.0) & (mass <= 14.0) & np.isfinite(r)
            xs, ys = mass[keep], r[keep]
            inside = ys <= y_hi
            if not inside.all():
                # stop where the curve leaves the top of the chart
                k = int(np.argmax(~inside))
                frac = (y_hi - ys[k - 1]) / (ys[k] - ys[k - 1])
                xe = xs[k - 1] + frac * (xs[k] - xs[k - 1])
                xs = np.append(xs[:k], xe)
                ys = np.append(ys[:k], y_hi)
            m = fr.curve(xs, ys, color, stroke_width=3)
            if dashed or dotted:
                m = DashedVMobject(m, num_dashes=60 if dashed else 110,
                                   dashed_ratio=0.55 if dashed else 0.3)
            return m, np.array([xs[-1], ys[-1]])

        # counter (right column)
        counter_head = layout.label("ruled out so far", font_size=19, color=P.FG)
        counter_head.move_to([5.1, 2.55, 0])
        counter = layout.label("0%", font_size=40, color=P.RED, weight="BOLD")
        counter.next_to(counter_head, DOWN, buff=0.15)

        def set_count(v):
            new = layout.label(f"{100 * v:.0f}%", font_size=40, color=P.RED, weight="BOLD")
            return counter.animate.become(new.move_to(counter))

        legend = VGroup(
            VGroup(Dot(radius=0.06, color=P.TEAL),
                   layout.label("still allowed", font_size=17, color=P.FG)).arrange(RIGHT,
                                                                                    buff=0.12),
            VGroup(Dot(radius=0.06, color=P.RED),
                   layout.label("would have been seen", font_size=17, color=P.FG)).arrange(
                RIGHT, buff=0.12),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        legend.move_to([5.1, -1.55, 0])
        self.play(FadeIn(counter_head), FadeIn(counter), FadeIn(legend), run_time=0.6)

        right_x = fr.p(14.0, 0)[0] + 0.15
        # right-hand labels for curves that end inside the chart, pushed apart
        ends = {key: reach_curve(key, P.PURPLE)[1] for key in ("ZTF", "Pan-STARRS", "infrared")}
        end_keys = [k for k in ends if ends[k][1] < y_hi - 1]
        end_y = dict(zip(end_keys, _spread([fr.p(0, ends[k][1])[1] for k in end_keys], 0.34)))
        steps = [
            ("ZTF", "ZTF (2021)", 0, "ZTF (2021): any dot bright enough and inside its sky "
                                     "would have been seen."),
            ("DES", "DES (2022)", 1, "DES (2022) looks deeper, but at only an eighth of the sky."),
            ("Pan-STARRS", "Pan-STARRS (2024)", 2, "Pan-STARRS (2024): deep and wide. Dots below a "
                                                   "line survive only if they hid."),
        ]
        for key, name, k, text in steps:
            curve, end = reach_curve(key, P.PURPLE)
            lab = layout.label(name, font_size=17, color=P.PURPLE)
            if key in end_y:
                lab.move_to([right_x, end_y[key], 0], aligned_edge=LEFT)
            else:
                lab.next_to(fr.p(*end), DOWN + RIGHT, buff=0.08)
            gone = [i for i, o in enumerate(draws) if o["u"] < o["p_cum"][k]]
            cap = _say(self, text, cap, hold=False)
            self.play(Create(curve), FadeIn(lab), run_time=1.2)
            self.play(*[dots[i].animate.set_color(P.RED).set_opacity(0.45) for i in gone],
                      set_count(d["excluded_frac"][k]), run_time=1.2)
            timing.hold_to_read(self, cap, settle=0.5)

        ir, _ = reach_curve("infrared", P.PURPLE, dotted=True)
        ir.set_stroke(width=4)
        ir_lab = layout.label("far-infrared", font_size=17, color=P.PURPLE)
        ir_lab.move_to([right_x, end_y["infrared"], 0], aligned_edge=LEFT)
        cap = _say(self, "Its own heat, seen by IRAS, AKARI and WISE, reaches about as far as ZTF.",
                   cap, hold=False)
        self.play(Create(ir), FadeIn(ir_lab), run_time=1.2)
        timing.hold_to_read(self, cap, settle=0.4)

        rubin, r_end = reach_curve("Rubin", P.TEAL, dashed=True)
        r_lab = layout.label("Rubin, one visit (forecast)", font_size=17, color=P.TEAL)
        r_lab.next_to(fr.p(*r_end), DOWN + LEFT, buff=0.08)
        cap = _say(self, f"Rubin (from 2025): {100 * d['rubin_bright_frac']:.1f}% of the "
                         "survivors are bright enough for a single visit.", cap, hold=False)
        self.play(Create(rubin), FadeIn(r_lab), run_time=1.4)
        timing.hold_to_read(self, cap, settle=1.2)
        cap = _say(self, "Each paper that follows moves one of these lines, or reshapes the cloud "
                         "itself.", cap, settle=1.2)
        self.play(FadeOut(cap))
        layout.show_takeaway(
            self, f"~{100 * d['excluded_frac'][-1]:.0f}% ruled out; Rubin can reach the rest.")
