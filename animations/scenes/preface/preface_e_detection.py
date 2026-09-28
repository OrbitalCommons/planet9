"""Preface E -- Seeing the unseen: reflected light, thermal glow, finding movers,
and the selection bias every survey carries.

Scenes: P09ReflectedLight, P10Thermal, P11Movers, P11SelectionBias.

Every plotted number comes from ``anim.json -> preface -> e_detection``
(``crates/p9-anim-data/src/preface/e_detection.rs``).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    PI,
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
    Line,
    MathTex,
    Rectangle,
    Scene,
    Text,
    Transform,
    ValueTracker,
    VGroup,
    VMobject,
    Write,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, timing, widgets


def _data():
    return dataio.section("preface")["e_detection"]


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


def _spread(ys, gap):
    """Push sorted label heights apart so neighbours sit at least ``gap`` apart."""
    order = np.argsort(ys)
    out = np.array(ys, dtype=float)
    for k in range(1, len(order)):
        a, b = order[k - 1], order[k]
        if out[b] - out[a] < gap:
            out[b] = out[a] + gap
    return out


# ---------------------------------------------------------------------------
# P09 -- reflected light falls as 1/r^4
# ---------------------------------------------------------------------------

class P09ReflectedLight(Scene):
    """Goal: reflected sunlight dims as 1/r^4, so doubling the distance costs
    3 magnitudes, and each survey's depth sets how far out it could see."""

    def construct(self):
        d = _data()["reflected"]
        _header(self, "PREFACE 09", "Seeing it I: borrowed sunlight")

        eq = layout.explain_equation(
            self,
            [r"F", r"\;\propto\;", r"p", r"\,R^2", r"\times", r"\frac{1}{r^2}", r"\times",
             r"\frac{1}{\Delta^2}"],
            [
                (2, "p, the albedo: the fraction of sunlight it reflects"),
                (3, "R, its radius: a bigger disk catches more light"),
                (5, "sunlight thins out on the way to it: 1/r²"),
                (7, "and thins again on the way back to us: 1/Δ², with Δ ≈ r"),
            ],
            scale=1.3,
            where=UP * 1.0,
        )
        four = MathTex(r"\Longrightarrow\quad F \propto \frac{1}{r^4}", color=P.TEAL).scale(1.1)
        four.next_to(eq, DOWN, buff=0.6)
        self.play(Write(four))
        cap = _say(self, "Twice as far means 2⁴ = 16 times fainter.")
        self.play(FadeOut(eq), FadeOut(four), FadeOut(cap))

        # ---- magnitude vs distance, with a slider ---------------------------
        dist = np.asarray(d["distance_au"])
        vmag = np.asarray(d["v_mag"])
        sun = np.asarray(d["sunlight_vs_earth"])
        lx = np.log10(dist)
        fr = _Frame(np.log10(20), np.log10(1500), -26.5, -5.5, 8.3, 4.6, (-1.3, 0.55))
        axes = fr.axes(
            xticks=[(np.log10(v), str(v)) for v in (20, 50, 100, 200, 500, 1000)],
            yticks=[(-m, str(m)) for m in (8, 12, 16, 20, 24)],
            xlabel="distance from the Sun (AU, log scale)",
            ylabel="V magnitude (fainter ↓)",
        )
        self.play(FadeIn(axes), run_time=0.8)
        cap = _say(self, "Astronomers count brightness in magnitudes: +5 mag is 100 times fainter.")

        nep, plu = d["neptune"], d["pluto"]
        anchors = VGroup()
        for body, name in ((nep, "Neptune"), (plu, "Pluto today")):
            pt = fr.p(np.log10(body["r_au"]), -body["v"])
            dot = Dot(pt, radius=0.07, color=P.GREEN)
            lab = layout.label(f"{name}  V {body['v']:.1f}", font_size=18, color=P.GREEN)
            if name == "Neptune":
                lab.next_to(dot, RIGHT, buff=0.12)
            else:
                lab.next_to(dot, DOWN, buff=0.1).align_to(dot, LEFT)
            anchors.add(VGroup(dot, lab))
        self.play(FadeIn(anchors, lag_ratio=0.4), run_time=1.0)

        curve = fr.curve(lx, -vmag, P.BLUE, stroke_width=3.5)
        model = VGroup(
            layout.label(f"Planet Nine model: {d['mass_earth']:.1f} M⊕,", font_size=18,
                         color=P.BLUE),
            layout.label(f"{d['radius_earth']:.1f} R⊕, albedo {d['albedo']:.2f}", font_size=18,
                         color=P.BLUE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        model.move_to(fr.p(np.log10(22), -18.3), aligned_edge=LEFT)
        cap = _say(self, "Now put a Planet Nine at Neptune's distance and slide it outward.", cap,
                   hold=False)
        self.play(Create(curve), FadeIn(model), run_time=1.6)

        r = ValueTracker(30.0)

        def v_at(x):
            return float(np.interp(np.log10(x), lx, vmag))

        def readout():
            x = r.get_value()
            s = float(np.interp(np.log10(x), lx, sun))
            rows = VGroup(
                layout.label(f"r = {x:,.0f} AU", font_size=21, color=P.BLUE),
                layout.label(f"sunlight there: 1/{1 / s:,.0f} of Earth's", font_size=19,
                             color=P.FG),
                layout.label(f"brightness V = {v_at(x):.1f}", font_size=19, color=P.FG),
            ).arrange(DOWN, aligned_edge=RIGHT, buff=0.12)
            rows.move_to(fr.p(np.log10(1450), -6.0), aligned_edge=UP + RIGHT)
            return rows

        probe = always_redraw(lambda: Dot(fr.p(np.log10(r.get_value()), -v_at(r.get_value())),
                                          radius=0.09, color=P.BLUE).set_z_index(4))
        panel = always_redraw(readout)
        self.add(probe)
        self.play(FadeIn(panel), run_time=0.4)
        self.play(r.animate.set_value(300.0), run_time=3.0, rate_func=rate_functions.smooth)
        mark300 = Dot(probe.get_center(), radius=0.06, color=P.MUTED).set_z_index(3)
        self.add(mark300)
        self.play(r.animate.set_value(600.0), run_time=2.0, rate_func=rate_functions.smooth)

        v300, v600 = d["v_300"], d["v_600"]
        p300 = fr.p(np.log10(300), -v300)
        corner = fr.p(np.log10(600), -v300)
        p600 = fr.p(np.log10(600), -v600)
        bracket = VGroup(
            DashedLine(p300, corner, color=P.FG, stroke_width=1.6),
            Arrow(corner, p600, buff=0.05, color=P.FG, stroke_width=2.5,
                  max_tip_length_to_length_ratio=0.18),
        )
        blab = VGroup(
            layout.label("2× farther:", font_size=18, color=P.FG),
            layout.label(f"+{v600 - v300:.1f} mag = 16× fainter", font_size=18, color=P.FG),
        ).arrange(DOWN, buff=0.08)
        blab.move_to(fr.p(np.log10(600), -(v300 - 3.0)))
        self.play(Create(bracket), FadeIn(blab))
        cap = _say(self, "Every doubling of distance costs another 3 magnitudes.", cap)

        # ---- survey depths -------------------------------------------------
        surveys = d["surveys"]
        heights = _spread([fr.p(0, -s["v_limit"])[1] for s in surveys], 0.34)
        depth_g = VGroup()
        for s, y in zip(surveys, heights):
            yv = -s["v_limit"]
            line = DashedLine(fr.p(fr.x0, yv), fr.p(fr.x1, yv), color=P.PURPLE, stroke_width=1.6,
                              dash_length=0.1).set_stroke(opacity=0.8)
            hit = Dot(fr.p(np.log10(s["reach_au"]), yv), radius=0.06, color=P.PURPLE).set_z_index(3)
            lab = layout.label(
                f"{s['name']} {s['band']}<{s['depth']:g}: {round(s['reach_au'], -1):.0f} AU",
                font_size=17, color=P.PURPLE)
            lab.move_to([fr.p(fr.x1, 0)[0] + 0.15, y, 0], aligned_edge=LEFT)
            depth_g.add(VGroup(line, hit, lab))
        cap = _say(self, "A survey sees it only out to where the curve crosses its depth limit.",
                   cap, hold=False)
        self.play(FadeOut(blab), FadeOut(bracket),
                  FadeIn(depth_g, lag_ratio=0.35), run_time=2.4)
        timing.hold_to_read(self, cap, settle=1.5)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, f"Reflected light fades as 1/r⁴: V ≈ {v600:.1f} at 600 AU.")


# ---------------------------------------------------------------------------
# P10 -- its own heat: ~40 K, peaking in the far infrared, fading as 1/r^2
# ---------------------------------------------------------------------------

class P10Thermal(Scene):
    """Goal: far from the Sun a giant planet is kept at ~40 K by its own
    internal heat; that glow peaks near 100 um and fades only as 1/r^2."""

    def construct(self):
        d = _data()["thermal"]
        _header(self, "PREFACE 10", "Seeing it II: its own faint heat")

        # ---- beat 1: what sets the temperature -----------------------------
        dist = np.asarray(d["distance_au"])
        teq = np.asarray(d["t_eq"])
        t_int = d["t_internal"]
        fr = _Frame(1.0, np.log10(1500), np.log10(5), np.log10(120), 8.6, 4.5, (-0.6, 0.55))
        axes = fr.axes(
            xticks=[(np.log10(v), str(v)) for v in (10, 30, 100, 300, 1000)],
            yticks=[(np.log10(v), f"{v}") for v in (5, 10, 20, 40, 80)],
            xlabel="distance from the Sun (AU, log scale)",
            ylabel="temperature (K)",
        )
        self.play(FadeIn(axes), run_time=0.8)
        eq_curve = fr.curve(np.log10(dist), np.log10(teq), P.SUN, stroke_width=3)
        eq_lab = layout.label("warmed by sunlight alone: T ∝ 1/√r", font_size=18, color=P.SUN)
        eq_lab.move_to(fr.p(np.log10(11), np.log10(105)), aligned_edge=LEFT)
        cap = _say(self, "Sunlight alone would leave a distant planet frigid.", hold=False)
        self.play(Create(eq_curve), FadeIn(eq_lab), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.4)

        t600 = d["t_eq_600"]
        p600 = Dot(fr.p(np.log10(600), np.log10(t600)), radius=0.07, color=P.SUN)
        l600 = layout.label(f"{t600:.0f} K at 600 AU", font_size=18, color=P.SUN)
        l600.next_to(p600, DOWN + LEFT, buff=0.1)
        self.play(FadeIn(p600), FadeIn(l600))

        nep = d["neptune"]
        n_obs = Dot(fr.p(np.log10(nep["r_au"]), np.log10(nep["t_obs"])), radius=0.07, color=P.GREEN)
        n_eq = Circle(radius=0.07, color=P.SUN, stroke_width=2).move_to(
            fr.p(np.log10(nep["r_au"]), np.log10(nep["t_eq"])))
        n_lab = VGroup(
            layout.label(f"Neptune: measured {nep['t_obs']:.0f} K,", font_size=18, color=P.GREEN),
            layout.label(f"sunlight alone gives {nep['t_eq']:.0f} K", font_size=18, color=P.GREEN),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        n_lab.move_to(fr.p(np.log10(48), np.log10(68)), aligned_edge=LEFT)
        cap = _say(self, "But giant planets still leak the heat of their formation.", cap,
                   hold=False)
        self.play(FadeIn(n_obs), Create(n_eq), FadeIn(n_lab))
        timing.hold_to_read(self, cap, settle=0.6)

        floor = DashedLine(fr.p(fr.x0, np.log10(t_int)), fr.p(fr.x1, np.log10(t_int)),
                           color=P.BLUE, stroke_width=2)
        teff = fr.curve(np.log10(dist), np.log10(np.maximum(teq, t_int)), P.BLUE, stroke_width=6)
        teff.set_stroke(opacity=0.45)
        f_lab = layout.label(f"internal heat ≈ {t_int:.0f} K", font_size=18, color=P.BLUE)
        f_lab.next_to(fr.p(np.log10(700), np.log10(t_int)), UP, buff=0.12)
        cap = _say(self, f"Beyond ~{d['crossover_au']:.0f} AU the planet's own heat wins: "
                         f"Planet Nine should sit near {t_int:.0f} K.", cap, hold=False)
        self.play(Create(floor), FadeIn(f_lab))
        self.play(Create(teff), run_time=1.4)
        timing.hold_to_read(self, cap, settle=1.0)
        beat1 = VGroup(axes, eq_curve, eq_lab, p600, l600, n_obs, n_eq, n_lab, floor, teff, f_lab)
        self.play(FadeOut(beat1), FadeOut(cap))

        # ---- beat 2: cool the Planck spectrum from 5778 K to 40 K ------------
        wl = np.asarray(d["wavelength_um"])
        temps = np.asarray(d["temps"])
        spectra = np.asarray(d["spectra"])
        peaks = np.asarray(d["peak_um"])
        lw = np.log10(wl)
        fr2 = _Frame(-1.0, np.log10(3000), 0.0, 1.15, 9.4, 3.9, (-0.1, 0.0))
        axes2 = fr2.axes(
            xticks=[(np.log10(v), f"{v:g}") for v in (0.1, 1, 10, 100, 1000)],
            yticks=[(0.0, "0"), (0.5, "0.5"), (1.0, "1")],
            xlabel="wavelength (µm, log scale)",
            ylabel="glow per unit frequency (peak = 1)",
        )
        bands = VGroup()
        for k, b in enumerate(d["bands"]):
            x = np.log10(b["um"])
            ln = DashedLine(fr2.p(x, 0), fr2.p(x, 1.08), color=P.PURPLE, stroke_width=1.4,
                            dash_length=0.08).set_stroke(opacity=0.7)
            lab = layout.label(b["name"], font_size=17, color=P.PURPLE)
            top = fr2.p(x, 1.08)
            if b["name"] == "IRAS":
                lab.next_to(top, UP + LEFT, buff=0.04)
            elif b["name"] == "AKARI":
                lab.next_to(top, UP + RIGHT, buff=0.04)
            else:
                lab.next_to(top, UP, buff=0.06)
            bands.add(VGroup(ln, lab))
        self.play(FadeIn(axes2), FadeIn(bands, lag_ratio=0.2), run_time=1.0)

        idx = ValueTracker(0.0)

        def state():
            f = idx.get_value()
            i = int(min(np.floor(f), len(temps) - 2))
            w = f - i
            spec = (1 - w) * spectra[i] + w * spectra[i + 1]
            t = float(np.exp((1 - w) * np.log(temps[i]) + w * np.log(temps[i + 1])))
            pk = float(np.exp((1 - w) * np.log(peaks[i]) + w * np.log(peaks[i + 1])))
            return spec, t, pk

        def spectrum():
            spec, _, _ = state()
            keep = spec > 1e-4
            return fr2.curve(lw[keep], spec[keep], P.SUN, stroke_width=3.5)

        def peak_marker():
            _, _, pk = state()
            x = np.log10(pk)
            return VGroup(Line(fr2.p(x, 0), fr2.p(x, 1.0), color=P.TEAL, stroke_width=2),
                          Dot(fr2.p(x, 1.0), radius=0.06, color=P.TEAL))

        def temp_readout():
            _, t, pk = state()
            s = layout.label(f"T = {t:,.0f} K     glow peaks at {pk:,.1f} µm" if pk < 10 else
                             f"T = {t:,.0f} K     glow peaks at {pk:,.0f} µm",
                             font_size=22, color=P.FG)
            return s.move_to([fr2.c[0], 2.62, 0])

        spec_m = always_redraw(spectrum)
        peak_m = always_redraw(peak_marker)
        read_m = always_redraw(temp_readout)
        cap = _say(self, "The Sun (5778 K) glows brightest in the visible, where our eyes work.",
                   hold=False)
        self.play(Create(spec_m), FadeIn(peak_m), FadeIn(read_m), run_time=1.2)
        timing.hold_to_read(self, cap, settle=0.5)
        i_earth = float(np.interp(np.log(288.0), np.log(temps[::-1]), np.arange(len(temps))[::-1]))
        cap = _say(self, "Cool it to Earth's 288 K and the glow slides into the infrared ...",
                   cap, hold=False)
        self.play(idx.animate.set_value(i_earth), run_time=3.0, rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap, settle=0.3)
        cap = _say(self, "... and at 40 K it peaks beyond 100 µm, the far infrared of IRAS and AKARI.",
                   cap, hold=False)
        self.play(idx.animate.set_value(len(temps) - 1), run_time=3.0,
                  rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap, settle=1.0)
        self.play(FadeOut(VGroup(axes2, bands, spec_m, peak_m, read_m)), FadeOut(cap))

        # ---- beat 3: reflected vs thermal fall-off --------------------------
        fd = np.asarray(d["flux_distance_au"])
        refl = np.asarray(d["reflected_rel"])
        therm = np.asarray(d["thermal_rel"])
        fr3 = _Frame(2.0, np.log10(1500), -5.2, 0.2, 7.6, 4.5, (-1.7, 0.55))
        axes3 = fr3.axes(
            xticks=[(np.log10(v), str(v)) for v in (100, 200, 500, 1000)],
            yticks=[(float(k), "1" if k == 0 else f"10{_sup(k)}") for k in range(0, -6, -1)],
            xlabel="distance from the Sun (AU, log scale)",
            ylabel="flux relative to 100 AU",
        )
        c_refl = fr3.curve(np.log10(fd), np.log10(refl), P.SUN, stroke_width=3.5)
        c_th = fr3.curve(np.log10(fd), np.log10(therm), P.BLUE, stroke_width=3.5)
        l_refl = VGroup(
            layout.label("reflected sunlight ∝ 1/r⁴", font_size=19, color=P.SUN),
            layout.label(f"÷ {1 / refl[-1]:,.0f} by {fd[-1]:.0f} AU", font_size=18, color=P.SUN),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        l_refl.next_to(c_refl.get_end(), RIGHT, buff=0.2)
        l_th = VGroup(
            layout.label("its own 40 K glow ∝ 1/r²", font_size=19, color=P.BLUE),
            layout.label(f"÷ {1 / therm[-1]:,.0f} by {fd[-1]:.0f} AU", font_size=18, color=P.BLUE),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.08)
        l_th.next_to(c_th.get_end(), RIGHT, buff=0.2)
        self.play(FadeIn(axes3), run_time=0.6)
        cap = _say(self, "Reflected light is thinned twice; the heat glow only once.", hold=False)
        self.play(Create(c_refl), Create(c_th), run_time=2.0)
        self.play(FadeIn(l_refl), FadeIn(l_th))
        timing.hold_to_read(self, cap, settle=0.6)
        cap = _say(self, "Far out, the heat signal holds up far better than the reflection.", cap)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Its own 40 K glow peaks near 100 µm and fades only as 1/r².")


def _sup(k):
    """Unicode superscript for a small signed integer exponent."""
    table = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")
    return str(k).translate(table)


# ---------------------------------------------------------------------------
# P11 -- finding a slow mover
# ---------------------------------------------------------------------------

class P11Movers(Scene):
    """Goal: distant bodies are found by their motion against the fixed stars;
    near opposition the drift is mostly Earth's parallax, so it scales as
    1/distance -- the crawl rate is a distance tag."""

    def construct(self):
        d = _data()["movers"]
        anchors = {a["name"]: a for a in d["anchors"]}
        sed, p9 = anchors["Sedna"], anchors["Planet Nine"]
        _header(self, "PREFACE 11", "Finding a slow mover")

        # ---- beat 1: blink three nights -----------------------------------
        scale = 0.1  # scene units per arcsec
        field = Rectangle(width=11.0, height=4.6, color=P.PURPLE, stroke_width=1.6)
        field.set_fill("#16171f", opacity=1.0).move_to([0, 0.35, 0])
        stars = widgets.star_field(n=70, seed=5, x=(-5.3, 5.3), y=(-1.7, 2.4), color=P.FG,
                                   radius=0.03, opacity=0.85)
        bar = Line([3.4, -1.5, 0], [3.4 + 20 * scale, -1.5, 0], color=P.FG, stroke_width=3)
        bar_l = layout.label("20 arcsec", font_size=17, color=P.FG).next_to(bar, UP, buff=0.08)
        self.play(FadeIn(field), FadeIn(stars, lag_ratio=0.01), FadeIn(bar), FadeIn(bar_l),
                  run_time=1.0)

        sed_step = sed["arcsec_per_day"] * scale
        p9_step = p9["arcsec_per_day"] * scale
        sed_pos = [np.array([-4.6 + k * sed_step, -0.75, 0]) for k in range(3)]
        p9_pos = [np.array([-0.9 + k * p9_step, 1.35, 0]) for k in range(3)]
        m_sed = Dot(sed_pos[0], radius=0.07, color=P.GREEN)
        m_p9 = Dot(p9_pos[0], radius=0.07, color=P.BLUE)
        night = layout.label("night 1", font_size=21, color=P.PURPLE)
        night.move_to(field.get_corner(UP + LEFT) + np.array([0.75, -0.32, 0]))
        self.add(m_sed, m_p9, night)
        cap = _say(self, "Take one patch of sky every night and blink between the images.",
                   hold=False)
        for k in (1, 2, 0, 1, 2):
            self.wait(0.7)
            m_sed.move_to(sed_pos[k])
            m_p9.move_to(p9_pos[k])
            new = layout.label(f"night {k + 1}", font_size=21, color=P.PURPLE).move_to(night)
            night.become(new)
        self.wait(0.7)

        ghosts = VGroup(*[Dot(p, radius=0.06, color=P.GREEN).set_opacity(0.45) for p in sed_pos],
                        *[Dot(p, radius=0.06, color=P.BLUE).set_opacity(0.45) for p in p9_pos])
        tracks = VGroup(
            Arrow(sed_pos[0], sed_pos[2], buff=0.12, color=P.GREEN, stroke_width=2.5),
            Arrow(p9_pos[0], p9_pos[2] + RIGHT * 0.35, buff=0.12, color=P.BLUE, stroke_width=2.5,
                  max_tip_length_to_length_ratio=0.3),
        )
        l_sed = layout.label(f"a body at {sed['r_au']:.0f} AU: {sed['arcsec_per_day']:.0f}″ per day",
                             font_size=19, color=P.GREEN).next_to(sed_pos[1], DOWN, buff=0.3)
        l_p9 = layout.label(f"Planet Nine at {p9['r_au']:.0f} AU: {p9['arcsec_per_day']:.1f}″ per day",
                            font_size=19, color=P.BLUE).next_to(p9_pos[1], UP, buff=0.3)
        cap = _say(self, "Stars stay put; anything that moves belongs to the solar system.", cap,
                   hold=False)
        self.play(stars.animate.set_opacity(0.25), FadeIn(ghosts), Create(tracks),
                  FadeIn(l_sed), FadeIn(l_p9), run_time=1.2)
        timing.hold_to_read(self, cap, settle=1.0)
        self.play(FadeOut(VGroup(field, stars, bar, bar_l, m_sed, m_p9, night, ghosts, tracks,
                                 l_sed, l_p9)), FadeOut(cap))

        # ---- beat 2: why it moves -- Earth's parallax -----------------------
        sun_pt = np.array([-4.2, -1.95, 0])
        orbit = Circle(radius=0.55, color=P.MUTED, stroke_width=1.5).move_to(sun_pt)
        sun = Dot(sun_pt, radius=0.1, color=P.SUN)
        y_star = 2.6
        backdrop = Line([-6.7, y_star, 0], [-1.7, y_star, 0], color=P.MUTED, stroke_width=1.2)
        bstars = VGroup(*[Dot([x, y_star + 0.1 * np.sin(7 * x), 0], radius=0.025, color=P.FG)
                          .set_opacity(0.7) for x in np.linspace(-6.5, -1.9, 15)])
        b_lab = layout.label("distant stars", font_size=17, color=P.MUTED)
        b_lab.next_to(backdrop, DOWN, buff=0.2).align_to(backdrop, RIGHT)
        scale_note = layout.label("(not to scale)", font_size=15, color=P.MUTED)
        scale_note.next_to(backdrop, DOWN, buff=0.2).align_to(backdrop, LEFT)
        near = np.array([-4.5, 0.0, 0])
        far = np.array([-4.2, 1.5, 0])
        n_dot = Dot(near, radius=0.08, color=P.GREEN)
        f_dot = Dot(far, radius=0.08, color=P.BLUE)
        n_lab = layout.label("nearer", font_size=18, color=P.GREEN).next_to(n_dot, RIGHT, buff=0.15)
        f_lab = layout.label("farther", font_size=18, color=P.BLUE).next_to(f_dot, RIGHT, buff=0.15)
        theta = ValueTracker(np.deg2rad(125))

        def earth():
            th = theta.get_value()
            return sun_pt + 0.55 * np.array([np.cos(th), np.sin(th), 0])

        def sight(obj, col):
            e = earth()
            t = (y_star - e[1]) / (obj[1] - e[1])
            hit = e + t * (obj - e)
            return VGroup(DashedLine(e, hit, color=col, stroke_width=1.5, dash_length=0.08),
                          Dot(hit, radius=0.07, color=col))

        e_dot = always_redraw(lambda: Dot(earth(), radius=0.08, color=P.TEAL))
        s_near = always_redraw(lambda: sight(near, P.GREEN))
        s_far = always_redraw(lambda: sight(far, P.BLUE))
        e_lab = layout.label("Earth", font_size=17, color=P.TEAL).next_to(orbit, RIGHT, buff=0.1)
        cap = _say(self, "Why they move: Earth itself races around the Sun at 30 km/s.",
                   hold=False)
        self.play(FadeIn(VGroup(orbit, sun, backdrop, bstars, b_lab, n_dot, f_dot, n_lab, f_lab,
                                e_lab, scale_note)), FadeIn(e_dot), run_time=1.0)
        self.play(Create(s_near), Create(s_far), run_time=0.8)
        self.play(theta.animate.set_value(np.deg2rad(55)), run_time=3.0,
                  rate_func=rate_functions.there_and_back_with_pause)
        cap = _say(self, "Seen from a moving Earth, nearer bodies swing farther against the stars.",
                   cap, hold=False)
        self.play(theta.animate.set_value(np.deg2rad(55)), run_time=2.5,
                  rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap, settle=0.3)

        # ---- beat 3: rate vs distance --------------------------------------
        dist = np.asarray(d["distance_au"])
        rate = np.asarray(d["arcsec_per_day"])
        fr = _Frame(np.log10(20), np.log10(1500), 0.0, np.log10(300), 5.4, 4.4, (3.75, 0.55))
        axes = fr.axes(
            xticks=[(np.log10(v), str(v)) for v in (20, 50, 100, 200, 500, 1000)],
            yticks=[(np.log10(v), str(v)) for v in (1, 3, 10, 30, 100, 300)],
            xlabel="distance (AU)",
            ylabel="drift (arcsec per day)",
        )
        c = fr.curve(np.log10(dist), np.log10(rate), P.FG, stroke_width=3)
        cap = _say(self, "So the drift rate falls as 1/distance: it is a distance tag.", cap,
                   hold=False)
        self.play(FadeIn(axes), Create(c), run_time=1.6)
        pts = VGroup()
        for a in d["anchors"]:
            col = P.BLUE if a["name"] == "Planet Nine" else P.GREEN
            dot = Dot(fr.p(np.log10(a["r_au"]), np.log10(a["arcsec_per_day"])), radius=0.07,
                      color=col)
            if a["name"] == "Planet Nine":
                lab = VGroup(
                    layout.label(f"Planet Nine, {a['r_au']:.0f} AU:", font_size=18, color=col),
                    layout.label(f"{a['arcsec_per_day']:.1f}″ per day", font_size=18, color=col),
                ).arrange(DOWN, aligned_edge=RIGHT, buff=0.06)
                lab.next_to(dot, DOWN + LEFT, buff=0.08)
            else:
                lab = layout.label(f"{a['name']}  {a['arcsec_per_day']:.0f}″/day", font_size=18,
                                   color=col).next_to(dot, RIGHT, buff=0.12)
            pts.add(VGroup(dot, lab))
        self.play(FadeIn(pts, lag_ratio=0.4), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.5)
        cap = _say(self, f"At that pace Planet Nine needs {p9['days_per_moon_width']:.0f} days "
                         "to crawl the width of the full Moon.", cap)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Its slow crawl against the stars gives away its distance.")


# ---------------------------------------------------------------------------
# P11b -- selection bias: surveys find what they point at
# ---------------------------------------------------------------------------

_R_MAX_AU = 1400.0
_R_DISP = 2.4


def _rho(r_au):
    """Square-root radial compression so ~60 AU and ~1000 AU both fit."""
    return _R_DISP * np.sqrt(np.asarray(r_au) / _R_MAX_AU)


def _orbit_path(o, centre, color, stroke_width, opacity):
    nu = np.linspace(-np.pi, np.pi, 160)
    a, e = o["a"], o["e"]
    r = a * (1 - e * e) / (1 + e * np.cos(nu))
    ang = np.deg2rad(o["varpi_deg"]) + nu
    rr = _rho(np.minimum(r, _R_MAX_AU))
    pts = [centre + np.array([q * np.cos(t), q * np.sin(t), 0]) for q, t in zip(rr, ang)]
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners(pts)
    return m.set_stroke(opacity=opacity)


def _body_point(o, centre):
    q = _rho(min(o["r_au"], _R_MAX_AU))
    t = np.deg2rad(o["lon_deg"])
    return centre + np.array([q * np.cos(t), q * np.sin(t), 0])


class P11SelectionBias(Scene):
    """Goal: a survey only finds distant objects where it looks and only near
    perihelion (when they are bright), so the perihelion directions it finds
    echo its own pointing -- even from an isotropic population."""

    def construct(self):
        d = _data()["bias"]
        _header(self, "PREFACE 11b", "Selection bias: you find what you point at")
        centre = np.array([-3.6, 0.05, 0])

        # ---- an isotropic population ---------------------------------------
        bg = d["background"]
        orbits_bg = VGroup(*[_orbit_path(o, centre, P.MUTED, 0.9, 0.35) for o in bg])
        bodies_bg = VGroup(*[Dot(_body_point(o, centre), radius=0.025, color=P.FG).set_opacity(0.6)
                             for o in bg])
        sun = Dot(centre, radius=0.09, color=P.SUN).set_z_index(5)
        note = VGroup(layout.label("distances", font_size=16, color=P.MUTED),
                      layout.label("√-scaled", font_size=16, color=P.MUTED)).arrange(DOWN, buff=0.06)
        note.move_to([-6.35, 2.45, 0])
        cap = _say(self, f"Simulate {d['n']:,} distant objects with perihelia pointing every which "
                         "way.", hold=False)
        self.play(FadeIn(sun), Create(orbits_bg, lag_ratio=0.01), FadeIn(note), run_time=2.0)
        self.play(FadeIn(bodies_bg), run_time=0.6)
        timing.hold_to_read(self, cap, settle=0.4)

        # ---- legend (top right) and the depth limit ---------------------------
        def key(swatch, text):
            return VGroup(swatch, layout.label(text, font_size=17, color=P.FG)).arrange(
                RIGHT, buff=0.18)

        r_lim = _rho(d["r_limit_au"])
        disc = Circle(radius=r_lim, color=P.PURPLE, stroke_width=2).move_to(centre)
        disc.set_fill(P.PURPLE, opacity=0.08)
        keys = VGroup(
            key(Circle(radius=0.1, color=P.PURPLE, stroke_width=2),
                f"bright enough to see (V < {d['depth']:.0f}): within {d['r_limit_au']:.0f} AU"),
            key(AnnularSector(inner_radius=0, outer_radius=0.24, angle=0.8, start_angle=-0.4,
                              color=P.PURPLE, fill_opacity=0.35, stroke_width=0),
                "where the survey looks"),
            key(Line(LEFT * 0.1, RIGHT * 0.1, color=P.GREEN, stroke_width=3).rotate(PI / 2),
                "perihelion direction of each find"),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.12)
        keys.move_to([3.5, 2.2, 0])
        cap = _say(self, "Faint with distance (1/r⁴), they can be seen only near the Sun,"
                         " so near perihelion.", cap, hold=False)
        self.play(Create(disc), FadeIn(keys[0]))
        timing.hold_to_read(self, cap, settle=0.4)

        # ---- right: histogram of perihelion longitude ------------------------
        bins = d["bins"]
        all_share = np.asarray(d["all_hist"], dtype=float) / sum(d["all_hist"])
        fr = _Frame(0.0, 360.0, 0.0, 0.3, 5.2, 3.1, (3.85, -0.2))
        axes = fr.axes(
            xticks=[(v, f"{v}°") for v in (0, 90, 180, 270, 360)],
            yticks=[(0.0, "0"), (0.1, "10%"), (0.2, "20%"), (0.3, "30%")],
            xlabel="longitude of perihelion ϖ",
            ylabel="share of objects",
        )
        edges = np.linspace(0, 360, bins + 1)

        def bars(share, color, opacity):
            g = VGroup()
            for k, s in enumerate(share):
                if s <= 0:
                    continue
                a, b = fr.p(edges[k], 0), fr.p(edges[k + 1], min(s, 0.3))
                g.add(Rectangle(width=b[0] - a[0], height=b[1] - a[1], stroke_width=0.6,
                                color=color).set_fill(color, opacity=opacity)
                      .move_to((a + b) / 2))
            return g

        h_all = bars(all_share, P.MUTED, 0.5)
        r_all = layout.label(f"all simulated: R̄ = {d['all_r_bar']:.2f}", font_size=18,
                             color=P.FG)
        r_all.move_to(fr.p(355, 0.29), aligned_edge=UP + RIGHT)
        self.play(FadeIn(axes), FadeIn(h_all), FadeIn(r_all), run_time=1.0)

        # ---- the survey looks in one direction --------------------------------
        def wedge_group(w):
            c, hw = w["centre_deg"], w["half_width_deg"]
            sector = AnnularSector(inner_radius=0, outer_radius=_R_DISP + 0.1,
                                   angle=np.deg2rad(2 * hw), start_angle=np.deg2rad(c - hw),
                                   color=P.PURPLE, fill_opacity=0.14, stroke_width=0,
                                   arc_center=centre)
            found = w["detected"]
            orbs = VGroup(*[_orbit_path(o, centre, P.GREEN, 1.1, 0.5) for o in found])
            dots = VGroup(*[Dot(_body_point(o, centre), radius=0.05, color=P.GREEN)
                            .set_z_index(4) for o in found])
            rug = VGroup()
            for o in found:
                t = np.deg2rad(o["varpi_deg"])
                u = np.array([np.cos(t), np.sin(t), 0])
                rug.add(Line(centre + u * (_R_DISP + 0.08), centre + u * (_R_DISP + 0.3),
                             color=P.GREEN, stroke_width=2))
            share = np.asarray(w["hist"], dtype=float) / max(1, sum(w["hist"]))
            hist = bars(share, P.GREEN, 0.7)
            rbar = layout.label(f"found ({len(found)}): R̄ = {w['r_bar']:.2f}", font_size=18,
                                color=P.GREEN)
            rbar.next_to(r_all, DOWN, buff=0.12).align_to(r_all, RIGHT)
            return sector, VGroup(orbs, dots, rug), hist, rbar

        w1, w2 = d["wedges"]
        sec1, found1, hist1, rb1 = wedge_group(w1)
        cap = _say(self, "A survey watches one patch of sky. What does it find?", cap, hold=False)
        self.play(FadeIn(sec1), FadeIn(keys[1]), run_time=0.8)
        self.play(orbits_bg.animate.set_stroke(opacity=0.15), bodies_bg.animate.set_opacity(0.3),
                  Create(found1[0], lag_ratio=0.02), FadeIn(found1[1]), run_time=1.6)
        self.play(Create(found1[2], lag_ratio=0.02), FadeIn(keys[2]), FadeIn(hist1), FadeIn(rb1),
                  run_time=1.2)
        cap = _say(self, "Caught near perihelion, their perihelia point into the field: a 'cluster'.",
                   cap)

        sec2, found2, hist2, rb2 = wedge_group(w2)
        cap = _say(self, "Point the survey elsewhere, and the cluster follows the telescope.", cap,
                   hold=False)
        self.play(Transform(sec1, sec2), FadeOut(found1),
                  Transform(hist1, hist2), Transform(rb1, rb2), run_time=1.6)
        self.play(Create(found2[0], lag_ratio=0.02), FadeIn(found2[1]),
                  Create(found2[2], lag_ratio=0.02), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.8)
        cap = _say(self, "Real tests must model each survey's pointing and depth before trusting "
                         "a cluster.", cap)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Surveys find what they point at: model the bias first.")
