"""Preface E -- Seeing the unseen: reflected light, thermal glow, and surveys.

Scenes: P09ReflectedLight, P10Thermal, P11Surveys.
"""
import numpy as np
from manim import (
    Axes,
    Create,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    LEFT,
    Line,
    RIGHT,
    Rectangle,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
    always_redraw,
    rate_functions,
    ValueTracker,
)

import p9_manim as P
from p9_manim import layout, orbits, dataio, timing


def _reflected_mag(mass_earth, r_au):
    """Opposition reflected V magnitude (delta ~ r-1). Documented fallback model
    consistent with figures/p9_magnitude_distance (Neptune albedo, m ~ 1/r^4)."""
    # H tuned so a ~10 M_earth body is ~mag 22 at 600 AU.
    H = -6.0 - 1.5 * np.log10(mass_earth / 10.0)
    delta = max(r_au - 1.0, 1.0)
    return H + 5.0 * np.log10(r_au * delta)


def _mag_curve(mass_earth):
    """Real (distance_au, v_mag) points for a mass from viability().mag_curves.

    Returns the curve whose mass_earth is closest to the request; falls back to
    the inline reflected-light model if data is unavailable.
    """
    try:
        v = dataio.viability()
        dist = v["distance_au"]
        curves = v["mag_curves"]
        if dist and curves:
            best = min(curves, key=lambda c: abs(c["mass_earth"] - mass_earth))
            return list(zip(dist, best["v_mag"]))
    except Exception:
        pass
    return [(r, _reflected_mag(mass_earth, r)) for r in np.arange(100, 1500.1, 10)]


class P09ReflectedLight(Scene):
    def construct(self):
        self.add(layout.concept_badge("PREFACE 09"))
        title = Text("Seeing it I: borrowed sunlight", color=P.FG, font_size=30, weight="BOLD")
        title.to_edge(UP, buff=0.55)
        self.play(Write(title))
        flux = layout.equation(
            r"F \propto \dfrac{p\,R^2}{r^2\Delta^2} \approx \dfrac{1}{r^4}",
            color=P.FG).scale(0.7).next_to(title, DOWN, buff=0.2)
        self.play(Write(flux))
        timing.hold_to_read(self, flux, settle=0.4)

        mag = layout.equation(r"m = H + 5\log_{10}(r\,\Delta)", color=P.TEAL).scale(0.6)
        mag.next_to(flux, DOWN, buff=0.15)
        self.play(Write(mag))
        timing.hold_to_read(self, mag, settle=0.4)
        self.play(FadeOut(mag))

        ax = Axes(x_range=[100, 1500, 200], y_range=[16, 28, 2],
                  x_length=9.0, y_length=4.0,
                  axis_config={"color": P.MUTED, "include_tip": False,
                               "font_size": 18, "decimal_number_config": {"num_decimal_places": 0}})
        ax.shift(DOWN * 0.6)
        xlab = ax.get_x_axis_label(Text("distance (AU)", font_size=18, color=P.FG))
        ylab = ax.get_y_axis_label(Text("apparent mag (fainter →)", font_size=18, color=P.FG))
        self.play(Create(ax), FadeIn(xlab), FadeIn(ylab))

        for mass, col in [(5, P.TEAL), (10, P.BLUE)]:
            pts = [(r, m_v) for r, m_v in _mag_curve(mass) if 100 <= r <= 1500]
            curve = ax.plot_line_graph(
                [r for r, _ in pts], [m_v for _, m_v in pts],
                line_color=col, add_vertex_dots=False)["line_graph"]
            lbl = Text(f"{mass} M⊕", font_size=18, color=col).next_to(curve.get_end(), UP, buff=0.1)
            self.play(Create(curve), FadeIn(lbl), run_time=1.0)

        for depth, name in [(20.5, "ZTF"), (23.8, "DES")]:
            y = ax.c2p(100, depth)[1]
            line = Line([ax.c2p(100, depth)[0], y, 0], [ax.c2p(1500, depth)[0], y, 0],
                        color=P.MUTED, stroke_width=1.5).set_stroke(opacity=0.7)
            tag = Text(f"{name} depth", font_size=14, color=P.MUTED).next_to(line, RIGHT, buff=0.05)
            self.play(Create(line), FadeIn(tag), run_time=0.6)

        layout.show_takeaway(self, "Brightness falls as ~1/r^4: at 600 AU, Planet Nine is ~mag 22.")


class P10Thermal(Scene):
    """Blackbody emission slides to the far-IR as the body cools (Wien's law)."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 10"))
        title = Text("Seeing it II: its own faint heat", color=P.FG, font_size=30, weight="BOLD")
        title.to_edge(UP, buff=0.55)
        self.play(Write(title))

        # Planck's law and Wien's displacement law govern the glow
        planck = layout.equation_card(
            r"B_\nu(T) = \dfrac{2h\nu^3/c^2}{e^{h\nu/kT}-1}"
            r"\qquad \lambda_{\rm peak} = \dfrac{b}{T}").scale(0.7)
        planck.next_to(title, DOWN, buff=0.2)
        layout.show_equation(self, planck)
        self.play(FadeOut(planck))

        # x axis = log10(wavelength / micron) from -1 (0.1) to 3 (1000)
        ax = Axes(x_range=[-1, 3, 1], y_range=[0, 1.1, 0.5], x_length=9.5, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.5)
        xlab = Text("wavelength (µm)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.3)
        self.play(Create(ax), FadeIn(xlab))
        # tick relabels
        for lx, txt in [(-1, "0.1"), (0, "1"), (1, "10"), (2, "100"), (3, "1000")]:
            t = Text(txt, font_size=14, color=P.MUTED).next_to(ax.c2p(lx, 0), DOWN, buff=0.12)
            self.add(t)

        # band markers
        bands = [(0.5, "optical", P.GREEN), (3.4, "WISE", P.ORANGE),
                 (75, "IRAS/AKARI", P.RED), (1000, "mm (ACT/SO)", P.PURPLE)]
        for wl, name, col in bands:
            x = np.log10(wl)
            if x < -1 or x > 3:
                continue
            ln = Line(ax.c2p(x, 0), ax.c2p(x, 1.05), color=col, stroke_width=1.5).set_stroke(opacity=0.5)
            tg = Text(name, font_size=13, color=col).next_to(ax.c2p(x, 1.05), UP, buff=0.05)
            self.add(ln, tg)

        def planck_norm(logwl_um, T):
            wl = 10 ** logwl_um * 1e-6
            h, c, k = 6.626e-34, 3.0e8, 1.381e-23
            x = h * c / (wl * k * T)
            x = min(x, 700.0)
            b = 1.0 / (wl ** 5) / (np.exp(x) - 1.0 + 1e-300)
            # normalise by peak (Wien): peak wl = 2.898e-3 / T
            wlp = 2.898e-3 / T
            xp = h * c / (wlp * k * T)
            bp = 1.0 / (wlp ** 5) / (np.exp(min(xp, 700)) - 1.0)
            return float(np.clip(b / bp, 0, 1.08))

        Tt = ValueTracker(5778.0)
        curve = always_redraw(lambda: ax.plot(
            lambda lx: planck_norm(lx, Tt.get_value()), x_range=[-1, 3, 0.05], color=P.YELLOW))
        Tlbl = always_redraw(lambda: Text(f"T = {Tt.get_value():.0f} K", font_size=22, color=P.YELLOW)
                             .to_corner(UP + RIGHT, buff=0.7))
        self.add(curve, Tlbl)
        self.wait(0.4)
        self.play(Tt.animate.set_value(40.0), run_time=4.5, rate_func=rate_functions.smooth)
        self.wait(0.3)

        data = dataio.section("thermal")
        if data:
            # overlay the REAL normalised 40 K Planck spectrum as the final state
            wl_um = np.asarray(data["wavelength_um"], dtype=float)
            pn = np.asarray(data["planck_norm_40k"], dtype=float)
            lx = np.log10(wl_um)
            mask = (lx >= -1) & (lx <= 3)
            real_curve = ax.plot_line_graph(lx[mask], pn[mask], add_vertex_dots=False, line_color=P.YELLOW)
            self.remove(curve)
            self.play(Create(real_curve), run_time=1.2)
            wien = float(data["wien_peak_um"])
            wln = Line(ax.c2p(np.log10(wien), 0), ax.c2p(np.log10(wien), 1.05),
                       color=P.TEAL, stroke_width=2).set_stroke(opacity=0.7)
            wlbl = Text(f"Wien peak {wien:.0f} µm", font_size=13, color=P.TEAL)
            wlbl.next_to(ax.c2p(np.log10(wien), 1.05), UP, buff=0.05)
            self.play(Create(wln), FadeIn(wlbl))
            self.wait(0.3)

        teq = layout.equation(
            r"T_{\rm eq} = T_\odot\sqrt{\dfrac{R_\odot}{2d}}\,(1-A)^{1/4}", color=P.TEAL).scale(0.65)
        teq.to_edge(DOWN, buff=1.1)
        self.play(Write(teq))
        timing.hold_to_read(self, teq, settle=0.4)
        self.play(FadeOut(teq))

        layout.show_takeaway(self, "A ~40 K world glows in the far-IR / mm -- a second way to find it.")


class P11Surveys(Scene):
    """Surveys: depth, footprint, finding movers, and selection bias."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 11"))
        title = Text("Finding a mover -- and fooling yourself",
                     color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.6)
        self.play(Write(title))

        rng = np.random.default_rng(11)
        stars = VGroup(*[Dot(np.array([rng.uniform(-6, 6), rng.uniform(-2.6, 2.4), 0]),
                             radius=0.02, color=P.FG).set_opacity(0.5) for _ in range(120)])
        self.play(FadeIn(stars, lag_ratio=0.005), run_time=1.2)

        # a survey footprint (covers part of the sky)
        foot = Rectangle(width=6.5, height=3.2, color=P.PURPLE, stroke_width=2).shift(LEFT * 2.3 + DOWN * 0.1)
        foot.set_fill(P.PURPLE, opacity=0.06)
        foot_lbl = Text("survey footprint (limited sky + depth)", font_size=16, color=P.PURPLE)
        foot_lbl.next_to(foot, UP, buff=0.08)
        self.play(Create(foot), FadeIn(foot_lbl))

        # a moving object detected across 3 epochs
        mover = Dot(LEFT * 4 + DOWN * 0.3, radius=0.08, color=P.GREEN)
        self.play(FadeIn(mover))
        for dx in (0.7, 0.7, 0.7):
            ghost = mover.copy().set_opacity(0.35)
            self.add(ghost)
            self.play(mover.animate.shift(RIGHT * dx), run_time=0.5)
        shift_txt = Text("a real planet shifts vs the stars -> link the epochs",
                         font_size=16, color=P.GREEN).next_to(mover, DOWN, buff=0.3)
        self.play(FadeIn(shift_txt))
        timing.hold_to_read(self, shift_txt, settle=0.4)

        # detection criterion and on-sky motion between epochs
        det_eq = layout.equation_card(
            r"\text{detect} \iff m < m_{\rm lim}"
            r"\qquad \Delta\theta = \mu\,\Delta t").scale(0.65)
        det_eq.to_edge(DOWN, buff=1.15)
        layout.show_equation(self, det_eq)
        self.play(FadeOut(det_eq), FadeOut(shift_txt))

        layout.show_takeaway(
            self, "You only find what your survey can see -- so coverage can fake (or hide) clustering.")
