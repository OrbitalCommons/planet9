"""Fortney et al. (2016) -- atmosphere, spectra, evolution and detectability.

Interior models leave Planet Nine at no more than 35-50 K, held there by its
own heat rather than by sunlight. At that temperature methane condenses, and
the clear hydrogen atmosphere lets 3-5 µm light escape some twenty orders of
magnitude above what a blackbody of the same temperature emits. Reproduced in
p9-2016-fortney-thermal, whose energy balance (with the heat flow scaled to the
paper's models) gives the temperatures and the blackbody spectrum shown here
(anim.json -> papers -> p9-2016-fortney-thermal).
"""
import numpy as np
from manim import (
    DOWN,
    RIGHT,
    UP,
    Arrow,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim.plot import Plot
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2016-fortney-thermal"


def power_of_ten(value):
    """'10⁻⁶' style tick text for a log axis."""
    digits = str(int(round(np.log10(value)))).translate(
        str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return "10" + digits


class FortneyThermal2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]

        self.add(paper.scene_header(CRATE))

        # 1. how warm: internal heat against sunlight
        tv = d["temp_vs_distance"]
        plot = Plot([0, 1000], [0, 60], [0, 200, 400, 600, 800, 1000], [0, 10, 20, 30, 40, 50, 60],
                    "distance from the Sun (AU)", "temperature (K)")
        band = plot.band(0, 1000, pub["teff_min_k"], pub["teff_max_k"], P.FG, opacity=0.10)
        band_lab = layout.label(
            f"paper: at most {pub['teff_min_k']:.0f}-{pub['teff_max_k']:.0f} K for "
            f"{pub['mass_min_earth']:.0f}-{pub['mass_max_earth']:.0f} M⊕",
            font_size=15, color=P.FG)
        band_lab.move_to(plot.p(700, pub["teff_max_k"] + 3.5))
        sun = plot.curve(tv["distance_au"], tv["solar_only_k"], P.SUN)
        sun_lab = layout.label("heated by sunlight alone", font_size=16, color=P.SUN)
        sun_lab.next_to(plot.p(500, d["t_solar_10me"]), DOWN, buff=0.18)
        self.play(FadeIn(plot))
        cap = layout.caption(
            f"Sunlight alone would leave Planet Nine at {d['t_solar_10me']:.0f} K "
            f"at {d['distance_au']:.0f} AU", font_size=22)
        self.play(Create(sun), FadeIn(sun_lab), FadeIn(cap), run_time=1.3)
        timing.hold_to_read(self, cap, settle=0.5)

        own = plot.curve(tv["distance_au"], tv["with_internal_heat_k"], P.BLUE)
        own_lab = layout.label("10 M⊕ with its internal heat", font_size=16, color=P.BLUE)
        own_lab.next_to(plot.p(260, d["teff_10me"]), DOWN, buff=0.14)
        others = VGroup()
        for key, mass, side in (("teff_5me", pub["mass_min_earth"], DOWN),
                                ("teff_20me", pub["mass_max_earth"], UP)):
            dot = Dot(plot.p(d["distance_au"], d[key]), radius=0.05, color=P.BLUE)
            lab = layout.label(f"{mass:.0f} M⊕: {d[key]:.0f} K", font_size=15, color=P.BLUE)
            lab.next_to(dot, side, buff=0.08).shift(RIGHT * 0.75)
            others.add(dot, lab)
        cap2 = layout.caption(
            f"Heat left from its formation holds it at {d['teff_10me']:.0f} K, "
            f"{d['internal_to_solar']:.0f} times the sunlight it absorbs", font_size=22)
        self.play(FadeIn(band), FadeIn(band_lab), Create(own), FadeIn(own_lab), FadeIn(others),
                  FadeOut(cap), FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, band_lab, settle=0.8)
        self.play(FadeOut(VGroup(plot, band, band_lab, sun, sun_lab, own, own_lab, others, cap2)))

        # 2. the spectrum of a body that cold
        sed = d["sed"]
        lam = 10 ** np.array(sed["log_wavelength_um"])
        flux = 10 ** np.array(sed["log_flux_jy"])
        plot2 = Plot([1, 3000], [1e-45, 1e3], [1, 10, 100, 1000],
                     [1e-40, 1e-30, 1e-20, 1e-10, 1], "wavelength (µm)", "flux density (Jy)",
                     x_log=True, y_log=True, y_fmt=power_of_ten)
        body = plot2.curve(lam, flux, P.BLUE)
        body_lab = layout.label(f"blackbody at {d['teff_10me']:.0f} K", font_size=16,
                                color=P.BLUE)
        body_lab.next_to(plot2.p(d["sed_peak_um"], d["far_ir_flux_jy"]), DOWN, buff=0.35)
        self.play(FadeIn(plot2))
        cap3 = layout.caption(
            f"A {d['teff_10me']:.0f} K blackbody peaks at {d['sed_peak_um']:.0f} µm "
            "and vanishes in the near-infrared", font_size=22)
        self.play(Create(body), FadeIn(body_lab), FadeIn(cap3), run_time=1.5)
        timing.hold_to_read(self, cap3, settle=0.5)

        wise = VGroup()
        for band_key, name in (("w1", "W1"), ("w2", "W2")):
            w = d[band_key]
            wise.add(Dot(plot2.p(w["wavelength_um"], 10 ** w["log_flux_jy"]), radius=0.06,
                         color=P.PURPLE))
        window = plot2.band(3.0, 5.2, 1e-45, 1e3, P.PURPLE, opacity=0.16)
        win_lab = layout.label("WISE W1, W2", font_size=15, color=P.PURPLE)
        win_lab.next_to(window, UP, buff=0.06)
        self.play(FadeIn(window), FadeIn(win_lab), FadeIn(wise))

        lift = pub["near_ir_excess_dex"]
        w2 = d["w2"]
        start = plot2.p(w2["wavelength_um"], 10 ** w2["log_flux_jy"])
        end = plot2.p(w2["wavelength_um"], 10 ** (w2["log_flux_jy"] + lift))
        arrow = Arrow(start, end, color=P.GREEN, buff=0.08, stroke_width=4,
                      max_tip_length_to_length_ratio=0.12)
        lift_lab = layout.label(
            f"paper: model atmospheres\nare ~{lift:.0f} orders of magnitude\nbrighter here",
            font_size=15, color=P.GREEN, line_spacing=0.9)
        lift_lab.move_to(plot2.p(200, 1e-27))
        cap4 = layout.caption(
            "With its methane frozen out, the real atmosphere leaks 3-5 µm light", font_size=22)
        self.play(Create(arrow), FadeIn(lift_lab), FadeOut(cap3), FadeIn(cap4))
        timing.hold_to_read(self, cap4, lift_lab, settle=1.0)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, "Cold but not dark: a methane-poor Planet Nine could just show up in WISE.")
