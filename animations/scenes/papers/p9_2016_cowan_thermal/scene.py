"""Cowan, Holder & Kaib (2016) -- the case for finding Planet Nine with CMB
experiments.

A Neptune-sized body at 40 K shines by reflected sunlight only in the optical
and near-infrared; beyond about 16 µm its own heat takes over, peaking near
70 µm and still delivering several mJy at 1-3 mm, where cosmology experiments
map large areas of sky again and again. Reproduced in p9-2016-cowan-thermal;
the spectrum and the flux-distance curves are the crate's own
(anim.json -> papers -> p9-2016-cowan-thermal).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
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

CRATE = "p9-2016-cowan-thermal"


def power_of_ten(value):
    """'10⁻⁶' style tick text for a log axis."""
    digits = str(int(round(np.log10(value)))).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return "10" + digits


def round_to_one_figure(value):
    scale = 10 ** np.floor(np.log10(value))
    return float(np.round(value / scale) * scale)


class CowanThermal2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        sed = d["sed"]
        lam = 10 ** np.array(sed["log_wavelength_um"])
        refl = 10 ** np.array(sed["log_reflected_mjy"])
        therm = 10 ** np.array(sed["log_thermal_mjy"])

        self.add(paper.scene_header(CRATE))

        # 1. the spectrum: reflected sunlight, then the planet's own heat
        plot = Plot([0.3, 1e4], [1e-8, 1e3], [1, 10, 100, 1000, 10000],
                    [1e-8, 1e-6, 1e-4, 1e-2, 1, 100], "wavelength (µm)", "flux density (mJy)",
                    x_log=True, y_log=True, y_fmt=power_of_ten)
        self.play(FadeIn(plot))

        c_refl = plot.curve(lam, refl, P.SUN)
        c_therm = plot.curve(lam, therm, P.BLUE)
        key = VGroup(
            layout.label("reflected sunlight", font_size=16, color=P.SUN),
            layout.label(f"its own heat, {d['temp_k']:.0f} K", font_size=16, color=P.BLUE),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        key.move_to(plot.p(2.2, 60.0))

        cap = layout.caption(
            f"A {d['mass_earth']:.0f} Earth-mass planet at {d['distance_au']:.0f} AU: "
            "in visible light it only reflects the Sun", font_size=22)
        self.play(Create(c_refl), FadeIn(key[0]), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.6)

        gain = round_to_one_figure(therm.max() / refl.max())
        cross = plot.vline(d["crossover_um"], P.MUTED, y_to=1.0)
        cross_lab = layout.label(f"heat wins beyond {d['crossover_um']:.0f} µm", font_size=16,
                                 color=P.MUTED)
        cross_lab.next_to(plot.p(d["crossover_um"], 1e-7), RIGHT, buff=0.12)
        peak = Dot(plot.p(d["wien_peak_um"], np.interp(d["wien_peak_um"], lam, therm)),
                   radius=0.06, color=P.BLUE)
        cap2 = layout.caption(
            f"Its own heat peaks near {d['wien_peak_um']:.0f} µm, "
            f"{gain:,.0f} times brighter than the reflected light", font_size=22)
        self.play(Create(c_therm), FadeIn(key[1]), FadeOut(cap), FadeIn(cap2), run_time=1.6)
        self.play(Create(cross), FadeIn(cross_lab), FadeIn(peak))
        timing.hold_to_read(self, cap2, cross_lab, settle=0.6)

        bands = d["cmb_bands"]
        w_lo = min(b["wavelength_um"] for b in bands) / 1.25
        w_hi = max(b["wavelength_um"] for b in bands) * 1.25
        window = plot.band(w_lo, w_hi, 1e-8, 1e3, P.PURPLE, opacity=0.18)
        win_lab = layout.label("CMB surveys", font_size=16, color=P.PURPLE)
        win_lab.next_to(window, UP, buff=0.06)
        band_dots = VGroup(*[
            Dot(plot.p(b["wavelength_um"], b["flux_mjy"]), radius=0.05, color=P.PURPLE)
            for b in bands])
        f_lo = min(b["flux_mjy"] for b in bands)
        f_hi = max(b["flux_mjy"] for b in bands)
        cap3 = layout.caption(
            f"CMB telescopes map the sky at 1-3 mm, where it is still {f_lo:.0f}-{f_hi:.0f} mJy",
            font_size=22)
        self.play(FadeIn(window), FadeIn(win_lab), FadeIn(band_dots), FadeOut(cap2),
                  FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.8)
        self.play(FadeOut(VGroup(plot, c_refl, c_therm, key, cross, cross_lab, peak, window,
                                 win_lab, band_dots, cap3)))

        # 2. how bright at 1 mm, and how far out
        fd = d["flux_vs_distance"]
        dist = np.array(fd["distance_au"])
        plot2 = Plot([200, 1400], [0.5, 300], [200, 400, 600, 800, 1000, 1200, 1400],
                     [1, 3, 10, 30, 100], "distance from the Sun (AU)",
                     "flux density at 1 mm (mJy)", y_log=True)
        self.play(FadeIn(plot2))

        fid = plot2.curve(dist, fd["fiducial_mjy"], P.BLUE)
        faint = plot2.curve(dist, fd["faint_mjy"], P.TEAL)
        fid_lab = layout.label(
            f"{d['mass_earth']:.0f} M⊕ at {d['temp_k']:.0f} K", font_size=16, color=P.BLUE)
        fid_lab.next_to(plot2.p(dist[48], fd["fiducial_mjy"][48]), UP, buff=0.25)
        faint_lab = layout.label(
            f"{fd['faint_mass_earth']:.0f} M⊕ at {fd['faint_temp_k']:.0f} K", font_size=16,
            color=P.TEAL)
        faint_lab.next_to(plot2.p(dist[24], fd["faint_mjy"][24]), DOWN, buff=0.28)

        def level(value, text, side):
            line = plot2.hline(value, P.FG)
            lab = layout.label(text, font_size=16, color=P.FG)
            lab.next_to(line.get_end(), side, buff=0.08).align_to(line.get_end(), RIGHT)
            return VGroup(line, lab)

        pub = level(d["published_flux_1mm_mjy"],
                    f"paper: {d['published_flux_1mm_mjy']:.0f} mJy at "
                    f"{d['distance_au']:.0f} AU", UP)
        floor = level(d["published_faint_flux_1mm_mjy"],
                      f"paper: as faint as {d['published_faint_flux_1mm_mjy']:.0f} mJy", DOWN)
        cap4 = layout.caption("At 1 mm the glow fades only as distance squared", font_size=22)
        self.play(Create(fid), Create(faint), FadeIn(fid_lab), FadeIn(faint_lab), FadeIn(cap4),
                  run_time=1.4)
        timing.hold_to_read(self, cap4, settle=0.5)

        here = Dot(plot2.p(d["distance_au"], d["flux_1mm_mjy"]), radius=0.08, color=P.GREEN)
        here_lab = layout.label(f"reproduced: {d['flux_1mm_mjy']:.0f} mJy", font_size=16,
                                color=P.GREEN)
        here_lab.next_to(here, UP, buff=0.14).shift(RIGHT * 0.9)
        cap5 = layout.caption(
            f"The blackbody model gives {d['flux_1mm_mjy']:.0f} mJy where the paper quotes "
            f"{d['published_flux_1mm_mjy']:.0f}", font_size=22)
        self.play(FadeIn(pub), FadeIn(floor), FadeIn(here), FadeIn(here_lab), FadeOut(cap4),
                  FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Planet Nine glows at a millimetre, where CMB surveys already scan the sky.")
