"""Rowan-Robinson (2021) -- a search for Planet 9 in the IRAS data.

IRAS scanned each patch of sky on passes weeks to months apart. A body a few
hundred AU away shifts by arcminutes between passes because the Earth itself
moves, so Planet Nine would appear as a 60 µm source seen on some passes and
a lone detection nearby on another. Of several hundred such pairings one
survives: a 0.57 Jy source that moved 20 arcmin in 12 weeks. Reproduced in
p9-2022-iras-candidate: the chance-pairing estimate, the flux-distance curves
and the predicted population are the crate's own
(anim.json -> papers -> p9-2022-iras-candidate).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arrow,
    Circle,
    Create,
    Dot,
    Ellipse,
    FadeIn,
    FadeOut,
    Scene,
    Star,
    SurroundingRectangle,
    VGroup,
)

import p9_manim as P
from p9_manim.plot import Plot
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2022-iras-candidate"


def _readout(title, value, note, colour):
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


class IrasCandidate2022(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        c = d["candidate"]
        self.add(paper.scene_header(CRATE))

        # 1. why a distant planet shows up as a pair of mismatched sources
        per_arcmin = 0.105
        a = d["parallax_semi_major_arcmin"]
        b = a * np.sin(np.radians(c["ecliptic_lat_deg"]))
        centre = np.array([-2.4, 0.45, 0])
        ell = Ellipse(width=2 * a * per_arcmin, height=2 * b * per_arcmin, color=P.MUTED,
                      stroke_width=1.6).move_to(centre)
        ell_lab = layout.label(f"its yearly parallax loop at {d['published_distance_au']:.0f} AU",
                               font_size=16, color=P.MUTED)
        ell_lab.next_to(ell, DOWN, buff=0.45)
        # sweep between the passes: the chord of the loop equals the measured motion
        sweep = 2 * np.arcsin(min(1.0, c["motion_arcmin"] / (2 * a)))
        th0 = np.radians(200.0)

        def on_loop(th):
            return centre + per_arcmin * np.array([a * np.cos(th), b * np.sin(th), 0])

        p1, p3 = on_loop(th0), on_loop(th0 + sweep)
        d1 = Dot(p1, radius=0.1, color=P.GREEN)
        d3 = Dot(p3, radius=0.1, color=P.GREEN)
        l1 = layout.label("passes 1 and 2", font_size=16, color=P.GREEN).next_to(d1, LEFT, buff=0.15)
        l3 = layout.label("pass 3", font_size=16, color=P.GREEN).next_to(d3, RIGHT, buff=0.15)
        hop = Arrow(p1, p3, buff=0.12, color=P.ORANGE, stroke_width=3)
        hop_lab = layout.label(f"{c['motion_arcmin']:.0f} arcmin in {c['motion_weeks']:.0f} weeks",
                               font_size=17, color=P.ORANGE)
        hop_lab.next_to(hop.get_center(), UP + RIGHT, buff=0.15)
        moon = Circle(radius=0.5 * 31.0 * per_arcmin, color=P.FG, stroke_width=1.4)
        moon.set_fill(P.FG, opacity=0.08).move_to([3.4, 0.45, 0])
        moon_lab = layout.label("the full Moon, same scale", font_size=16, color=P.FG)
        moon_lab.next_to(moon, DOWN, buff=0.45)
        cap = layout.caption("As the Earth circles the Sun, a distant body traces a small loop",
                             font_size=22)
        self.play(Create(ell), FadeIn(ell_lab), FadeIn(moon), FadeIn(moon_lab), FadeIn(cap),
                  run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.2)
        cap2 = layout.caption("IRAS would log it in one place, then as a lone source nearby",
                              font_size=22)
        self.play(FadeIn(d1), FadeIn(l1), FadeOut(cap), FadeIn(cap2))
        self.play(Create(hop), FadeIn(d3), FadeIn(l3), FadeIn(hop_lab), run_time=1.2)
        timing.hold_to_read(self, cap2, settle=0.6)
        self.play(FadeOut(VGroup(ell, ell_lab, d1, d3, l1, l3, hop, hop_lab, moon, moon_lab,
                                 cap2)))

        # 2. hundreds of such pairings, one survivor
        n = int(d["published_associations"])
        cols = 38
        dots = VGroup(*[Dot(radius=0.045, color=P.GREEN) for _ in range(n)])
        dots.arrange_in_grid(cols=cols, buff=0.11).move_to(UP * 0.55)
        count = layout.label(
            f"{n} pairings examined by eye ({d['published_pairs']:.0f} pairs, "
            f"{d['published_triplets']:.0f} triplets, {d['published_close_pairs']:.0f} close pairs)",
            font_size=18, color=P.FG)
        count.next_to(dots, UP, buff=0.3)
        cap3 = layout.caption(
            f"Chance alone predicts {d['chance_associations']:,.0f} coincidences: "
            "most pairs are unrelated sources", font_size=22)
        self.play(FadeIn(dots, lag_ratio=0.002), FadeIn(count), FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, settle=0.3)
        keep = dots[n // 2 + cols // 2]
        others = VGroup(*[x for x in dots if x is not keep])
        cap4 = layout.caption(
            "Checking every pass in the raw scans leaves a single candidate", font_size=22)
        self.play(others.animate.set_color(P.RED).set_opacity(0.25),
                  keep.animate.scale(2.2).set_color(P.GREEN), FadeOut(cap3), FadeIn(cap4),
                  run_time=1.6)
        timing.hold_to_read(self, cap4, settle=0.6)
        self.play(FadeOut(VGroup(dots, count, cap4)))

        # 3. how far away a 0.57 Jy source would be
        dist = np.array(d["distance_au"])
        plot = Plot([100, 450], [0.05, 3.0], [100, 150, 200, 250, 300, 350, 400, 450],
                    [0.1, 0.3, 1, 3], "distance from the Sun (AU)", "60 µm flux density (Jy)",
                    y_log=True, centre=(-0.4, 0.5), width=8.8, height=4.3)
        curves = VGroup()
        for k, cv in enumerate(d["curves"]):
            curves.add(plot.curve(dist, cv["flux_60um_jy"], P.BLUE,
                                  stroke_width=[1.8, 3.0, 1.8][k]))
        cl = layout.label(f"{d['published_mass_lo']:.0f}-{d['published_mass_hi']:.0f} Earth "
                          f"masses at {d['model_temp_k']:.0f} K", font_size=16, color=P.BLUE)
        cl.next_to(plot.p(330, 0.2), UP + RIGHT, buff=0.1)
        flux = plot.hline(c["flux_60um_jy"], P.GREEN)
        flux_lab = layout.label(f"candidate {c['flux_60um_jy']:.2f} Jy", font_size=16,
                                color=P.GREEN)
        flux_lab.next_to(plot.p(450, c["flux_60um_jy"]), UP + LEFT, buff=0.08)
        pub = plot.band(d["published_distance_au"] - d["published_distance_err_au"],
                        d["published_distance_au"] + d["published_distance_err_au"], 0.05, 3.0,
                        P.FG, opacity=0.1)
        pub_lab = layout.label(f"paper: {d['published_distance_au']:.0f} ± "
                               f"{d['published_distance_err_au']:.0f} AU", font_size=15,
                               color=P.FG)
        pub_lab.next_to(plot.p(d["published_distance_au"], 3.0), DOWN, buff=0.12)
        self.play(FadeIn(plot))
        cap5 = layout.caption("A body's heat fades as distance squared", font_size=22)
        self.play(Create(curves), FadeIn(cl), FadeIn(cap5), run_time=1.4)
        timing.hold_to_read(self, cap5, settle=0.2)
        hits = VGroup(*[Dot(plot.p(cv["implied_distance_au"], c["flux_60um_jy"]), radius=0.07,
                            color=P.GREEN) for cv in d["curves"]])
        lo = min(cv["implied_distance_au"] for cv in d["curves"])
        hi = max(cv["implied_distance_au"] for cv in d["curves"])
        box = _readout("implied distance", f"{lo:.0f}-{hi:.0f} AU",
                       f"paper: {d['published_distance_au']:.0f} ± "
                       f"{d['published_distance_err_au']:.0f} AU", P.GREEN)
        box.to_edge(RIGHT, buff=0.35).shift(UP * 1.2)
        temps = [cv["temperature_at_published_distance_k"] for cv in d["curves"]]
        cap6 = layout.caption(
            f"At {d['model_temp_k']:.0f} K it lies nearer than the paper's fit; "
            f"{min(temps):.0f}-{max(temps):.0f} K would put it at {d['published_distance_au']:.0f} AU",
            font_size=22)
        self.play(Create(flux), FadeIn(flux_lab), FadeIn(hits), FadeIn(pub), FadeIn(pub_lab),
                  FadeIn(box), FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, box, settle=0.6)
        self.play(FadeOut(VGroup(plot, curves, cl, flux, flux_lab, hits, pub, pub_lab, box,
                                 cap6)))

        # 4. where it is, against where Planet Nine is predicted to be
        m = sky.SkyMap(width=11.0, dec_range=(-80, 80), centre=(0.0, 0.5, 0.0))
        ecl, gal = m.reference_curves()
        pop = m.dots(d["population"], color=P.TEAL, radius=0.028)
        star = Star(n=5, outer_radius=0.17, color=P.GREEN).set_fill(P.GREEN, 1.0)
        star.move_to(m.p(c["ra_deg"], c["dec_deg"]))
        star_lab = layout.label(
            f"the IRAS candidate: {c['ecliptic_lat_deg']:.0f}° from the ecliptic", font_size=16,
            color=P.GREEN)
        star_lab.next_to(star, RIGHT, buff=0.15)
        self.play(FadeIn(m), Create(ecl), Create(gal), run_time=1.0)
        cap7 = layout.caption(
            f"Predicted Planet Nines stay within {d['population_max_ecliptic_lat_deg']:.0f}° "
            "of the ecliptic", font_size=22)
        self.play(FadeIn(pop, lag_ratio=0.02), FadeIn(cap7), run_time=1.2)
        timing.hold_to_read(self, cap7, settle=0.2)
        cap8 = layout.caption("The candidate is far off every predicted orbit, and unconfirmed",
                              font_size=22)
        self.play(FadeIn(star, scale=2.0), FadeIn(star_lab), FadeOut(cap7), FadeIn(cap8))
        timing.hold_to_read(self, cap8, settle=1.0)
        self.play(FadeOut(cap8))

        layout.show_takeaway(
            self, "One faint IRAS mover survives, but not where Planet Nine should be.")
