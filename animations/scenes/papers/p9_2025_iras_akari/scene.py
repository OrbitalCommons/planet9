"""Phan et al. (2025) -- a search for Planet Nine with IRAS and AKARI data.

Two far-infrared all-sky catalogues taken 23 years apart. A bound planet at
500-700 AU moves 42'-69.6' between them, so the search pairs an IRAS source
with an AKARI source that far away and with compatible fluxes. One pair
survives. The scene draws the crate's two-epoch motion model, the candidate's
real catalogue positions, and the far-infrared brightness the thermal model
predicts. Everything is from anim.json -> papers -> p9-2025-iras-akari.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Annulus,
    Arrow,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    Square,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2025-iras-akari"


def rows(items, font_size=16, buff=0.2):
    """Stacked (text, colour) lines, left aligned."""
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


class IrasAkari2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        cand = d["candidate"]
        pub = d["published"]
        win_lo, win_hi = d["window_arcmin"]
        d_lo, d_hi = d["window_distance_au"]

        self.add(paper.scene_header(CRATE))

        # 1. how far a bound planet moves in 23 years
        s = d["separation"]
        dist = np.array(s["distance_au"])
        ax, labels = widgets.labeled_axes(
            [350, 950, 100], [0, 120, 20], x_label="heliocentric distance (AU)",
            y_label="motion between the surveys (arcmin)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.3)
        top = [ax.c2p(x, min(y, 120)) for x, y in zip(dist, s["max_arcmin"])]
        bottom = [ax.c2p(x, y) for x, y in zip(dist, s["min_arcmin"])]
        band = Polygon(*top, *bottom[::-1], stroke_width=0).set_fill(P.BLUE, opacity=0.18)
        circ = widgets.curve(ax, dist, s["circular_arcmin"], color=P.BLUE)
        window = Polygon(ax.c2p(350, win_lo), ax.c2p(950, win_lo), ax.c2p(950, win_hi),
                         ax.c2p(350, win_hi), stroke_width=0).set_fill(P.PURPLE, opacity=0.2)
        key = rows([
            ("bound orbits, e ≤ 0.7, any season", P.BLUE),
            ("circular orbit", P.BLUE),
            (f"search window {win_lo:.0f}′ to {win_hi:.1f}′", P.PURPLE),
        ], font_size=15)
        key.move_to([4.6, 2.0, 0])
        cap = layout.caption(
            f"IRAS 1983, AKARI 2006: in {d['baseline_years']:.1f} years a distant planet "
            "drifts most of a degree", font_size=21)
        self.play(Create(ax), FadeIn(labels), run_time=1.0)
        self.play(FadeIn(band), Create(circ), FadeIn(key[0]), FadeIn(key[1]), FadeIn(cap),
                  run_time=1.2)
        timing.hold_to_read(self, cap, settle=0.5)
        edges = VGroup(
            DashedLine(ax.c2p(d_lo, 0), ax.c2p(d_lo, win_hi), color=P.PURPLE, stroke_width=2),
            DashedLine(ax.c2p(d_hi, 0), ax.c2p(d_hi, win_lo), color=P.PURPLE, stroke_width=2))
        note = rows([
            ("fastest orbits reach the window at", P.FG),
            (f"{d_lo:.0f} to {d_hi:.0f} AU", P.PURPLE),
            (f"paper: {pub['distance_au'][0]:.0f} to {pub['distance_au'][1]:.0f} AU", P.MUTED),
        ], font_size=15, buff=0.12)
        note.next_to(key, DOWN, buff=0.5, aligned_edge=LEFT)
        cap2 = layout.caption("The search keeps IRAS-AKARI source pairs separated by that much",
                              font_size=21)
        self.play(FadeIn(window), FadeIn(key[2]), Create(edges), FadeIn(note), FadeOut(cap),
                  FadeIn(cap2), run_time=1.2)
        timing.hold_to_read(self, cap2, note, settle=0.8)
        self.play(FadeOut(VGroup(ax, labels, band, circ, window, key, edges, note, cap2)))

        # 2. the pair that survives, at its catalogue positions
        centre = np.array([-3.4, 0.1, 0.0])
        half, span = 2.35, 80.0            # panel half-size (scene units), arcmin
        k = half / span
        frame = Square(side_length=2 * half, color=P.MUTED, stroke_width=1.4)
        frame.set_fill("#16171f", opacity=1.0).move_to(centre)
        ring = Annulus(inner_radius=win_lo * k, outer_radius=win_hi * k, color=P.PURPLE,
                       fill_opacity=0.16, stroke_width=0).move_to(centre)
        ring_edges = VGroup(
            Circle(radius=win_lo * k, color=P.PURPLE, stroke_width=1.2).move_to(centre),
            Circle(radius=win_hi * k, color=P.PURPLE, stroke_width=1.2).move_to(centre))
        p_iras = centre
        p_akari = centre + k * np.array([-cand["delta_ra_arcmin"], cand["delta_dec_arcmin"], 0])
        dot_i = Dot(p_iras, radius=0.08, color=P.GREEN)
        dot_a = Dot(p_akari, radius=0.08, color=P.ORANGE)
        lab_i = layout.label("IRAS 1983", font_size=14, color=P.GREEN).next_to(dot_i, UP, buff=0.1)
        lab_a = layout.label("AKARI 2006", font_size=14, color=P.ORANGE)
        lab_a.next_to(dot_a, DOWN, buff=0.1)
        move = Arrow(p_iras, p_akari, buff=0.1, color=P.FG, stroke_width=2.5,
                     max_tip_length_to_length_ratio=0.12)
        bar = Line(centre + [-half + 0.25, -half + 0.3, 0],
                   centre + [-half + 0.25 + 30 * k, -half + 0.3, 0], color=P.FG, stroke_width=2)
        bar_lab = layout.label("30′", font_size=13, color=P.FG).next_to(bar, UP, buff=0.06)
        compass = layout.label("north up, east left", font_size=12, color=P.MUTED)
        compass.move_to(centre + [half - 0.95, -half + 0.2, 0])
        panel = VGroup(frame, ring, ring_edges, bar, bar_lab, compass)

        m = sky.SkyMap(width=5.4, dec_range=(-90, 90), centre=(3.55, 1.45, 0.0))
        ecl, gal = m.reference_curves()
        here = Circle(radius=0.13, color=P.GREEN, stroke_width=2.5).move_to(
            m.p(cand["iras_ra_deg"], cand["iras_dec_deg"]))
        facts = rows([
            (f"separation {cand['separation_arcmin']:.1f}′   "
             f"(paper: {pub_sep(d):.1f}′)", P.FG),
            (f"{cand['rate_arcmin_yr']:.2f}′ per year, toward the south-west", P.FG),
            (f"{cand['distance_circular_au']:.0f} AU if circular, up to "
             f"{cand['distance_au']:.0f} AU if eccentric", P.BLUE),
        ], font_size=15, buff=0.16)
        facts.move_to([3.55, -1.75, 0])
        cap3 = layout.caption(
            f"{pub['pairs']} pairs pass the cuts; {pub['good']} survives image inspection",
            font_size=21)
        self.play(FadeIn(panel), FadeIn(m), Create(ecl), Create(gal), run_time=1.0)
        self.play(FadeIn(dot_i), FadeIn(lab_i), Create(here), FadeIn(cap3), run_time=0.8)
        self.play(Create(move), FadeIn(dot_a), FadeIn(lab_a), run_time=1.0)
        self.play(FadeIn(facts, lag_ratio=0.2))
        timing.hold_to_read(self, cap3, facts, settle=1.0)
        chance = layout.caption(
            f"Unrelated sources alone give about {d['chance_pairs']:.0f} pairs this far apart "
            f"({d['post_cut_iras']} × {d['post_cut_akari']} sources)", font_size=21)
        self.play(FadeOut(cap3), FadeIn(chance))
        timing.hold_to_read(self, chance, settle=0.8)
        self.play(FadeOut(VGroup(panel, m, ecl, gal, here, dot_i, dot_a, lab_i, lab_a, move,
                                 facts, chance)))

        # 3. is a planet bright enough to be in IRAS at all?
        f = d["flux"]
        fd = np.array(f["distance_au"])
        ax3, labels3 = widgets.labeled_axes(
            [300, 900, 100], [0, 0.8, 0.2], x_label="heliocentric distance (AU)",
            y_label="flux at 60 µm (Jy)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=-0.3)
        VGroup(ax3, labels3).shift(LEFT * 2.3)
        shades = {30.0: 0.45, 40.0: 0.7, 50.0: 1.0}
        curves, tags = VGroup(), VGroup()
        for c in f["curves"]:
            if c["mass_earth"] != 10.0:
                continue
            pts = [(x, y) for x, y in zip(fd, c["f60_jy"]) if y <= 0.8]
            line = widgets.curve(ax3, [p[0] for p in pts], [p[1] for p in pts], color=P.BLUE)
            line.set_stroke(opacity=shades[c["t_eff"]])
            curves.add(line)
            tag = layout.label(f"{c['t_eff']:.0f} K", font_size=15, color=P.BLUE)
            tag.set_opacity(shades[c["t_eff"]])
            tag.next_to(ax3.c2p(pts[-1][0], pts[-1][1]), RIGHT, buff=0.12)
            tags.add(tag)
        # the cold curves end close together: keep their tags at least a line apart
        stack = sorted(tags, key=lambda t: t.get_y())
        for lower, upper in zip(stack, stack[1:]):
            if upper.get_y() - lower.get_y() < 0.22:
                upper.set_y(lower.get_y() + 0.22)
        names = [("10 M⊕ planet, three temperatures", P.BLUE)]
        lim = d["iras"]["limit_jy"]
        limit = DashedLine(ax3.c2p(300, lim), ax3.c2p(900, lim), color=P.PURPLE, stroke_width=2)
        seen = Line(ax3.c2p(cand["distance_circular_au"], cand["iras_flux_60"]),
                    ax3.c2p(cand["distance_au"], cand["iras_flux_60"]), color=P.GREEN,
                    stroke_width=6)
        key3 = rows(names + [
            (f"IRAS limit {lim:.1f} Jy", P.PURPLE),
            (f"candidate: {cand['iras_flux_60']:.2f} Jy", P.GREEN),
        ], font_size=15)
        key3.move_to([4.6, 1.6, 0])
        cap4 = layout.caption(
            "Far-infrared light is the planet's own heat: only a body near 50 K reaches IRAS",
            font_size=21)
        self.play(Create(ax3), FadeIn(labels3), run_time=0.9)
        self.play(*[Create(c) for c in curves], Create(limit), FadeIn(key3[:2]), FadeIn(cap4),
                  run_time=1.4)
        self.play(FadeIn(tags))
        self.play(Create(seen), FadeIn(key3[2]))
        timing.hold_to_read(self, cap4, key3, settle=1.0)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, "One pair fits a planet at 500-700 AU, but two positions cannot confirm it.")


def pub_sep(d):
    """The separation printed in the paper's Table 2, from the ledger entry."""
    from p9_manim import ledger

    return float(ledger.entry(CRATE)["result"]["published"])
