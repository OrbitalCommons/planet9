"""Belyakov, Bernardinelli & Brown (2022) -- limits on the detection of Planet
Nine in the Dark Energy Survey.

DES imaged 5,000 deg² of the southern sky ten times over six years to r = 23.8,
three magnitudes deeper than ZTF. Deep enough for nearly every predicted
Planet Nine -- but only the few whose paths cross its footprint can be tested.
Reproduced in p9-2022-des: the footprint, the population scored by both the
DES and ZTF survey models, the recovery rates and the exclusion bookkeeping
are the crate's own (anim.json -> papers -> p9-2022-des).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    FadeIn,
    FadeOut,
    Polygon,
    Scene,
    SurroundingRectangle,
    VGroup,
)

import p9_manim as P
from p9_manim.plot import Plot
from p9_manim import dataio, layout, paper, sky, timing

CRATE = "p9-2022-des"


def _readout(title, value, note, colour):
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


def _bars(plot, edges, counts, base, colour, opacity):
    g = VGroup()
    for k, n in enumerate(counts):
        if n <= 0:
            continue
        a = plot.p(edges[k], base[k])
        b = plot.p(edges[k + 1], base[k] + n)
        bar = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=colour, stroke_width=0.6)
        g.add(bar.set_fill(colour, opacity=opacity))
    return g


class Des2022(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pop = d["population"]
        ztf = [s for s in pop if s["p_ztf"] >= 0.5]
        new = [s for s in pop if s["p_ztf"] < 0.5 and s["p_des"] >= 0.5]
        rest = [s for s in pop if s["p_ztf"] < 0.5 and s["p_des"] < 0.5]
        self.add(paper.scene_header(CRATE))

        # 1. the predicted planets, and what ZTF already removed
        m = sky.SkyMap(width=9.6, dec_range=(-80, 80), centre=(-1.55, 0.55, 0.0), ra_centre=0.0)
        ecl, gal = m.reference_curves()
        dz = m.dots(ztf, color=P.TEAL, radius=0.028)
        dn = m.dots(new, color=P.TEAL, radius=0.028)
        dr = m.dots(rest, color=P.TEAL, radius=0.028)
        self.play(FadeIn(m), Create(ecl), Create(gal), run_time=1.0)
        cap = layout.caption(f"{len(pop)} predicted Planet Nines on tonight's sky", font_size=22)
        self.play(FadeIn(VGroup(dz, dn, dr), lag_ratio=0.02), FadeIn(cap), run_time=1.3)
        timing.hold_to_read(self, cap, settle=0.2)
        cap2 = layout.caption(
            f"ZTF had already ruled out the bright northern ones: {100 * d['ztf']:.0f}% of them",
            font_size=22)
        self.play(dz.animate.set_color(P.RED).set_opacity(0.3), FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.3)

        # 2. the DES footprint
        foot = VGroup(*[m.box(b["ra_lo"], b["ra_hi"], b["dec_lo"], b["dec_hi"], opacity=0.3,
                              stroke_width=0) for b in d["footprint"]])
        foot_lab = _readout("DES wide survey", f"{d['footprint_area_deg2']:,.0f} deg²",
                            f"{100 * d['footprint_sky_fraction']:.0f}% of the sky, r ≈ {d['depth_r']:.1f}",
                            P.PURPLE)
        foot_lab.next_to(m.frame, RIGHT, buff=0.25).align_to(m.frame, UP)
        cap3 = layout.caption(
            f"DES goes to r = {d['depth_r']:.1f}, three magnitudes deeper, but only in the south",
            font_size=22)
        self.play(FadeIn(foot), FadeIn(foot_lab), FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=0.3)
        cross = _readout("orbits crossing it", f"{100 * d['crossing_fraction']:.0f}%",
                         f"paper: {100 * d['published_crossing_fraction']:.0f}%", P.TEAL)
        cross.next_to(foot_lab, DOWN, buff=0.25)
        cap4 = layout.caption(
            "Those inside that DES would have linked are newly ruled out", font_size=22)
        self.play(*[x.animate.set_color(P.RED).scale(1.8) for x in dn], FadeIn(cross), FadeOut(cap3),
                  FadeIn(cap4), run_time=1.4)
        timing.hold_to_read(self, cap4, settle=0.8)
        self.play(FadeOut(VGroup(m, ecl, gal, dz, dn, dr, foot, foot_lab, cross, cap4)))

        # 3. deep enough for nearly all of them: the footprint is the limit
        r = np.array([s["r_mag"] for s in pop])
        is_z = np.array([s["p_ztf"] >= 0.5 for s in pop])
        is_n = np.array([s["p_ztf"] < 0.5 and s["p_des"] >= 0.5 for s in pop])
        edges = np.arange(16.0, 24.51, 0.5)
        nz, _ = np.histogram(r[is_z], bins=edges)
        nn, _ = np.histogram(r[is_n], bins=edges)
        na, _ = np.histogram(r, bins=edges)
        top = int(np.ceil(na.max() / 25.0) * 25)
        plot = Plot([16, 24.5], [0, top], [16, 17, 18, 19, 20, 21, 22, 23, 24],
                    list(range(0, top + 1, top // 5)), "brightness r  (fainter →)",
                    "predicted planets", centre=(-0.5, 0.45), width=9.2, height=4.2)
        b_z = _bars(plot, edges, nz, np.zeros_like(nz), P.RED, 0.35)
        b_n = _bars(plot, edges, nn, nz, P.RED, 0.9)
        b_r = _bars(plot, edges, na - nz - nn, nz + nn, P.TEAL, 0.5)
        comp = d["completeness"]
        c_line = plot.curve(comp["r_mag"], np.array(comp["fraction"]) * top, P.PURPLE,
                            stroke_width=2.5)
        c_lab = layout.label("DES detection efficiency", font_size=15, color=P.PURPLE)
        c_lab.next_to(plot.p(21.2, 0.97 * top), UP, buff=0.05)
        self.play(FadeIn(plot))
        cap5 = layout.caption("Every predicted planet by brightness", font_size=22)
        self.play(FadeIn(VGroup(b_z, b_n, b_r), lag_ratio=0.05), FadeIn(cap5), run_time=1.3)
        timing.hold_to_read(self, cap5, settle=0.2)
        cap6 = layout.caption(
            "DES could see nearly all of them; most just never cross its patch of sky",
            font_size=22)
        self.play(Create(c_line), FadeIn(c_lab), FadeOut(cap5), FadeIn(cap6))
        key = VGroup(
            layout.label("ruled out by ZTF", font_size=16, color=P.RED).set_opacity(0.6),
            layout.label("newly ruled out by DES", font_size=16, color=P.RED),
            layout.label("still possible", font_size=16, color=P.TEAL),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        rec = _readout("recovered if it crosses", f"{100 * d['recovery']:.0f}%",
                       f"paper: {100 * d['published_recovery']:.0f}%", P.PURPLE)
        uniq = _readout("new exclusion from DES", f"+{100 * d['des_unique']:.1f}%",
                        f"paper: +{100 * d['published_des_unique']:.0f}%", P.RED)
        side = VGroup(key, rec, uniq).arrange(DOWN, buff=0.3)
        side.to_edge(RIGHT, buff=0.3).align_to(plot.p(16, top), UP)
        self.play(FadeIn(key), FadeIn(rec))
        timing.hold_to_read(self, cap6, key, settle=0.4)
        cap7 = layout.caption(
            f"Together with ZTF, {100 * d['cumulative']:.0f}% of the predicted orbits are now gone",
            font_size=22)
        self.play(FadeIn(uniq), FadeOut(cap6), FadeIn(cap7))
        timing.hold_to_read(self, cap7, uniq, settle=1.0)
        self.play(FadeOut(cap7))

        layout.show_takeaway(
            self, "Deep but narrow: DES adds only a few percent to what ZTF excluded.")
