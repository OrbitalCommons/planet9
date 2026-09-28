"""Holman et al. (2025) -- a Pan-STARRS search for distant planets, part 1.

A second, deeper pass through Pan-STARRS1: every exposure, linked into orbits
out to 1600 AU, with the detection efficiency measured by injecting synthetic
detections into the source catalogues. No planet. The scene scores the same
Brown & Batygin (2021) population used for ZTF, DES and the 2024 Pan-STARRS1
search through the crate's survey model (north of -30 deg, 50% complete at
r = 22.5) and shows which members only this search reaches.
Everything is from anim.json -> papers -> p9-2025-ps1-holman.
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
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2025-ps1-holman"


def dot_key(items, font_size=16):
    row = VGroup()
    for text, col in items:
        row.add(VGroup(Dot(radius=0.06, color=col),
                       layout.label(text, font_size=font_size, color=col)).arrange(RIGHT, buff=0.1))
    return row.arrange(RIGHT, buff=0.5)


def rows(items, font_size=15, buff=0.16):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


class Ps1Holman2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pop = d["population"]
        pub = d["published"]
        sv = d["survey"]

        self.add(paper.scene_header(CRATE))

        # 1. what three searches had left
        m = sky.SkyMap(width=11.0, dec_range=(-75, 75), centre=(0.0, 0.42, 0.0))
        ecl, gal = m.reference_curves()
        old = [s for s in pop if s["p_before"] >= 0.5]
        new = [s for s in pop if s["p_before"] < 0.5 and s["p_after"] >= 0.5]
        left = [s for s in pop if s["p_after"] < 0.5]
        dots_old = m.dots(old, color=P.MUTED, opacity=0.45)
        dots_new = m.dots(new, color=P.TEAL)
        dots_left = m.dots(left, color=P.TEAL)
        key1 = dot_key([(f"ruled out by ZTF, DES, Pan-STARRS1 2024  {100 * d['before']:.0f}%",
                         P.MUTED),
                        (f"still viable  {100 * (1 - d['before']):.0f}%", P.TEAL)])
        key1.next_to(m.frame, DOWN, buff=0.62)
        cap = layout.caption(
            f"The predicted Planet Nines that three searches had not reached "
            f"({len(pop)} drawn)", font_size=21)
        self.play(FadeIn(m), run_time=0.7)
        self.play(Create(ecl), Create(gal), run_time=1.0)
        self.play(FadeIn(dots_old), FadeIn(dots_new), FadeIn(dots_left), FadeIn(key1),
                  FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, key1, settle=0.6)

        # 2. the same telescope, every exposure, a magnitude deeper
        foot = m.dec_band(sv["dec_limit_deg"], 90)
        cap2 = layout.caption(
            f"Pan-STARRS1 again: all {pub['exposures']:,} exposures, linked to "
            f"{pub['distance_au'][1]:.0f} AU, to w ≈ {pub['depth_w']:.1f}",
            font_size=21)
        self.play(FadeIn(foot), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        timing.hold_to_read(self, cap2, settle=0.6)

        key2 = dot_key([(f"earlier searches  {100 * d['before']:.0f}%", P.MUTED),
                        (f"new with this search  +{100 * d['new']:.1f}%", P.RED),
                        (f"still viable  {100 * d['remaining']:.1f}%", P.TEAL)])
        key2.next_to(m.frame, DOWN, buff=0.62)
        cap3 = layout.caption(
            f"No planet. Alone it would have found {100 * d['alone']:.0f}% of them "
            f"(paper: {100 * pub['alone']:.0f}%)", font_size=21)
        self.play(dots_new.animate.set_color(P.RED).set_opacity(1.0), FadeOut(key1),
                  FadeIn(key2), FadeOut(cap2), FadeIn(cap3), run_time=1.5)
        timing.hold_to_read(self, cap3, key2, settle=0.9)

        # 3. where the survivors are
        cap4 = layout.caption(
            f"Survivors crowd the far south and the galactic plane "
            f"(±{d['plane_b_deg']:.0f}°: {100 * d['plane_share_of_survivors']:.0f}%, "
            f"up from {100 * d['plane_share_of_population']:.0f}%)",
            font_size=21)
        self.play(dots_old.animate.set_opacity(0.12), dots_new.animate.set_opacity(0.15),
                  FadeOut(foot), gal.animate.set_stroke(opacity=1.0), FadeOut(cap3),
                  FadeIn(cap4), run_time=1.4)
        timing.hold_to_read(self, cap4, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, dots_old, dots_new, dots_left, key2, cap4)))

        # 4. why: depth
        r = np.array([s["r_mag"] for s in pop])
        w_old = np.array([s["p_before"] for s in pop])
        w_all = np.array([s["p_after"] for s in pop])
        edges = np.arange(17.0, 26.01, 0.5)
        h_old, _ = np.histogram(r, bins=edges, weights=w_old)
        h_all, _ = np.histogram(r, bins=edges, weights=w_all)
        tot, _ = np.histogram(r, bins=edges)
        top = max(20, int(np.ceil(tot.max() / 20.0) * 20))
        ax, labels = widgets.labeled_axes(
            [17, 26, 1], [0, top, top // 4], x_label="apparent magnitude r  (fainter →)",
            y_label="synthetic planets", y_rotate=True, numbers=True,
            x_length=7.0, y_length=4.0, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.3)
        bars_old = widgets.histogram(ax, edges, h_old, color=P.MUTED, opacity=0.6)
        bars_new = widgets.histogram(ax, edges, h_all - h_old, color=P.RED, opacity=0.8,
                                     base=h_old)
        bars_left = widgets.histogram(ax, edges, tot - h_all, color=P.TEAL, opacity=0.55,
                                      base=h_all)
        depth = widgets.marker_line(ax, sv["effective_depth"], (0, top),
                                    f"half complete at r = {sv['effective_depth']:.1f}",
                                    side=UP)
        found = rows([
            ("found along the way (paper)", P.FG),
            (f"{pub['objects']} solar-system objects", P.GREEN),
            (f"{pub['tnos']} beyond Neptune", P.GREEN),
            (f"{pub['dwarf_planets']} dwarf planets", P.GREEN),
            (f"{pub['new']} previously unknown", P.GREEN),
        ])
        legend = rows([
            ("ruled out by earlier searches", P.MUTED),
            ("newly ruled out", P.RED),
            ("still viable", P.TEAL),
        ])
        side = VGroup(legend, found).arrange(DOWN, buff=0.5, aligned_edge=LEFT)
        side.move_to([4.5, 0.9, 0])
        cap5 = layout.caption(
            "A magnitude deeper than the 2024 search of the same survey", font_size=21)
        self.play(Create(ax), FadeIn(labels), run_time=0.9)
        self.play(FadeIn(bars_old, lag_ratio=0.1), FadeIn(bars_left, lag_ratio=0.1),
                  FadeIn(bars_new, lag_ratio=0.1), Create(depth), FadeIn(legend),
                  FadeIn(cap5), run_time=1.5)
        self.play(FadeIn(found, lag_ratio=0.2))
        timing.hold_to_read(self, cap5, legend, found, settle=0.9)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, f"Four searches rule out {100 * d['cumulative']:.0f}%; the rest is "
                  "in the galactic plane or far south.")
