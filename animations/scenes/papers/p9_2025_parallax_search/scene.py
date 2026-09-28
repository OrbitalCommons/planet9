"""Socas-Navarro & Trujillo (2025) -- a targeted, parallax-based search for
Planet Nine.

Two images of the same 98 square-degree field on consecutive nights. Earth
moves 0.017 AU in a day, so a planet several hundred AU away shifts by a few
arcseconds against the stars, 20-30 times more than its own orbital drift. The
crate computes that shift; the field size, depth and seeing are the paper's.
Everything is from anim.json -> papers -> p9-2025-parallax-search.
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
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Rectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-parallax-search"


def rows(items, font_size=15, buff=0.16):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


class ParallaxSearch2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        lo, hi = pub["nightly_arcsec"]
        d_near, d_far = d["distance_for_shift"]

        self.add(paper.scene_header(CRATE))

        # 1. one night of motion against distance
        dist = np.array(d["distance_au"])
        ax, labels = widgets.labeled_axes(
            [300, 1000, 100], [0, 12, 2], x_label="heliocentric distance (AU)",
            y_label="shift in one night (arcsec)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.3)
        par = widgets.curve(ax, dist, d["net_arcsec"], color=P.ORANGE)
        orb = widgets.curve(ax, dist, d["orbital_arcsec"], color=P.BLUE)
        band = Polygon(ax.c2p(300, lo), ax.c2p(1000, lo), ax.c2p(1000, hi), ax.c2p(300, hi),
                       stroke_width=0).set_fill(P.PURPLE, opacity=0.2)
        edges = VGroup(
            DashedLine(ax.c2p(d_near, 0), ax.c2p(d_near, hi), color=P.PURPLE, stroke_width=2),
            DashedLine(ax.c2p(d_far, 0), ax.c2p(d_far, lo), color=P.PURPLE, stroke_width=2))
        r_lo, r_hi = pub["ratio"]
        key = rows([
            ("apparent shift, seen from Earth", P.ORANGE),
            (f"{d['shift_500']:.1f}″ at 500 AU, {d['shift_700']:.1f}″ at 700 AU", P.ORANGE),
            ("the planet's own orbital drift", P.BLUE),
            (f"{d['ratio_500']:.0f} times smaller at 500 AU", P.BLUE),
            (f"paper: {r_lo:.0f} to {r_hi:.0f} times", P.MUTED),
        ])
        key[2:].shift(DOWN * 0.25)
        key.move_to([4.6, 1.6, 0])
        cap = layout.caption(
            f"In one night Earth moves {d['baseline_au']:.3f} AU: nearby things shift, "
            "the stars do not", font_size=21)
        self.play(Create(ax), FadeIn(labels), run_time=1.0)
        self.play(Create(par), FadeIn(key[:2]), FadeIn(cap), run_time=1.2)
        self.play(Create(orb), FadeIn(key[2:]), run_time=1.0)
        timing.hold_to_read(self, cap, key, settle=0.6)
        note = rows([
            (f"paper: {lo:.0f}″ to {hi:.0f}″ per night", P.PURPLE),
            (f"here that is {d_near:.0f} to {d_far:.0f} AU", P.PURPLE),
        ])
        note.next_to(key, DOWN, buff=0.45, aligned_edge=LEFT)
        cap_b = layout.caption("The search wants a source displaced that far, along the "
                               "direction Earth moved", font_size=21)
        self.play(FadeIn(band), Create(edges), FadeIn(note), FadeOut(cap), FadeIn(cap_b))
        timing.hold_to_read(self, cap_b, note, settle=0.8)
        self.play(FadeOut(VGroup(ax, labels, par, orb, band, edges, key, note, cap_b)))

        # 2. the two nights, to scale
        scale = 0.62                      # scene units per arcsec
        seeing = pub["seeing_arcsec"]
        field = Rectangle(width=12.4, height=4.6, color=P.MUTED, stroke_width=1.2)
        field.set_fill("#16171f", opacity=1.0).move_to([0, 0.2, 0])
        pairs = VGroup()
        for au, shift, y in ((500, d["shift_500"], 1.2), (700, d["shift_700"], -0.8)):
            a = np.array([-3.4, y, 0.0])
            b = a + RIGHT * shift * scale
            first = Circle(radius=0.5 * seeing * scale, color=P.BLUE, stroke_width=2)
            first.set_fill(P.BLUE, opacity=0.25).move_to(a)
            second = Circle(radius=0.5 * seeing * scale, color=P.BLUE, stroke_width=2)
            second.set_fill(P.BLUE, opacity=0.8).move_to(b)
            hop = Arrow(a, b, buff=0.35, color=P.ORANGE, stroke_width=3,
                        max_tip_length_to_length_ratio=0.1)
            tag = layout.label(f"{shift:.1f}″", font_size=16, color=P.ORANGE, weight="BOLD")
            tag.next_to(hop, UP, buff=0.08)
            who = layout.label(f"planet at {au} AU", font_size=15, color=P.FG)
            who.next_to(first, LEFT, buff=0.5)
            pairs.add(VGroup(first, second, hop, tag, who))
        n1 = layout.label("night 1", font_size=13, color=P.MUTED)
        n1.next_to(pairs[0][0], UP, buff=0.35)
        n2 = layout.label("night 2", font_size=13, color=P.MUTED)
        n2.next_to(pairs[0][1], UP, buff=0.35)
        star = Dot([3.6, 0.2, 0], radius=0.5 * seeing * scale, color=P.FG).set_opacity(0.8)
        star_lab = layout.label(
            timing.wrap("a star: same place both nights, image 1″ wide in good seeing",
                        width=26), font_size=14, color=P.FG, line_spacing=0.9)
        star_lab.set_opacity(0.75)
        star_lab.next_to(star, DOWN, buff=0.3)
        bar = Line([3.6 - 2.5 * scale, -1.75, 0], [3.6 + 2.5 * scale, -1.75, 0], color=P.FG,
                   stroke_width=2)
        bar_lab = layout.label("5″", font_size=13, color=P.FG).next_to(bar, DOWN, buff=0.08)
        err = d["tolerable_error_arcsec"]
        cap2 = layout.caption(
            f"The shift is several image widths: positions good to {err['at_700']:.1f}″ "
            f"give {err['snr']:.0f}σ even at 700 AU", font_size=21)
        self.play(FadeIn(field), FadeIn(star), FadeIn(star_lab), Create(bar), FadeIn(bar_lab),
                  run_time=0.9)
        for pair in pairs:
            self.play(FadeIn(pair[0]), FadeIn(pair[4]), run_time=0.5)
            self.play(Create(pair[2]), FadeIn(pair[1]), FadeIn(pair[3]), run_time=0.9)
        self.play(FadeIn(n1), FadeIn(n2), FadeIn(cap2))
        timing.hold_to_read(self, cap2, star_lab, settle=0.8)
        self.play(FadeOut(VGroup(field, pairs, n1, n2, star, star_lab, bar, bar_lab, cap2)))

        # 3. how deep the two nights reach
        r = np.array(d["r_mags"])
        edges_h = np.arange(17.0, 26.01, 0.5)
        counts, _ = np.histogram(r, bins=edges_h)
        depth = pub["depth_r"]
        top = int(np.ceil(counts.max() / 100.0) * 100)
        ax3, labels3 = widgets.labeled_axes(
            [17, 26, 1], [0, top, top // 4], x_label="predicted r magnitude  (fainter →)",
            y_label="synthetic planets", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.0, shift_down=-0.3)
        VGroup(ax3, labels3).shift(LEFT * 2.3)
        centres = 0.5 * (edges_h[:-1] + edges_h[1:])
        seen = np.where(centres < depth, counts, 0)
        bars_seen = widgets.histogram(ax3, edges_h, seen, color=P.PURPLE, opacity=0.7)
        bars_rest = widgets.histogram(ax3, edges_h, counts - seen, color=P.TEAL, opacity=0.55)
        mark = widgets.marker_line(ax3, depth, (0, top), f"depth reached  r = {depth:.1f}",
                                   side=UP)
        w, h = pub["field_deg"]
        key3 = rows([
            (f"{100 * d['bright_fraction']:.0f}% of the 2021 prediction", P.PURPLE),
            ("is brighter than the limit", P.PURPLE),
            ("fainter: out of reach", P.TEAL),
            (f"field: {w:.1f}° by {h:.1f}°, {pub['area_deg2']:.0f} square degrees", P.FG),
            (f"{100 * d['sky_fraction']:.2f}% of the sky", P.FG),
        ])
        key3[3:].shift(DOWN * 0.25)
        key3.move_to([4.6, 1.3, 0])
        cap3 = layout.caption(
            f"No candidate: nothing brighter than r = {depth:.1f} in the field moves "
            "like Planet Nine", font_size=21)
        self.play(Create(ax3), FadeIn(labels3), run_time=0.9)
        self.play(FadeIn(bars_seen, lag_ratio=0.1), FadeIn(bars_rest, lag_ratio=0.1),
                  Create(mark), FadeIn(key3), FadeIn(cap3), run_time=1.5)
        timing.hold_to_read(self, cap3, key3, settle=1.0)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, "Two nights suffice, but only for a bright planet in a small field.")
