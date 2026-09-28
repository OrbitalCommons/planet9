"""Finale -- where a small space telescope should image.

WhereToImage walks the sky: the probability nobody else will collect, the three
regions the reference campaign images and why each is there. WhatItBuys shows
the return on telescope time under each planning stance and the campaign
calendar. Everything is read from figures/space_strategy.json, written by
``cargo run -p p9-space-strategy``.
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
    Line,
    Polygon,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, sky, timing, widgets

ZONE_COLOURS = {
    "Anticentre crossing": P.ORANGE,
    "North of Rubin": P.TEAL,
    "Galactic-centre crossing": P.RED,
}
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def _season(months):
    """'Jul Jun' -> 'Jun–Jul': the span of a zone's opposition months, read
    around the year starting in September."""
    idx = sorted((MONTHS.index(m) for m in months.split()), key=lambda k: (k - 8) % 12)
    names = [MONTHS[k] for k in idx]
    return names[0] if len(names) == 1 else f"{names[0]}–{names[-1]}"


def _header(text):
    badge = layout.concept_badge("FINALE")
    title = layout.label(text, font_size=24, color=P.FG, weight="BOLD")
    title.to_edge(UP, buff=0.32)
    return VGroup(badge, title)


def _tile(m, t, colour, opacity=0.0, stroke=1.4):
    a = m.p(t["ra_deg"] - t["ra_width_deg"] / 2, t["dec_deg"] - t["dec_height_deg"] / 2)
    b = m.p(t["ra_deg"] + t["ra_width_deg"] / 2, t["dec_deg"] + t["dec_height_deg"] / 2)
    poly = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], color=colour, stroke_width=stroke)
    return poly.set_fill(colour, opacity=opacity)


class WhereToImage(Scene):
    def construct(self):
        d = dataio.space_strategy()
        self.add(_header("Where a small space telescope should image"))

        m = sky.SkyMap(width=12.2, dec_range=(-60, 60), centre=(0.0, 0.05, 0.0))
        ecl, gal = m.reference_curves()
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), run_time=1.0)

        # 1. what is left for anyone
        cells = d["sky"]
        dens = np.array([c["unique"] / (c["ra_width_deg"] * c["dec_height_deg"]) for c in cells])
        peak = dens.max()
        shade = VGroup()
        for c, f in zip(cells, np.sqrt(dens / peak)):
            if f >= 0.12 and m.inside(c["dec_deg"]):
                shade.add(_tile(m, c, P.ORANGE, opacity=0.75 * f, stroke=0))
        stance = d["policies"][2]
        cap = layout.caption(
            f"Ground surveys have found {100 * stance['found_by_ground']:.0f}% of the prediction. "
            f"Shaded: the {100 * stance['unique']:.0f}% nobody else will search.",
            font_size=20)
        self.play(FadeIn(shade, lag_ratio=0.002), FadeIn(cap), run_time=1.8)
        timing.hold_to_read(self, cap, settle=0.8)

        # 2. Rubin's reach
        rubin = Line(m.p(180.0 + 179.9, 12.0), m.p(180.0 - 179.9, 12.0),
                     color=P.BLUE, stroke_width=2.2)
        rubin_lab = layout.label("Rubin takes everything south of +12°", font_size=14,
                                 color=P.BLUE)
        rubin_lab.next_to(m.p(308, 12), UP, buff=0.08)
        cap2 = layout.caption(
            "Rubin will reach V ≈ 25 here. Planet Nine is V ≈ 19–23: leave that sky to Rubin.",
            font_size=20)
        self.play(Create(rubin), FadeIn(rubin_lab), FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, settle=0.6)
        self.play(FadeOut(cap2))

        # 3. the three regions, one at a time
        for zone in d["zones"]:
            colour = ZONE_COLOURS[zone["zone"]]
            tiles = VGroup(*[_tile(m, t, colour, opacity=0.18)
                             for t in d["tiles"] if t["zone"] == zone["zone"]])
            ra0, ra1 = zone["ra_range_h"]
            name = layout.label(zone["zone"], font_size=17, color=colour, weight="BOLD")
            where = layout.label(
                f"RA {ra0:.1f}h–{ra1:.1f}h   Dec {zone['dec_range_deg'][0]:+.0f}° to "
                f"{zone['dec_range_deg'][1]:+.0f}°   {zone['area_deg2']:.0f} deg²",
                font_size=13, color=P.FG)
            card = VGroup(name, where).arrange(DOWN, buff=0.07, aligned_edge=LEFT)
            anchor = tiles.get_center()
            card.next_to(tiles, UP if anchor[1] < 1.0 else DOWN, buff=0.18)
            if card.get_right()[0] > 6.6:
                card.shift(LEFT * (card.get_right()[0] - 6.6))
            if card.get_left()[0] < -6.6:
                card.shift(RIGHT * (-6.6 - card.get_left()[0]))
            line = layout.caption(
                f"{zone['why'][0].upper()}{zone['why'][1:]}: "
                f"{zone['median_integration_s']:.0f} s visits, "
                f"{zone['hours']:.0f} h, best in {_season(zone['opposition_months'])}.",
                font_size=19)
            self.play(FadeIn(tiles, lag_ratio=0.01), FadeIn(card), FadeIn(line), run_time=1.3)
            timing.hold_to_read(self, line, settle=1.2)
            self.play(FadeOut(line), FadeOut(card), tiles.animate.set_fill(opacity=0.06),
                      run_time=0.6)

        layout.show_takeaway(
            self,
            f"Three regions, {d['reference_area_deg2']:.0f} deg², "
            f"{d['reference_hours']:.0f} hours: sky only a sharp telescope can search.")


class WhatItBuys(Scene):
    def construct(self):
        d = dataio.space_strategy()
        self.add(_header("What telescope time buys"))

        # 1. return against hours, by planning stance
        budgets = d["budgets_h"]
        top = max(max(p["captured"]) for p in d["policies"]) * 100
        ymax = int(np.ceil(top / 5.0) * 5)
        ax, labels = widgets.labeled_axes(
            [0, 8000, 1000], [0, ymax, 5], x_label="wall-clock hours",
            y_label="chance of finding Planet Nine (%)", numbers=True, x_length=8.2, y_length=4.4, shift_down=-0.25)
        VGroup(ax, labels).shift(LEFT * 2.2)
        self.play(Create(ax), FadeIn(labels))
        colours = [P.RED, P.ORANGE, P.TEAL]
        legend = VGroup()
        for pol, col in zip(d["policies"][:3], colours):
            xs = [0.0] + list(budgets)
            ys = [0.0] + [100 * c for c in pol["captured"]]
            curve = widgets.curve(ax, xs, ys, color=col, stroke_width=3.2)
            row = VGroup(
                layout.label(pol["policy"], font_size=15, color=col, weight="BOLD"),
                layout.label(f"{100 * pol['unique']:.0f}% of the prediction left to find",
                             font_size=12, color=P.MUTED),
            ).arrange(DOWN, buff=0.05, aligned_edge=LEFT)
            legend.add(row)
            legend.arrange(DOWN, buff=0.28, aligned_edge=LEFT)
            legend.move_to([4.6, 1.2, 0])
            self.play(Create(curve), FadeIn(row), run_time=1.1)
        ref = Dot(ax.c2p(d["reference_hours"], 100 * d["reference_captured"]),
                  radius=0.09, color=P.FG)
        ref_lab = layout.label(
            f"reference campaign: {d['reference_hours']:.0f} h, "
            f"{100 * d['reference_captured']:.1f}%", font_size=14, color=P.FG)
        ref_lab.next_to(ref, DOWN + RIGHT, buff=0.1)
        cap = layout.caption(
            "Most of what is left is Rubin's. The campaign takes the part that is not.",
            font_size=20)
        self.play(FadeIn(ref), FadeIn(ref_lab), FadeIn(cap))
        timing.hold_to_read(self, cap, legend, settle=1.2)
        self.play(*[FadeOut(x) for x in self.mobjects[1:]])

        # 2. the calendar
        hours = d["hours_by_month"]
        order = [(k + 6) % 12 for k in range(12)]
        by_zone = {z: [0.0] * 12 for z in ZONE_COLOURS}
        for t in d["tiles"]:
            by_zone[t["zone"]][t["opposition_month"] - 1] += t["hours"]
        ymax = int(np.ceil(max(hours) / 100.0) * 100)
        ax2, lab2 = widgets.labeled_axes(
            [0, 12, 1], [0, ymax, 100], y_label="telescope hours",
            x_length=9.6, y_length=3.9, shift_down=-0.1)
        ax2.get_y_axis().add_numbers()
        self.play(Create(ax2), FadeIn(lab2))
        bars, ticks = VGroup(), VGroup()
        for slot, month in enumerate(order):
            base = 0.0
            for zone, col in ZONE_COLOURS.items():
                h = by_zone[zone][month]
                if h > 0:
                    a, b = ax2.c2p(slot + 0.15, base), ax2.c2p(slot + 0.85, base + h)
                    bar = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0)
                    bars.add(bar.set_fill(col, opacity=0.85))
                    base += h
            ticks.add(layout.label(MONTHS[month], font_size=13, color=P.FG)
                      .next_to(ax2.c2p(slot + 0.5, 0), DOWN, buff=0.12))
        key = VGroup(*[layout.label(z, font_size=14, color=c, weight="BOLD")
                       for z, c in ZONE_COLOURS.items()]).arrange(DOWN, buff=0.12,
                                                                   aligned_edge=LEFT)
        key.move_to(ax2.c2p(2.6, 0.78 * ymax))
        cap2 = layout.caption(
            "Each region is imaged at opposition, when the planet moves fastest.",
            font_size=20)
        self.play(FadeIn(ticks), FadeIn(bars, lag_ratio=0.05), FadeIn(key), FadeIn(cap2),
                  run_time=1.6)
        timing.hold_to_read(self, cap2, settle=1.4)
        self.play(FadeOut(cap2))

        layout.show_takeaway(
            self, "Four short visits per field, hours apart: the motion is the detection.")
