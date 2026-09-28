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
    "Galactic-centre crossing": P.PURPLE,
}
# where each zone's name card sits relative to its tiles
CARD_SIDE = {
    "Anticentre crossing": UP,
    "North of Rubin": UP,
    "Galactic-centre crossing": DOWN,
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


def _tick_numbers(ax, xs=(), ys=(), x_fmt="{:,.0f}", y_fmt="{:.0f}", font_size=14):
    """Legible tick numbers (the Axes defaults render too small at 480p)."""
    x0, y0 = ax.x_range[0], ax.y_range[0]
    g = VGroup()
    for x in xs:
        g.add(layout.label(x_fmt.format(x), font_size=font_size, color=P.FG)
              .next_to(ax.c2p(x, y0), DOWN, buff=0.12))
    for y in ys:
        g.add(layout.label(y_fmt.format(y), font_size=font_size, color=P.FG)
              .next_to(ax.c2p(x0, y), LEFT, buff=0.12))
    return g


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
        rb = dataio.section("preface")["finale"]["rubin"]
        rubin = Line(m.p(180.0 + 179.9, rb["dec_max"]), m.p(180.0 - 179.9, rb["dec_max"]),
                     color=P.PURPLE, stroke_width=2.2)
        rubin_lab = layout.label(f"Rubin's limit, {rb['dec_max']:+.0f}°", font_size=14,
                                 color=P.PURPLE)
        rubin_lab.next_to(m.p(345, rb["dec_max"]), UP, buff=0.08).align_to(m.frame, LEFT)
        rubin_lab.shift(RIGHT * 0.15)
        cap2 = layout.caption(
            f"Rubin will search south of {rb['dec_max']:+.0f}° anyway: "
            "spend the space telescope where it can't.", font_size=20)
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
            card.add_background_rectangle(color=P.BG, opacity=0.85, buff=0.06)
            card.next_to(tiles, CARD_SIDE[zone["zone"]], buff=0.18)
            if card.get_right()[0] > 6.6:
                card.shift(LEFT * (card.get_right()[0] - 6.6))
            if card.get_left()[0] < -6.6:
                card.shift(RIGHT * (-6.6 - card.get_left()[0]))
            line = VGroup(
                layout.label(f"{zone['why'][0].upper()}{zone['why'][1:]}.", font_size=19,
                             color=P.FG),
                layout.label(
                    f"median planet V {zone['median_v']:.1f} at "
                    f"{zone['median_dist_au']:.0f} AU  ·  "
                    f"{zone['median_integration_s']:.0f} s visits  ·  {zone['hours']:.0f} h  ·  "
                    f"best in {_season(zone['opposition_months'])}",
                    font_size=16, color=colour),
            ).arrange(DOWN, buff=0.1).to_edge(DOWN, buff=0.3)
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
        ax = widgets.axes([0, 8000, 1000], [0, ymax, 5], x_length=8.0, y_length=4.3,
                          font_size=14, shift_down=0.0)
        ax.move_to([-1.7, 0.2, 0])
        nums = _tick_numbers(ax, range(2000, 8001, 2000), range(5, ymax + 1, 5))
        labels = VGroup(
            nums,
            layout.label("space-telescope hours", font_size=16, color=P.FG)
            .next_to(ax, DOWN, buff=0.45),
            layout.label("chance of finding Planet Nine (%)", font_size=16, color=P.FG)
            .rotate(np.pi / 2).next_to(ax, LEFT, buff=0.45),
        )
        self.play(Create(ax), FadeIn(labels))
        colours = [P.PURPLE, P.ORANGE, P.TEAL]
        legend = VGroup()
        for pol, col in zip(d["policies"][:3], colours):
            xs = [0.0] + list(budgets)
            ys = [0.0] + [100 * c for c in pol["captured"]]
            curve = widgets.curve(ax, xs, ys, color=col, stroke_width=3.2)
            row = VGroup(
                layout.label(pol["policy"], font_size=15, color=col, weight="BOLD"),
                layout.label(f"{100 * pol['unique']:.0f}% of the prediction in play",
                             font_size=13, color=P.FG),
            ).arrange(DOWN, buff=0.05, aligned_edge=LEFT)
            legend.add(row)
            legend.arrange(DOWN, buff=0.28, aligned_edge=LEFT)
            legend.move_to([4.75, 1.3, 0])
            self.play(Create(curve), FadeIn(row), run_time=1.1)
        ref = Dot(ax.c2p(d["reference_hours"], 100 * d["reference_captured"]),
                  radius=0.09, color=P.FG)
        ref_lab = layout.label(
            f"reference campaign: {d['reference_hours']:.0f} h, "
            f"{100 * d['reference_captured']:.1f}%", font_size=14, color=P.FG)
        ref_lab.next_to(ref, DOWN + RIGHT, buff=0.1)
        cap0 = layout.caption(
            "Racing Rubin everywhere buys the most, but mostly by duplicating Rubin.",
            font_size=20)
        self.play(FadeIn(cap0))
        timing.hold_to_read(self, cap0, legend, settle=0.8)
        cap = layout.caption(
            f"The reference plan concedes Rubin's sky and still covers "
            f"{100 * d['reference_share_of_unique']:.0f}% of what only it can reach.",
            font_size=20)
        self.play(FadeIn(ref), FadeIn(ref_lab), FadeOut(cap0), FadeIn(cap))
        timing.hold_to_read(self, cap, settle=1.2)
        self.play(*[FadeOut(x) for x in self.mobjects[1:]])

        # 2. the calendar
        hours = d["hours_by_month"]
        order = [(k + 6) % 12 for k in range(12)]
        by_zone = {z: [0.0] * 12 for z in ZONE_COLOURS}
        for t in d["tiles"]:
            by_zone[t["zone"]][t["opposition_month"] - 1] += t["hours"]
        ymax = int(np.ceil(max(hours) / 100.0) * 100)
        ax2 = widgets.axes([0, 12, 1], [0, ymax, 100], x_length=9.6, y_length=3.9,
                           shift_down=-0.1)
        lab2 = VGroup(
            _tick_numbers(ax2, (), range(100, ymax + 1, 100)),
            layout.label("telescope hours per month", font_size=16, color=P.FG)
            .rotate(np.pi / 2).next_to(ax2, LEFT, buff=0.6),
        )
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
            "Each region is imaged at opposition, when the planet drifts fastest.",
            font_size=20)
        self.play(FadeIn(ticks), FadeIn(bars, lag_ratio=0.05), FadeIn(key), FadeIn(cap2),
                  run_time=1.6)
        timing.hold_to_read(self, cap2, settle=1.0)
        cap3 = layout.caption(
            f"{d['telescope']['n_epochs']} visits per field, hours apart: "
            "a star stays put, Planet Nine creeps.", font_size=20)
        self.play(FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=1.0)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, f"{d['reference_hours']:.0f} hours in one year, all on sky "
                  "no other survey will search.")
