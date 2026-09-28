"""Shankman et al. (2017) -- OSSOS VI: striking biases in the detection of large
semimajor axis trans-Neptunian objects.

OSSOS recorded where it pointed and how deep it went, so its selection effects
are known. Put a population with perihelia in every direction through a survey
that looks in four patches of sky and the objects it finds look clustered. The
synthetic population, the survey blocks, the detections and the statistics are
the crate's own (anim.json -> papers -> p9-2017-ossos-bias). The four real OSSOS
objects beyond 230 AU come from the p9-2019-clustering sample.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Create,
    FadeIn,
    FadeOut,
    Line,
    Rectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2017-ossos-bias"


class OssosBias2017(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        lost = d["lost"]
        hist = d["varpi"]

        self.add(paper.scene_header(CRATE))

        # 1. a population with no preferred direction
        m = sky.SkyMap(width=12.0, dec_range=(-60, 60), centre=(0.0, 0.45, 0.0))
        ecl, gal = m.reference_curves()
        everyone = m.dots(d["sky_bright"], color=P.TEAL, radius=0.03, opacity=0.7)
        cap = layout.caption(
            "A made-up population: perihelia spread evenly around the sky", font_size=22)
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), run_time=1.0)
        self.play(FadeIn(everyone, lag_ratio=0.003), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.8)

        # 2. the survey looks in four patches
        blocks = VGroup(*[
            m.polyline(b["outline"], color=P.PURPLE, stroke_width=2.2) for b in d["blocks"]])
        found = m.dots(d["sky_detected"], color=P.GREEN, radius=0.04, opacity=0.95)
        key = m.legend([("ecliptic", P.ORANGE), ("galactic plane ±10°", P.PURPLE),
                        ("survey blocks", P.PURPLE), ("found", P.GREEN)])
        key.next_to(m.frame, DOWN, buff=0.75)
        cap2 = layout.caption(
            f"The survey looks in {len(d['blocks'])} blocks along the ecliptic, to magnitude "
            f"{d['limiting_mag']:.1f}", font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(Create(blocks), FadeIn(key), FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, settle=0.6)
        cap3 = layout.caption(
            f"It finds {100 * lost['detected']:.0f}%: the ones at perihelion inside a block; "
            f"{100 * lost['outside_blocks']:.0f}% were simply elsewhere", font_size=22)
        self.play(FadeOut(cap2), run_time=0.4)
        self.play(everyone.animate.set_opacity(0.2), FadeIn(found, lag_ratio=0.01),
                  FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, everyone, blocks, found, key, cap3)))

        # 3. what that does to the perihelion directions
        edges = np.array(hist["edges_deg"])
        parent = 100 * np.array(hist["parent"])
        seen = 100 * np.array(hist["detected"])
        top = float(np.ceil(seen.max() / 2.0) * 2 + 2)
        ax, labels = widgets.labeled_axes(
            [0, 360, 45], [0, top, 2], x_label="longitude of perihelion (degrees)",
            y_label="share of objects (%)", y_rotate=True, numbers=True,
            x_length=10.5, y_length=4.2, shift_down=-0.3)
        bars_parent = widgets.histogram(ax, edges, parent, color=P.TEAL, opacity=0.35)
        bars_seen = widgets.histogram(ax, edges, seen, color=P.GREEN, opacity=0.6)
        legend = VGroup(
            layout.label(f"the whole population: alignment {d['r_bar_parent']:.2f}",
                         font_size=16, color=P.TEAL),
            layout.label(f"the {d['n_detected']:,} found: alignment {d['r_bar_detected']:.2f}",
                         font_size=16, color=P.GREEN),
            layout.label("(0 = every direction, 1 = one direction)", font_size=13,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT)
        legend.move_to(ax.c2p(130, 0.82 * top))
        cap4 = layout.caption(
            "A population with no preferred direction comes back looking clustered",
            font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(bars_parent, lag_ratio=0.03), run_time=1.0)
        self.play(FadeIn(bars_seen, lag_ratio=0.05), FadeIn(legend), FadeIn(cap4), run_time=1.4)
        timing.hold_to_read(self, cap4, legend, settle=1.2)

        # 4. the real survey's objects against that expectation
        real = VGroup(*[
            Line(ax.c2p(o["varpi_deg"], 0), ax.c2p(o["varpi_deg"], 0.55 * top), color=P.FG,
                 stroke_width=3.5) for o in d["ossos"]])
        names = VGroup()
        for k, o in enumerate(sorted(d["ossos"], key=lambda o: o["varpi_deg"])):
            lab = layout.label(o["name"], font_size=14, color=P.FG)
            lab.next_to(ax.c2p(o["varpi_deg"], 0.55 * top), UP, buff=0.08)
            lab.shift(UP * 0.3 * (k % 2))
            plate = Rectangle(width=lab.width + 0.1, height=lab.height + 0.08, stroke_width=0)
            plate.set_fill(P.BG, opacity=0.9).move_to(lab)
            names.add(VGroup(plate, lab))
        stat = VGroup(
            layout.label(f"alignment of the {d['n_ossos']} real objects: "
                         f"{d['r_bar_ossos']:.2f}", font_size=16, color=P.FG),
            layout.label(f"a uniform sky seen through OSSOS does at least this well "
                         f"{100 * d['p_chance']:.0f}% of the time", font_size=14,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        stat.move_to(ax.c2p(95, 0.9 * top), aligned_edge=LEFT)
        cap5 = layout.caption(
            f"The {d['n_ossos']} real OSSOS orbits beyond 230 AU look no more aligned "
            f"than chance", font_size=22)
        self.play(FadeOut(legend), FadeOut(cap4), run_time=0.5)
        self.play(FadeIn(real, lag_ratio=0.2), FadeIn(names), FadeIn(cap5), run_time=1.4)
        self.play(FadeIn(stat))
        timing.hold_to_read(self, cap5, settle=1.4)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Where a survey looks shapes what it finds; OSSOS alone shows no clustering.")
