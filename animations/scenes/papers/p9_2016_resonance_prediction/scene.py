"""Malhotra, Volk & Wang (2016) -- corralling a distant planet with resonances.

The four longest-period objects known in 2016 have periods close to simple
ratios of one another. If each sits in an N:1 or N:2 mean-motion resonance with
an unseen outer planet, the planet's period follows from theirs. The scene
slides a candidate planet outward, carrying its ladder of resonant distances
across the four real orbits, then shows the crate's scan of the misfit against
the planet's semi-major axis. Reproduced in p9-2016-resonance-prediction; every
number is the crate's own (anim.json -> papers -> p9-2016-resonance-prediction).
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
    NumberLine,
    Scene,
    ValueTracker,
    VGroup,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2016-resonance-prediction"

A_LO, A_HI = 140.0, 720.0
LINE_Y = -0.4


class ResonancePrediction2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        best = d["best_a9_au"]
        matches = d["matches_best"]
        ladder = d["ladder"]

        self.add(paper.scene_header(CRATE))

        # 1. the four longest-period objects on a distance line
        nl = NumberLine(x_range=[A_LO, A_HI, 100], length=12.4, color=P.MUTED,
                        include_tip=False, stroke_width=2)
        nl.move_to([0, LINE_Y, 0])
        nums = VGroup(*[
            layout.label(f"{a}", font_size=15, color=P.MUTED).next_to(nl.n2p(a), DOWN, buff=1.05)
            for a in range(200, 701, 100)])
        axis_name = layout.label("semi-major axis (AU)", font_size=16, color=P.MUTED)
        axis_name.next_to(nums, DOWN, buff=0.15)
        self.play(Create(nl), FadeIn(nums), FadeIn(axis_name))

        etnos = VGroup()
        order = sorted(matches, key=lambda m: m["a_au"])
        for k, m in enumerate(order):
            p = nl.n2p(m["a_au"])
            dot = Dot(p, radius=0.09, color=P.GREEN).set_z_index(4)
            top = p + UP * (1.2 + 0.95 * (k % 2))
            lead = Line(p, top, color=P.GREEN, stroke_width=1.2).set_stroke(opacity=0.6)
            name = layout.label(m["name"], font_size=19, color=P.GREEN, weight="BOLD")
            per = layout.label(f"{m['period_yr']:,.0f} yr", font_size=16, color=P.FG)
            tag = VGroup(name, per).arrange(DOWN, buff=0.06).next_to(top, UP, buff=0.06)
            etnos.add(VGroup(dot, lead, tag))
        cap = layout.caption("The four longest-period objects known in 2016, by orbital period",
                             font_size=22)
        self.play(FadeIn(etnos, lag_ratio=0.2), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.8)

        # 2. a candidate planet and the resonant distances it implies
        a9 = ValueTracker(A_HI - 20.0)

        def x_of(a):
            return nl.n2p(a)[0]

        planet = Dot(radius=0.13, color=P.BLUE).set_z_index(5)
        planet_lab = layout.label("candidate planet", font_size=17, color=P.BLUE, weight="BOLD")

        def place_planet(_):
            planet.move_to([x_of(a9.get_value()), LINE_Y, 0])
            planet_lab.next_to(planet, UP, buff=0.22)

        planet.add_updater(place_planet)
        place_planet(None)

        rungs = VGroup()
        ladder = sorted(ladder, key=lambda r: r["a_res_au"])
        for k, r in enumerate(ladder):
            frac = (r["q"] / r["p"]) ** (2.0 / 3.0)
            tick = Line([0, LINE_Y - 0.28, 0], [0, LINE_Y + 0.28, 0], color=P.ORANGE,
                        stroke_width=2.5)
            lab = layout.label(f"{r['p']}:{r['q']}", font_size=14, color=P.ORANGE)
            lab.next_to(tick, DOWN, buff=0.06 + 0.28 * (k % 2))
            rung = VGroup(tick, lab)

            def place(mob, frac=frac):
                a = frac * a9.get_value()
                mob.set_x(x_of(a))
                mob.set_opacity(1.0 if a >= A_LO + 5 else 0.0)

            rung.add_updater(place)
            place(rung)
            rungs.add(rung)

        cap2 = layout.caption(
            "An outer planet holds objects at fixed period ratios: N:1 and N:2",
            font_size=22)
        self.play(FadeIn(planet), FadeIn(planet_lab), FadeIn(rungs), FadeOut(cap), FadeIn(cap2),
                  FadeOut(nums[:]), run_time=1.0)
        nums2 = VGroup(*[
            layout.label(f"{a}", font_size=15, color=P.MUTED).next_to(nl.n2p(a), DOWN, buff=1.05)
            for a in range(200, 701, 100)])
        self.add(nums2)
        self.play(a9.animate.set_value(560.0), run_time=3.0, rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap2, settle=0.4)

        cap3 = layout.caption(
            f"At {best:.0f} AU the ladder lands on all four at once", font_size=22)
        self.play(a9.animate.set_value(best), FadeOut(cap2), FadeIn(cap3), run_time=2.4)
        planet.clear_updaters()
        hit = {(m["p"], m["q"]) for m in matches}
        dim = []
        for rung, r in zip(rungs, ladder):
            rung.clear_updaters()
            if (r["p"], r["q"]) not in hit:
                dim.append(rung.animate.set_opacity(0.3))
        self.play(*dim, run_time=0.8)
        timing.hold_to_read(self, cap3, settle=0.6)
        sedna = next(m for m in matches if m["name"] == "Sedna")
        vp = next(m for m in matches if m["name"] == "2012 VP113")
        cap3b = layout.caption(
            f"{sedna['p']}:{sedna['q']}: Sedna laps {sedna['p']} times while the planet laps "
            f"{sedna['q']};  {vp['p']}:{vp['q']}: VP113 laps {vp['p']} times per planet orbit",
            font_size=21)
        self.play(FadeOut(cap3), FadeIn(cap3b))
        timing.hold_to_read(self, cap3b, settle=1.0)
        line_group = VGroup(nl, nums2, axis_name, etnos, planet, planet_lab, rungs, cap3b)
        self.play(FadeOut(line_group))

        # 3. the scan: how sharply do the four orbits pick the planet?
        s = d["scan"]
        ax, labels = widgets.labeled_axes(
            [500, 900, 100], [0, 0.3, 0.1],
            x_label="candidate planet semi-major axis (AU)",
            y_label="average miss from an N:1 or N:2 ratio", y_rotate=True, numbers=True,
            x_length=10.0, y_length=4.3, shift_down=-0.35)
        curve = widgets.curve(ax, s["a9_au"], s["residual"], color=P.ORANGE, stroke_width=2.5)
        self.play(Create(ax), FadeIn(labels))
        self.play(Create(curve), run_time=2.0)
        here = widgets.marker_line(ax, best, (0, 0.27), f"best fit here  {best:.0f} AU",
                                   color=P.BLUE, side=LEFT)
        pub = widgets.marker_line(ax, d["published_a9_au"], (0, 0.22),
                                  f"paper  {d['published_a9_au']:.0f} AU", color=P.FG,
                                  side=RIGHT)
        self.play(Create(here), Create(pub))
        tally = paper.result_readout(
            "implied orbital period", f"{d['best_period_yr']:,.0f} yr", color=P.BLUE).scale(0.8)
        tally.move_to(ax.c2p(820, 0.055))
        cap4 = layout.caption(
            f"Paper: {d['published_period_yr']:,.0f} yr.  Only "
            f"{100 * d['fraction_as_good_as_published']:.0f}% of candidate distances fit "
            f"as well as {d['published_a9_au']:.0f} AU",
            font_size=22)
        self.play(FadeIn(tally), FadeIn(cap4))
        timing.hold_to_read(self, cap4, tally, settle=1.2)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, f"Four orbital periods point to one planet near {best:.0f} AU.")
