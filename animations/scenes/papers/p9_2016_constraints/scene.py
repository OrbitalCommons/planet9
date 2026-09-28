"""Brown & Batygin (2016) -- observational constraints on the orbit and location
of Planet Nine.

Two months after the hypothesis, the follow-up asks two practical questions:
which planets actually reproduce the alignment, and where on the sky could such
a planet still be hiding. The grid of trial planets, the nominal orbit's path
across the sky and its brightness along that path are the crate's own
(anim.json -> papers -> p9-2016-constraints). The allowed region of the grid
comes from the paper's 4-Gyr simulations, which the crate does not rerun at
film scale; it is drawn as the published range and labelled as such.
"""
import numpy as np
from manim import (
    DOWN,
    RIGHT,
    UP,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Polygon,
    Rectangle,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2016-constraints"

# Published in the paper's abstract; drawn only as a labelled comparison.
PAPER_A = (380.0, 980.0)
PAPER_Q = (150.0, 350.0)
PAPER_MASS = (5.0, 20.0)


def in_paper_range(g):
    return PAPER_A[0] <= g["a"] <= PAPER_A[1] and PAPER_Q[0] <= g["q"] <= PAPER_Q[1]


def plated(text, color, font_size=16):
    """A label on a background plate, readable where map curves cross it."""
    lab = layout.label(text, font_size=font_size, color=color)
    box = Rectangle(width=lab.width + 0.12, height=lab.height + 0.12, stroke_width=0)
    box.set_fill("#16171f", opacity=0.9).move_to(lab)
    return VGroup(box, lab)


class Constraints2016(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        grid = d["grid"]
        orbit = d["orbit"]
        depth = d["wide_depth"]

        self.add(paper.scene_header(CRATE))

        # 1. the grid of trial planets, and the part of it that works
        ax, labels = widgets.labeled_axes(
            [200, 2000, 200], [0, 1.1, 0.2], x_label="semi-major axis of the trial planet (AU)",
            y_label="eccentricity", y_rotate=True, numbers=True,
            x_length=10.5, y_length=4.5, shift_down=-0.3)
        dots = VGroup(*[Dot(ax.c2p(g["a"], g["e"]), radius=0.04, color=P.MUTED) for g in grid])
        cap = layout.caption(
            f"{d['grid_total']} trial planets: every dot is one orbit, tried at five masses",
            font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(dots, lag_ratio=0.01), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.8)

        a_run = np.linspace(PAPER_A[0], PAPER_A[1], 40)
        upper = [ax.c2p(a, 1.0 - PAPER_Q[0] / a) for a in a_run]
        lower = [ax.c2p(a, max(1.0 - PAPER_Q[1] / a, 0.0)) for a in a_run[::-1]]
        region = Polygon(*upper, *lower, color=P.TEAL, stroke_width=2)
        region.set_fill(P.TEAL, opacity=0.12)
        good = VGroup(*[Dot(ax.c2p(g["a"], g["e"]), radius=0.05, color=P.TEAL)
                        for g in grid if in_paper_range(g)])
        nominal = Dot(ax.c2p(d["a"], d["e"]), radius=0.1, color=P.BLUE).set_z_index(4)
        key = VGroup(
            layout.label("paper: planets that reproduce the alignment", font_size=15,
                         color=P.TEAL),
            layout.label(f"nominal orbit:  a = {d['a']:.0f} AU,  e = {d['e']:.1f},  "
                         f"{d['mass']:.0f} Earth masses", font_size=15, color=P.BLUE),
        ).arrange(RIGHT, buff=0.7)
        key.move_to(ax.c2p(1100, 1.03))
        cap2 = layout.caption(
            f"Paper: only a = {PAPER_A[0]:.0f}-{PAPER_A[1]:.0f} AU with perihelion "
            f"{PAPER_Q[0]:.0f}-{PAPER_Q[1]:.0f} AU and "
            f"{PAPER_MASS[0]:.0f}-{PAPER_MASS[1]:.0f} Earth masses works", font_size=22)
        # The edges of the allowed band are curves of constant perihelion q = a(1 - e).
        iso = VGroup()
        for q in PAPER_Q:
            a_iso = np.linspace(max(q / 0.95, 200.0), 2000.0, 60)
            iso.add(widgets.curve(ax, a_iso, 1.0 - q / a_iso, color=P.MUTED, stroke_width=1.6))
            tag = layout.label(f"perihelion {q:.0f} AU", font_size=15, color=P.FG)
            tag.next_to(ax.c2p(1500, 1.0 - q / 1500), DOWN, buff=0.08)
            iso.add(tag)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(FadeIn(region), FadeIn(good, lag_ratio=0.05), FadeIn(nominal), FadeIn(key),
                  Create(iso), FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, key, settle=1.0)
        self.play(FadeOut(VGroup(ax, labels, dots, region, good, nominal, key, iso, cap2)))


        # 2. where that planet would be on the sky, and how bright
        m = sky.SkyMap(width=10.0, dec_range=(-60, 60), centre=(-1.45, 0.4, 0.0))
        ecl, gal = m.reference_curves()
        bright = [s for s in orbit if s["v_mag"] <= depth]
        faint = [s for s in orbit if s["v_mag"] > depth]
        path_b = m.dots(bright, color=P.RED, radius=0.035, opacity=0.95)
        path_f = m.dots(faint, color=P.TEAL, radius=0.035, opacity=0.95)
        peri = min(orbit, key=lambda s: s["r_au"])
        apo = max(orbit, key=lambda s: s["r_au"])
        peri_pt = m.p(peri["ra_deg"], peri["dec_deg"])
        apo_pt = m.p(apo["ra_deg"], apo["dec_deg"])
        peri_lab = plated(f"perihelion {peri['r_au']:.0f} AU,  V = {peri['v_mag']:.1f}",
                          P.RED, font_size=14)
        peri_lab.next_to(peri_pt, DOWN + RIGHT, buff=0.1)
        apo_lab = plated(f"aphelion {apo['r_au']:,.0f} AU,  V = {apo['v_mag']:.1f}",
                         P.TEAL, font_size=14)
        apo_lab.next_to(apo_pt, UP + RIGHT, buff=0.1)
        key2 = m.legend([("ecliptic", P.ORANGE), ("galactic plane", P.PURPLE),
                         (f"brighter than V = {depth:.1f}", P.RED),
                         ("fainter", P.TEAL)])
        key2.next_to(m.frame, DOWN, buff=0.72)
        cap3 = layout.caption(
            f"The nominal orbit across the sky, coloured by the depth of {d['wide_survey']}",
            font_size=22)
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), run_time=1.0)
        self.play(FadeIn(path_b, lag_ratio=0.02), FadeIn(path_f, lag_ratio=0.02),
                  FadeIn(key2), FadeIn(cap3), run_time=1.8)
        self.play(FadeIn(peri_lab), FadeIn(apo_lab))
        timing.hold_to_read(self, cap3, key2, settle=0.8)

        col_x = 5.35
        tally = paper.result_readout("path already surveyable",
                                     f"{100 * d['orbit_fraction_wide']:.0f}%",
                                     color=P.RED).scale(0.78)
        tally.move_to([col_x, 1.55, 0])
        paper_note = layout.label("paper: about two-thirds", font_size=15, color=P.MUTED)
        paper_note.next_to(tally, DOWN, buff=0.1)
        cap4 = layout.caption(
            "Wide surveys could already see it along the red two-thirds of the path",
            font_size=22)
        self.play(FadeOut(cap3), run_time=0.4)
        self.play(FadeIn(tally), FadeIn(paper_note), FadeIn(cap4))
        timing.hold_to_read(self, cap4, settle=0.8)

        # 3. but the planet does not spend its time evenly along the path
        snaps = d["snapshots"]
        step = d["period_yr"] / len(snaps)
        snap_dots = VGroup(*[
            Dot(m.p(s["ra_deg"], s["dec_deg"]), radius=0.06,
                color=P.RED if s["v_mag"] <= depth else P.TEAL).set_z_index(3)
            for s in snaps])
        n_bright = sum(1 for s in snaps if s["v_mag"] <= depth)
        cap5 = layout.caption(
            f"Now one dot every {step:,.0f} years: it races past perihelion "
            f"and crawls near aphelion", font_size=22)
        self.play(FadeOut(cap4), run_time=0.4)
        self.play(path_b.animate.set_opacity(0.25), path_f.animate.set_opacity(0.25),
                  FadeOut(peri_lab), FadeOut(apo_lab), FadeIn(cap5))
        self.play(LaggedStart(*[FadeIn(x, scale=1.8) for x in snap_dots], lag_ratio=0.12),
                  run_time=5.0)
        timing.hold_to_read(self, cap5, settle=0.6)

        tally2 = paper.result_readout("time spent fainter",
                                      f"{100 * d['time_share_too_faint']:.0f}%",
                                      color=P.TEAL).scale(0.78)
        tally2.move_to([col_x, -0.45, 0])
        count = layout.label(f"{len(snaps) - n_bright} of {len(snaps)} dots", font_size=15,
                             color=P.MUTED)
        count.next_to(tally2, DOWN, buff=0.1)
        faint_v = [s["v_mag"] for s in snaps if s["v_mag"] > depth]
        cap6 = layout.caption(
            f"So it most likely sits near aphelion today, at V ≈ {min(faint_v):.1f}-"
            f"{max(faint_v):.1f}  (paper: 22-25)", font_size=22)
        self.play(FadeOut(cap5), run_time=0.4)
        self.play(FadeIn(tally2), FadeIn(count), FadeIn(cap6))
        timing.hold_to_read(self, cap6, tally2, settle=1.2)
        self.play(FadeOut(cap6), FadeOut(key2))

        layout.show_takeaway(
            self, "A narrow family of orbits works, and the planet is most likely faint, near aphelion.")
