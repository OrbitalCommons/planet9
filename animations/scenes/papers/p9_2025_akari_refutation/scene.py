"""Chen et al. (2025) -- a far-infrared search for Planet Nine using the AKARI
all-sky survey.

AKARI scanned each patch of sky again after six months. Over that interval a
planet at 300-800 AU shifts by 9-23 arcminutes of parallax, far more than its
orbital drift, so it shows up as a source that is steady for a day and absent
half a year later. The crate computes that motion; the search region, the
selection counts and the two candidates are the paper's, exported as labelled
published values. Everything is from anim.json -> papers ->
p9-2025-akari-refutation.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Rectangle,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2025-akari-refutation"


def rows(items, font_size=15, buff=0.16):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


class Zoom(VGroup):
    """A plate-carree close-up of an RA/Dec box, east to the left."""

    def __init__(self, region, height, centre):
        super().__init__()
        self.region = region
        self.h = height
        self.w = height * (region["ra_hi"] - region["ra_lo"]) / (region["dec_hi"] - region["dec_lo"])
        self.c = np.array(centre, dtype=float)
        frame = Rectangle(width=self.w, height=self.h, color=P.PURPLE, stroke_width=1.6)
        frame.set_fill("#16171f", opacity=1.0).move_to(self.c)
        self.add(frame)
        for ra in (region["ra_lo"], region["ra_hi"]):
            lab = layout.label(f"{ra / 15:.1f}h", font_size=12, color=P.MUTED)
            self.add(lab.move_to(self.p(ra, region["dec_lo"]) + DOWN * 0.2))
        for dec in (region["dec_lo"], 0.0, region["dec_hi"]):
            lab = layout.label(f"{dec:+.0f}°", font_size=12, color=P.MUTED)
            self.add(lab.move_to(self.p(region["ra_hi"], dec) + LEFT * 0.35))

    def p(self, ra, dec):
        r = self.region
        u = (ra - 0.5 * (r["ra_lo"] + r["ra_hi"])) / (r["ra_hi"] - r["ra_lo"])
        v = (dec - 0.5 * (r["dec_lo"] + r["dec_hi"])) / (r["dec_hi"] - r["dec_lo"])
        return self.c + np.array([-u * self.w, v * self.h, 0.0])

    def inside(self, ra, dec):
        r = self.region
        return r["ra_lo"] <= ra <= r["ra_hi"] and r["dec_lo"] <= dec <= r["dec_hi"]

    def curve(self, radec, color):
        pts = [self.p(ra, dec) for ra, dec in radec if self.inside(ra, dec)]
        pts.sort(key=lambda q: q[0])
        m = VMobject(color=color, stroke_width=1.6).set_stroke(opacity=0.8)
        if len(pts) > 1:
            m.set_points_as_corners(pts)
        return m


class AkariRefutation2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        mo = d["motion"]
        near, far = d["search_distance_au"]

        self.add(paper.scene_header(CRATE))

        # 1. six months of motion: parallax dwarfs the orbital drift
        dist = np.array(mo["distance_au"])
        ax, labels = widgets.labeled_axes(
            [200, 800, 100], [0, 30, 5], x_label="heliocentric distance (AU)",
            y_label="motion in six months (arcmin)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.3)
        par = widgets.curve(ax, dist, mo["parallax_arcmin"], color=P.ORANGE)
        drift = Polygon(
            *[ax.c2p(x, y) for x, y in zip(dist, mo["pm_perihelion_arcmin"])],
            *[ax.c2p(x, y) for x, y in zip(dist[::-1], mo["pm_circular_arcmin"][::-1])],
            stroke_width=1.0, color=P.BLUE).set_fill(P.BLUE, opacity=0.45)
        lo, hi = d["parallax_range_arcmin"]
        pm_lo, pm_hi = d["proper_motion_range_arcmin"]
        key = rows([
            ("parallax: Earth crosses its orbit", P.ORANGE),
            (f"{lo:.0f}′ to {hi:.0f}′ at {far:.0f} to {near:.0f} AU", P.ORANGE),
            (f"paper: {pub['parallax_arcmin'][0]:.0f}′ to {pub['parallax_arcmin'][1]:.0f}′",
             P.MUTED),
            ("the planet's own orbital drift", P.BLUE),
            (f"{pm_lo:.1f}′ to {pm_hi:.1f}′", P.BLUE),
            (f"paper: {pub['proper_motion_arcmin'][0]:.1f}′ to "
             f"{pub['proper_motion_arcmin'][1]:.1f}′", P.MUTED),
        ])
        key[3:].shift(DOWN * 0.25)
        key.move_to([4.6, 1.2, 0])
        cap = layout.caption(
            "AKARI returned to each field after six months; by then a distant planet "
            "has jumped aside", font_size=21)
        self.play(Create(ax), FadeIn(labels), run_time=1.0)
        self.play(Create(par), FadeIn(key[:3]), FadeIn(cap), run_time=1.2)
        self.play(FadeIn(drift), FadeIn(key[3:]), run_time=1.0)
        timing.hold_to_read(self, cap, key, settle=0.8)
        cap_b = layout.caption(
            "So look for a source that is steady for a day and gone half a year later",
            font_size=21)
        self.play(FadeOut(cap), FadeIn(cap_b))
        timing.hold_to_read(self, cap_b, settle=0.6)
        self.play(FadeOut(VGroup(ax, labels, par, drift, key, cap_b)))

        # 2. where they looked, what was left
        region = d["region"]
        m = sky.SkyMap(width=5.4, dec_range=(-90, 90), centre=(-3.75, 0.55, 0.0))
        ecl, gal = m.reference_curves()
        box = m.box(region["ra_lo"], region["ra_hi"], region["dec_lo"], region["dec_hi"],
                    opacity=0.3)
        old = d["earlier_pair"]
        old_mark = Circle(radius=0.09, color=P.MUTED, stroke_width=2).move_to(
            m.p(old["iras_ra_deg"], old["iras_dec_deg"]))
        old_lab = layout.label("earlier IRAS-AKARI pair", font_size=12, color=P.MUTED)
        old_lab.next_to(old_mark, DOWN, buff=0.06)
        z = Zoom(region, height=4.6, centre=(0.75, 0.2, 0.0))
        s = dataio.section("sky")
        z_ecl = VGroup(z.curve(s["ecliptic"], P.ORANGE),
                       layout.label("ecliptic", font_size=12, color=P.ORANGE)
                       .move_to(z.c + np.array([0.0, 1.35, 0.0])))
        ties = VGroup(*[
            Line(m.p(region["ra_lo"], dec), z.p(region["ra_hi"], dec), color=P.PURPLE,
                 stroke_width=1.0).set_stroke(opacity=0.5)
            for dec in (region["dec_lo"], region["dec_hi"])])
        steps = d["funnel"]
        table = VGroup()
        for k, st in enumerate(steps):
            last = k == len(steps) - 1
            col = P.GREEN if last else P.FG
            n = layout.label(f"{st['n']:,}", font_size=18 if last else 16, color=col,
                             weight="BOLD")
            t = layout.label(st["label"], font_size=14, color=col)
            t.set_opacity(1.0 if last else 0.7)
            table.add(VGroup(n, t).arrange(DOWN, buff=0.04, aligned_edge=LEFT))
        table.arrange(DOWN, buff=0.17, aligned_edge=LEFT).move_to([4.75, 0.15, 0])
        area = ((region["ra_hi"] - region["ra_lo"]), (region["dec_hi"] - region["dec_lo"]))
        cap2 = layout.caption(
            f"The search region, {area[0]:.0f}° by {area[1]:.0f}°, where simulations "
            "placed the planet (counts: paper)", font_size=21)
        self.play(FadeIn(m), Create(ecl), Create(gal), run_time=0.9)
        self.play(FadeIn(box), FadeIn(z), Create(z_ecl), Create(ties), FadeIn(cap2),
                  run_time=1.0)
        self.play(FadeIn(table[:-1], lag_ratio=0.25), run_time=2.0)
        timing.hold_to_read(self, cap2, table[:-1], settle=0.5)

        cands = d["candidates"]
        dots = VGroup(*[Dot(z.p(c["ra_deg"], c["dec_deg"]), radius=0.07, color=P.GREEN)
                        for c in cands])
        tags = VGroup()
        for k in range(len(cands)):
            tag = layout.label(f"seen {cands[k]['epoch']}", font_size=13, color=P.GREEN)
            tags.add(tag.next_to(dots[k], RIGHT, buff=0.12))
        sep = d["pair_separation_arcmin"]
        cap3 = layout.caption(
            f"Two candidates, each seen once. They lie {sep / 60:.1f}° apart; parallax "
            f"moves a planet {hi / 60:.1f}° at most", font_size=21)
        self.play(FadeIn(dots), FadeIn(tags), FadeIn(table[-1]), FadeIn(old_mark),
                  FadeIn(old_lab), FadeOut(cap2), FadeIn(cap3), run_time=1.2)
        timing.hold_to_read(self, cap3, settle=1.2)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, "Two new far-infrared candidates, neither confirmed by a second sighting.")
