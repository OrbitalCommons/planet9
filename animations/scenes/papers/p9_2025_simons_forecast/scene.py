"""Simons Observatory Collaboration (2025) -- science goals and forecasts for the
enhanced Large Aperture Telescope.

The Planet Nine entry of the forecast: a 5 Earth-mass planet radiating at ~40 K
is a millimetre point source whose flux falls only as 1/d^2. The crate computes
that flux (shared p9-core thermal model), turns each survey depth into a
signal-to-noise curve and a reach in distance, and sets the reach against the
Brown & Batygin (2021) reference population on the mass-distance plane.
Everything drawn is from anim.json -> papers -> p9-2025-simons-forecast.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Polygon,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-simons-forecast"

SURVEYS = [("act", P.PURPLE), ("so_shallow", P.ORANGE), ("so_deep", P.TEAL)]


def clipped(ax, xs, ys, y_max, color, stroke_width=3.0):
    """A polyline through the points that lie inside the axes."""
    pts = [ax.c2p(x, y) for x, y in zip(xs, ys) if y <= y_max]
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners(pts)
    return m


def line_key(rows, font_size=15, buff=0.16):
    """Stacked (colour, text) keys with a line swatch."""
    g = VGroup()
    for col, text in rows:
        sw = Line(LEFT * 0.2, RIGHT * 0.2, color=col, stroke_width=4)
        g.add(VGroup(sw, layout.label(text, font_size=font_size, color=col))
              .arrange(RIGHT, buff=0.12))
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


class SimonsForecast2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        body = d["body"]

        self.add(paper.scene_header(CRATE))

        # 1. signal-to-noise of a 5 Earth-mass planet against distance
        dist = np.array(d["distance_au"])
        ax, labels = widgets.labeled_axes(
            [300, 1300, 100], [0, 10, 2], x_label="heliocentric distance (AU)",
            y_label="signal-to-noise", y_rotate=True, numbers=True,
            x_length=7.6, y_length=4.2, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.2)
        limit = DashedLine(ax.c2p(300, 2), ax.c2p(1300, 2), color=P.MUTED, stroke_width=2)
        limit_lab = layout.label("95% detection limit", font_size=13, color=P.MUTED)
        limit_lab.next_to(ax.c2p(1300, 2), UP, buff=0.06).shift(LEFT * 0.75)
        cap = layout.caption(
            f"A 5 M⊕ planet at {body['temp_k']:.0f} K, {body['radius_earth']:.1f} R⊕: "
            f"{body['flux_500_mjy']:.0f} mJy at 500 AU, falling as 1/d²", font_size=21)
        self.play(Create(ax), FadeIn(labels), run_time=1.0)
        self.play(Create(limit), FadeIn(limit_lab), FadeIn(cap), run_time=0.8)

        curves, drops, rows = VGroup(), VGroup(), []
        for key, col in SURVEYS:
            s = d[key]
            curves.add(clipped(ax, dist, s["snr_5"], 10, col))
            reach = s["reach_5"]
            drops.add(VGroup(
                DashedLine(ax.c2p(reach, 2), ax.c2p(reach, 0), color=col, stroke_width=2),
                Dot(ax.c2p(reach, 2), radius=0.06, color=col)))
            rows.append((col, f"{s['label']}  {s['flux_limit_mjy']:.1f} mJy at "
                              f"{s['nu_ghz']:.0f} GHz:  {reach:.0f} AU"))
        key = line_key(rows)
        head = layout.label("depth, and reach for 5 M⊕", font_size=15, color=P.FG,
                            weight="BOLD")
        legend = VGroup(head, key).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        legend.move_to([4.55, 1.6, 0])
        note = layout.label(
            timing.wrap(
                f"paper: SO {pub['so_shallow_5']:.0f} to {pub['so_deep_5']:.0f} AU, "
                f"ACT {pub['act_5_min']:.0f} to {pub['act_5_max']:.0f} AU. The SO depths "
                "here are set from the quoted reach; the ACT depth is its published median.",
                width=38),
            font_size=14, color=P.FG, line_spacing=0.9).set_opacity(0.75)
        note.next_to(legend, DOWN, buff=0.35, aligned_edge=LEFT)
        for k in range(3):
            self.play(Create(curves[k]), run_time=0.8)
            self.play(Create(drops[k]), FadeIn(key[k]), *([FadeIn(head)] if k == 0 else []),
                      run_time=0.6)
        self.play(FadeIn(note))
        timing.hold_to_read(self, cap, legend, note, settle=0.8)
        self.play(FadeOut(VGroup(ax, labels, limit, limit_lab, curves, drops, legend, note,
                                 cap)))

        # 2. reach on the mass-distance plane, against the predicted planets
        mass = np.array(d["mass_earth"])
        ax2, labels2 = widgets.labeled_axes(
            [2, 15, 1], [200, 1400, 200], x_label="mass (M⊕)",
            y_label="heliocentric distance (AU)", y_rotate=True, numbers=True,
            x_length=7.6, y_length=4.2, shift_down=-0.3)
        VGroup(ax2, labels2).shift(LEFT * 2.2)
        pop = [p for p in d["population"]
               if 2 <= p["mass_earth"] <= 15 and 200 <= p["dist_au"] <= 1400]
        dots = VGroup(*[Dot(ax2.c2p(p["mass_earth"], p["dist_au"]), radius=0.022,
                            color=P.BLUE).set_opacity(0.6) for p in pop])
        cap2 = layout.caption(
            f"{len(d['population'])} predicted Planet Nines (Brown & Batygin 2021): "
            "mass and present distance", font_size=21)
        self.play(Create(ax2), FadeIn(labels2), run_time=0.9)
        self.play(FadeIn(dots, lag_ratio=0.01), FadeIn(cap2), run_time=1.2)
        timing.hold_to_read(self, cap2, settle=0.5)

        deep = d["so_deep"]["reach_au"]
        keep = [(m, min(r, 1400)) for m, r in zip(mass, deep)]
        shade = Polygon(*[ax2.c2p(m, r) for m, r in keep],
                        ax2.c2p(mass[-1], 200), ax2.c2p(mass[0], 200), stroke_width=0)
        shade.set_fill(P.TEAL, opacity=0.10)
        reach_curves, rows2 = VGroup(), []
        for key_, col in SURVEYS:
            s = d[key_]
            reach_curves.add(clipped(ax2, mass, s["reach_au"], 1400, col))
            rows2.append((col, f"{s['label']}:  {100 * s['box_fraction']:.0f}%"))
        key2 = line_key(rows2)
        box = d["reference_box"]
        frame = Polygon(ax2.c2p(box["mass_min"], box["dist_min"]),
                        ax2.c2p(box["mass_max"], box["dist_min"]),
                        ax2.c2p(box["mass_max"], box["dist_max"]),
                        ax2.c2p(box["mass_min"], box["dist_max"]),
                        color=P.FG, stroke_width=1.6).set_stroke(opacity=0.7)
        head2 = layout.label(
            timing.wrap(
                f"share of the {box['mass_min']:.0f}-{box['mass_max']:.0f} M⊕, "
                f"{box['dist_min']:.0f}-{box['dist_max']:.0f} AU box within reach", width=30),
            font_size=15, color=P.FG, weight="BOLD", line_spacing=0.9)
        legend2 = VGroup(head2, key2).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        legend2.move_to([4.55, 1.5, 0])
        marks = VGroup(
            Line(ax2.c2p(5, pub["act_5_min"]), ax2.c2p(5, pub["act_5_max"]), color=P.FG,
                 stroke_width=5).set_stroke(opacity=0.5),
            Dot(ax2.c2p(5, pub["so_shallow_5"]), radius=0.07, color=P.FG),
            Dot(ax2.c2p(5, pub["so_deep_5"]), radius=0.07, color=P.FG))
        marks_lab = layout.label(
            timing.wrap(
                f"white, at 5 M⊕: the paper's {pub['so_shallow_5']:.0f} and "
                f"{pub['so_deep_5']:.0f} AU, and the ACT range", width=36),
            font_size=14, color=P.FG, line_spacing=0.9).set_opacity(0.75)
        marks_lab.next_to(legend2, DOWN, buff=0.35, aligned_edge=LEFT)
        cap3 = layout.caption("Below each line the planet is detected: heavier planets are "
                              "larger, so seen farther", font_size=21)
        self.play(FadeIn(shade), *[Create(c) for c in reach_curves], Create(frame),
                  FadeIn(legend2), FadeOut(cap2), FadeIn(cap3), run_time=1.6)
        self.play(FadeIn(marks), FadeIn(marks_lab))
        timing.hold_to_read(self, cap3, legend2, marks_lab, settle=1.0)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, f"Simons Observatory reaches a 5 M⊕ Planet Nine out to "
                  f"{d['so_deep_reach_5']:.0f} AU.")
