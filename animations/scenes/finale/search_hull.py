"""Finale -- the search-hull "where to point" conclusion.

SearchReachHull shows the distance x apparent-brightness reach hull: how far
each survey depth can see a 6.2 Mearth Planet Nine. Numbers come from
p9-search-hull via figures/search_hull.json (dataio.search_hull()). The sky
strategy that follows it is in space_strategy.py.
"""
import numpy as np
from manim import (
    Circle,
    Create,
    DashedLine,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    LEFT,
    Polygon,
    Rectangle,
    RIGHT,
    Scene,
    SurroundingRectangle,
    Text,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, widgets, timing, dataio


# ---- documented fallbacks (match figures/search_hull.json) -------------------
_STUDY_COLORS = [P.BLUE, P.GREEN, P.ORANGE, P.PURPLE, P.TEAL]

_FALLBACK_ALLSKY = [
    {"name": "CRTS", "depth": 19.5, "reach_au": 379.1},
    {"name": "ZTF", "depth": 20.5, "reach_au": 477.2},
    {"name": "PS1 3pi", "depth": 21.5, "reach_au": 600.6},
]
_FALLBACK_SPACE_DEPTH = 24.5
_FALLBACK_SPACE_REACH = 1197.9
_FALLBACK_STUDY_NAMES = [
    "2016 Batygin & Brown (nominal)",
    "2016 Batygin & Brown (inclined-TNO variant)",
    "2019 Batygin et al. (review best-fit)",
    "2021 Brown & Batygin (MCMC median)",
    "2024 Siraj, Chyba & Tremaine (independent)",
]
# short legend labels keyed by study order
_STUDY_SHORT = ["2016 BB nominal", "2016 BB inclined", "2019 review",
                "2021 MCMC", "2024 Siraj"]


def _fallback_hull():
    """A smooth V(d) reflected-light curve standing in for the data."""
    d = np.linspace(80.0, 1500.0, 120)
    # V = H + 5 log10(d * Delta) ~ a + 10 log10(d) for outer-solar-system geometry
    v = -2.0 + 10.0 * np.log10(d)
    return d, v


def _fallback_clouds():
    """Five synthetic posterior clouds: an arc across RA with a distance ladder."""
    rng = np.random.default_rng(9)
    medians = [836.0, 665.0, 501.0, 392.0, 301.0]  # per-study current distance
    clouds = []
    for name, med in zip(_FALLBACK_STUDY_NAMES, medians):
        n = 60
        ra = (rng.normal(180.0, 70.0, n)) % 360.0
        dec = 9.0 + rng.normal(0.0, 12.0, n)
        dist = np.clip(med + rng.normal(0.0, med * 0.18, n), 120.0, 1500.0)
        vmag = -2.0 + 10.0 * np.log10(dist) + rng.normal(0.0, 0.4, n)
        samples = [{"ra_deg": float(a), "dec_deg": float(c),
                    "dist_au": float(s), "v_mag": float(v)}
                   for a, c, s, v in zip(ra, dec, dist, vmag)]
        clouds.append({"name": name, "samples": samples})
    return clouds


class SearchReachHull(Scene):
    """Distance x apparent brightness: how deep a survey must see, by distance."""

    def construct(self):
        self.add(layout.concept_badge("FINALE"))
        title = Text("Where could it still be -- and could we see it?",
                     color=P.FG, font_size=32, weight="BOLD").to_edge(UP, buff=0.55)
        self.play(Write(title))

        data = dataio.search_hull()
        if data:
            hull = data["hull"]
            dist = np.array(hull["distance_au"])
            vcurve = np.array(hull["v_curve"])
            allsky = hull["allsky_depths"]
            space_depth = hull["space_depth"]
            space_reach = hull["space_reach_au"]
            clouds = data["study_clouds"]
        else:  # documented fallback
            dist, vcurve = _fallback_hull()
            allsky = _FALLBACK_ALLSKY
            space_depth = _FALLBACK_SPACE_DEPTH
            space_reach = _FALLBACK_SPACE_REACH
            clouds = _fallback_clouds()

        ax = widgets.axes([80, 1500, 200], [16, 26, 2], x_length=9.6, y_length=4.2,
                          font_size=16, shift_down=0.35)
        xlab = layout.label("heliocentric distance (AU)", font_size=18,
                            color=P.FG).next_to(ax, DOWN, buff=0.25)
        ylab = layout.label("apparent magnitude (V)", font_size=18,
                            color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.12)
        fainter = layout.label("fainter v", font_size=14, color=P.MUTED)
        fainter.next_to(ax.c2p(80, 26), RIGHT, buff=0.1).shift(DOWN * 0.1)
        self.play(Create(ax), FadeIn(xlab), FadeIn(ylab), FadeIn(fainter))

        # deepest all-sky depth defines the "already searched (bright)" band
        deepest = max(d["depth"] for d in allsky)
        searched = Polygon(
            ax.c2p(80, 16), ax.c2p(1500, 16),
            ax.c2p(1500, deepest), ax.c2p(80, deepest),
            color=P.MUTED, fill_opacity=0.14, stroke_width=0,
        )
        searched_lbl = layout.label("already searched (all-sky)", font_size=15,
                                    color=P.MUTED).move_to(ax.c2p(1100, 18.4))
        self.play(FadeIn(searched), FadeIn(searched_lbl))

        # fiducial reflected-light V(d) curve
        m = (dist >= 80) & (dist <= 1500) & (vcurve >= 16) & (vcurve <= 26)
        curve = ax.plot_line_graph(dist[m], vcurve[m], line_color=P.FG,
                                   add_vertex_dots=False, stroke_width=3)
        curve_lbl = layout.label("V(d) for 6.2 Mearth", font_size=14,
                                 color=P.FG).move_to(ax.c2p(1280, 24.3))
        self.play(Create(curve), FadeIn(curve_lbl))
        timing.hold_to_read(self, curve_lbl, settle=0.3)

        # per-study posterior clouds, one beat each, with a compact legend
        legend_rows = VGroup()
        for i, study in enumerate(clouds[:5]):
            color = _STUDY_COLORS[i % len(_STUDY_COLORS)]
            sub = study["samples"][::25]
            dots = VGroup()
            for s in sub:
                d, v = s["dist_au"], s["v_mag"]
                if 80 <= d <= 1500 and 16 <= v <= 26:
                    dots.add(Dot(ax.c2p(d, v), radius=0.026, color=color,
                                 fill_opacity=0.8))
            chip = Dot(radius=0.06, color=color)
            name = _STUDY_SHORT[i] if i < len(_STUDY_SHORT) else study["name"]
            txt = layout.label(name, font_size=13, color=P.FG).next_to(chip, RIGHT, buff=0.12)
            row = VGroup(chip, txt)
            legend_rows.add(row)
            legend_rows.arrange(DOWN, aligned_edge=LEFT, buff=0.1)
            legend_rows.to_corner(UP + RIGHT, buff=0.35).shift(DOWN * 1.1)
            self.play(FadeIn(dots, lag_ratio=0.02), FadeIn(row), run_time=0.55)

        # survey depth lines (all-sky, orange dashed) + space depth (green solid)
        depth_group = VGroup()
        for d in allsky:
            y = d["depth"]
            ln = DashedLine(ax.c2p(80, y), ax.c2p(1500, y), color=P.ORANGE,
                            stroke_width=2, dash_length=0.12)
            lbl = layout.label(f"{d['name']} ({y:.1f})", font_size=13,
                               color=P.ORANGE).next_to(ax.c2p(160, y), UP, buff=0.04)
            depth_group.add(ln, lbl)
        self.play(*[Create(ln) for ln in depth_group if isinstance(ln, DashedLine)],
                  *[FadeIn(lbl) for lbl in depth_group if not isinstance(lbl, DashedLine)])

        space_line = ax.plot_line_graph([80, 1500], [space_depth, space_depth],
                                        line_color=P.GREEN, add_vertex_dots=False,
                                        stroke_width=3)
        space_lbl = layout.label(f"space telescope (V={space_depth:.1f})",
                                 font_size=15, color=P.GREEN, weight="BOLD")
        space_lbl.next_to(ax.c2p(420, space_depth), UP, buff=0.06)
        self.play(Create(space_line), FadeIn(space_lbl))

        # vertical guide at the space reach
        reach = DashedLine(ax.c2p(space_reach, 16), ax.c2p(space_reach, 26),
                           color=P.GREEN, stroke_width=2, dash_length=0.1)
        reach_lbl = layout.label(f"space reach ~{round(space_reach / 100) * 100:.0f} AU",
                                 font_size=14, color=P.GREEN)
        reach_lbl.next_to(ax.c2p(space_reach, 16.4), LEFT, buff=0.1)
        self.play(Create(reach), FadeIn(reach_lbl))
        timing.hold_to_read(self, reach_lbl, settle=0.4)

        eq = layout.explain_equation(
            self,
            [r"F", r"\propto", r"d^{-4}", r"\Rightarrow",
             r"V = H + 5\log_{10}(d\,\Delta)"],
            [
                (0, "reflected sunlight we'd receive"),
                (2, "it falls off as 1/distance^4 -- distant means faint"),
                (4, "so a survey's depth sets a hard maximum distance"),
                (4, "V=24.5 reaches ~1200 AU -- past every all-sky survey"),
            ],
            color=P.FG, scale=0.7, where=ax.c2p(790, 24.8),
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self,
            "All-sky surveys reach ~600 AU; V=24.5 from space reaches ~1200 AU.")
