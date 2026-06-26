"""Finale -- the search-hull "where to point" conclusion.

Two data-driven scenes that close the film: SearchReachHull shows the
distance x apparent-brightness reach hull (how far each survey depth can see a
6.2 Mearth Planet Nine), and WhereToPoint maps the un-searched, space-reachable
sliver of the posterior onto the sky. Numbers come from p9-search-hull via
figures/search_hull.json (dataio.search_hull()).
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


class WhereToPoint(Scene):
    """The sky map: the un-searched, space-reachable sliver of the posterior."""

    def construct(self):
        self.add(layout.concept_badge("FINALE"))
        title = Text("Where to point next", color=P.FG, font_size=32,
                     weight="BOLD").to_edge(UP, buff=0.55)
        self.play(Write(title))

        data = dataio.search_hull()
        if data:
            clouds = data["study_clouds"]
            targets = data["targets"]
            surveys = data["surveys"]
            lsst = data["lsst"]
            resid_reach = data["summary"]["residual_space_reachable_prob"]
        else:  # documented fallback
            clouds = _fallback_clouds()
            targets = [{"rank": 1, "ra_deg": 61.5, "dec_deg": 10.5,
                        "v_mean": 21.6, "dist_mean_au": 667.0,
                        "best_survey": "PS1 3pi"}]
            surveys = [
                {"name": "PS1 3pi", "depth": 21.5, "dec_min_deg": -30.0, "dec_max_deg": 90.0},
                {"name": "DES", "depth": 23.8, "dec_min_deg": -65.0, "dec_max_deg": 5.0},
            ]
            lsst = {"dec_min_deg": -75.0, "dec_max_deg": 12.0}
            resid_reach = 0.29

        ax = widgets.axes([0, 360, 60], [-60, 60, 30], x_length=10.0, y_length=4.0,
                          font_size=16, shift_down=0.35)
        xlab = layout.label("right ascension (deg)", font_size=18,
                            color=P.FG).next_to(ax, DOWN, buff=0.25)
        ylab = layout.label("declination (deg)", font_size=18,
                            color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.12)
        self.play(Create(ax), FadeIn(xlab), FadeIn(ylab))

        def band(dec_lo, dec_hi, color, opacity):
            dec_lo = max(dec_lo, -60.0)
            dec_hi = min(dec_hi, 60.0)
            if dec_hi <= dec_lo:
                return None
            return Polygon(
                ax.c2p(0, dec_lo), ax.c2p(360, dec_lo),
                ax.c2p(360, dec_hi), ax.c2p(0, dec_hi),
                color=color, fill_opacity=opacity, stroke_width=0,
            )

        # schematic survey coverage bands ("already searched")
        survey_colors = [P.BLUE, P.GREEN, P.MUTED, P.ORANGE]
        cov_group = VGroup()
        for i, s in enumerate(surveys):
            b = band(s.get("dec_min_deg", -60.0), s.get("dec_max_deg", 60.0),
                     survey_colors[i % len(survey_colors)], 0.07)
            if b is not None:
                cov_group.add(b)
        cov_lbl = layout.label("already searched (survey dec bands)", font_size=14,
                               color=P.MUTED).move_to(ax.c2p(180, -52))
        self.play(FadeIn(cov_group), FadeIn(cov_lbl))

        # P9 posterior swarm across the sky (pool a couple of studies)
        swarm = VGroup()
        for i, study in enumerate(clouds[:3]):
            for s in study["samples"][::25]:
                ra, dec = s["ra_deg"], s["dec_deg"]
                if 0 <= ra <= 360 and -60 <= dec <= 60:
                    swarm.add(Dot(ax.c2p(ra, dec), radius=0.024, color=P.TEAL,
                                  fill_opacity=0.7))
        swarm_lbl = layout.label("P9 posterior (study clouds)", font_size=14,
                                 color=P.TEAL).move_to(ax.c2p(260, 48))
        self.play(FadeIn(swarm, lag_ratio=0.01), FadeIn(swarm_lbl))
        timing.hold_to_read(self, swarm_lbl, settle=0.4)

        # Rubin/LSST footprint outline
        lo = max(lsst.get("dec_min_deg", -75.0), -60.0)
        hi = min(lsst.get("dec_max_deg", 12.0), 60.0)
        lsst_rect = Rectangle(
            width=ax.c2p(360, 0)[0] - ax.c2p(0, 0)[0],
            height=ax.c2p(0, hi)[1] - ax.c2p(0, lo)[1],
            color=P.PURPLE, stroke_width=2.5,
        ).set_fill(P.PURPLE, opacity=0.05)
        lsst_rect.move_to(ax.c2p(180, (lo + hi) / 2))
        lsst_lbl = layout.label("Rubin/LSST reach (r<=24.5)", font_size=14,
                                color=P.PURPLE, weight="BOLD").move_to(ax.c2p(70, -22))
        self.play(Create(lsst_rect), FadeIn(lsst_lbl))

        # mark top targets; #1 gets a bright ring + readout
        ring_group = VGroup()
        for t in targets[:4]:
            ra, dec = t["ra_deg"], t["dec_deg"]
            if t.get("rank", 99) == 1 or t is targets[0]:
                ring = Circle(radius=0.22, color=P.RED, stroke_width=4).move_to(ax.c2p(ra, dec))
            else:
                ring = Circle(radius=0.13, color=P.ORANGE, stroke_width=2.5).move_to(ax.c2p(ra, dec))
            ring_group.add(ring)
        top = targets[0]
        t1_lbl = layout.label(
            f"#1: RA {top['ra_deg']:.0f}, Dec {top['dec_deg']:+.0f}, "
            f"V~{top['v_mean']:.1f}, ~{round(top['dist_mean_au']):.0f} AU",
            font_size=15, color=P.RED, weight="BOLD")
        t1_lbl.next_to(ax.c2p(top["ra_deg"], top["dec_deg"]), UP, buff=0.28)
        self.play(*[Create(r) for r in ring_group], FadeIn(t1_lbl))
        timing.hold_to_read(self, t1_lbl, settle=0.5)

        # boxed headline readout
        pct = f"{resid_reach * 100:.0f}%"
        rd_lab = layout.label("un-searched AND space-reachable", font_size=16, color=P.MUTED)
        rd_val = layout.label(f"{pct} of the P9 posterior", font_size=24,
                              color=P.GREEN, weight="BOLD").next_to(rd_lab, DOWN, buff=0.1)
        rd = VGroup(rd_lab, rd_val)
        rd_box = SurroundingRectangle(rd, color=P.GREEN, buff=0.2, corner_radius=0.1)
        rd_box.set_fill(P.GREEN, opacity=0.06)
        readout = VGroup(rd_box, rd).scale(0.85).to_corner(UP + RIGHT, buff=0.4).shift(DOWN * 0.9)
        self.play(FadeIn(readout))
        timing.hold_to_read(self, rd_val, settle=0.6)

        layout.show_takeaway(
            self,
            "Point at the northern arc past PS1's depth: 29% un-searched, reachable.")
