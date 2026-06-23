"""Preface F -- Indirect fingerprints, and the state of the search.

Scenes: P12Indirect, P13Map.
"""
import numpy as np
from manim import (
    Arc,
    Arrow,
    Axes,
    Circle,
    Create,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    LEFT,
    Line,
    Polygon,
    RIGHT,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, orbits, dataio, timing


class P12Indirect(Scene):
    """Three ways a planet you cannot see still gives itself away."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 12"))
        title = Text("Fingerprints of an unseen planet", color=P.FG, font_size=30, weight="BOLD")
        title.to_edge(UP, buff=0.55)
        self.play(Write(title))

        col_x = [-4.4, 0.0, 4.4]

        # (1) spacecraft ranging
        c1 = VGroup()
        sun1 = Dot(LEFT * 4.4 + DOWN * 0.2, radius=0.1, color=P.SUN)
        sat = Dot(LEFT * 4.4 + UP * 1.0, radius=0.07, color=P.ORANGE)
        wig = Line(sun1.get_center(), sat.get_center(), color=P.ORANGE, stroke_width=2)
        c1lbl = Text("Cassini ranging\n(Earth–Saturn distance)", font_size=15, color=P.ORANGE)
        c1lbl.next_to(sat, UP, buff=0.2)
        self.play(FadeIn(sun1), FadeIn(sat), Create(wig), FadeIn(c1lbl))

        # (2) anomalous precession
        ell = orbits.ellipse_orbit(0.9, 0.5, color=P.BLUE, varpi=0.0).shift(DOWN * 0.2)
        ell2 = orbits.ellipse_orbit(0.9, 0.5, color=P.BLUE, varpi=0.5).shift(DOWN * 0.2).set_stroke(opacity=0.4)
        prec_lbl = layout.equation(r"\dot{\varpi}_{\rm anom}", color=P.ORANGE).scale(0.7).shift(UP * 1.2)
        self.play(Create(ell), Create(ell2), FadeIn(prec_lbl))

        # (3) solar obliquity
        sun3 = Dot(RIGHT * 4.4 + DOWN * 0.2, radius=0.32, color=P.SUN)
        axis = Line(RIGHT * 4.4 + DOWN * 1.0, RIGHT * 4.4 + UP * 0.9, color=P.FG, stroke_width=2)
        axis.rotate(np.deg2rad(6), about_point=sun3.get_center())
        vert = Line(RIGHT * 4.4 + DOWN * 1.0, RIGHT * 4.4 + UP * 0.9, color=P.MUTED, stroke_width=1).set_stroke(opacity=0.5)
        obl_lbl = Text("Sun tilted ~6°", font_size=15, color=P.ORANGE).next_to(sun3, DOWN, buff=0.5)
        self.play(FadeIn(sun3), Create(vert), Create(axis), FadeIn(obl_lbl))
        self.wait(0.4)

        # the ranging perturbation scale and the solar obliquity it implies
        sig_eq = layout.explain_equation(
            self,
            [r"\Delta\rho", r"\propto", r"\dfrac{GM_9}{d^3}",
             r"\qquad", r"\varepsilon_\odot", r"\approx", r"6^\circ"],
            [
                (0, "the tiny range tug a spacecraft would feel"),
                (2, "Planet Nine's pull (GM₉), falling off as 1/distance³"),
                (4, "the Sun's tilt versus the planets"),
                (6, "about 6 degrees -- a possible fingerprint"),
            ],
            scale=0.72,
            where=DOWN * 1.4,
        )
        self.play(FadeOut(sig_eq))

        layout.show_takeaway(
            self, "Even without a photo, dynamics constrain where (and how massive) Planet Nine can be.")


class P13Map(Scene):
    """The state of the search: the mass x distance viable region (bridges to the film)."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 13"))
        title = Text("Where could it still be?", color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.55)
        self.play(Write(title))

        data = dataio.viability()
        if data:
            mass = np.array(data["mass_earth"])
            floor = np.array([v if v else np.nan for v in data["allsky_envelope_au"]])
            deep = np.array([v if v else np.nan for v in data["detection_envelope_au"]])
            noms = [(o["mass_earth"], o["current_distance_au"], o["name"]) for o in data["nominal_orbits"]]
        else:  # documented fallback (matches the p9-viability result)
            mass = np.linspace(1, 20, 40)
            floor = 560 + 8.5 * mass
            deep = 950 + 13 * mass
            noms = [(10, 628, "BB 2016"), (5, 532, "review 2019"), (6.2, 404, "MCMC 2021")]

        ax = Axes(x_range=[1, 20, 5], y_range=[100, 2000, 400], x_length=9.5, y_length=4.4,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.5)
        xlab = Text("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        ylab = Text("distance (AU)", font_size=18, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.15)
        self.play(Create(ax), FadeIn(xlab), FadeIn(ylab))

        def poly_between(low, high, color, opacity):
            pts = []
            xs = mass
            for x, y in zip(xs, np.clip(low, 100, 2000)):
                pts.append(ax.c2p(x, y))
            for x, y in zip(xs[::-1], np.clip(high, 100, 2000)[::-1]):
                pts.append(ax.c2p(x, y))
            return Polygon(*pts, color=color, fill_opacity=opacity, stroke_width=0)

        ruled = poly_between(np.full_like(mass, 100.0), floor, P.RED, 0.22)
        partial = poly_between(floor, deep, P.ORANGE, 0.16)
        viable = poly_between(deep, np.full_like(mass, 2000.0), P.TEAL, 0.16)
        self.play(FadeIn(ruled), FadeIn(partial), FadeIn(viable))

        floor_line = ax.plot_line_graph(mass, np.clip(floor, 100, 2000), line_color=P.RED,
                                        add_vertex_dots=False, stroke_width=3)
        self.play(Create(floor_line))

        t_ruled = Text("ruled out (any sky)", font_size=16, color=P.RED).move_to(ax.c2p(6, 280))
        t_viable = Text("still viable", font_size=18, color=P.TEAL, weight="BOLD").move_to(ax.c2p(7, 1600))
        self.play(FadeIn(t_ruled), FadeIn(t_viable))

        for m, d, name in noms:
            dot = Dot(ax.c2p(m, d), radius=0.07, color=P.GREEN)
            lbl = Text(name, font_size=13, color=P.GREEN).next_to(dot, UP, buff=0.08)
            self.play(FadeIn(dot), FadeIn(lbl), run_time=0.5)

        layout.show_takeaway(
            self, "Rubin/LSST will sweep most of the viable box within a few years.")
