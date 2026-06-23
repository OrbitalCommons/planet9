"""Siraj et al. (2024) -- a refined Planet Nine orbit estimate.

Combines the dynamical and survey constraints into a posterior on P9's orbit
(mass, a, e, i), narrowing where to point. Reproduced in p9-2024-siraj-orbit.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, UR, Write, Dot, VGroup
import p9_manim as P
from p9_manim import dataio, layout, paper

CRATE = "p9-2024-siraj-orbit"


class SirajOrbit2024(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Narrowing the posterior")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[2, 14, 4], y_range=[200, 800, 200], x_length=9.0, y_length=4.0,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.4)
        self.play(Create(ax),
                  FadeIn(Text("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("semi-major axis (AU)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))

        rng = np.random.default_rng(2024)
        # a posterior cloud around ~(6 M_earth, 460 AU)
        pts = VGroup(*[Dot(ax.c2p(np.clip(6 + rng.normal(0, 1.6), 2, 14),
                                  np.clip(460 + rng.normal(0, 70), 200, 800)),
                           radius=0.03, color=P.TEAL).set_opacity(0.55) for _ in range(150)])
        self.play(Create(pts, lag_ratio=0.003), run_time=2.0)

        orbits = dataio.section("nominal_orbits")
        if orbits:
            # The three published P9 orbit fits the posterior narrows toward.
            short = {"Batygin & Brown 2016": "BB2016",
                     "Batygin et al. 2019 (review)": "review2019",
                     "Brown & Batygin 2021 (MCMC)": "MCMC2021"}
            marks = VGroup()
            labels = VGroup()
            for o in orbits:
                m = float(np.clip(o["mass_earth"], 2, 14))
                a = float(np.clip(o["a"], 200, 800))
                d = Dot(ax.c2p(m, a), radius=0.08, color=P.GREEN)
                tag = short.get(o["name"], o["name"])
                lab = Text(f"{tag}: {o['mass_earth']:.1f} M⊕, {o['a']:.0f} AU",
                           font_size=14, color=P.GREEN).next_to(d, UP, buff=0.08)
                marks.add(d)
                labels.add(lab)
            self.play(FadeIn(marks), FadeIn(labels))
        else:
            best = Dot(ax.c2p(6, 460), radius=0.08, color=P.GREEN)
            self.play(FadeIn(best))
            self.add(Text("best estimate", font_size=15, color=P.GREEN).next_to(best, UP, buff=0.1))

        eq = layout.equation_card(
            r"P(m,\,a,\,e,\,i \mid \mathcal{D})", scale=0.85).to_corner(UR, buff=0.5)
        layout.show_equation(self, eq)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Folding every constraint together points to a few-Earth-mass body near ~450 AU.")
