"""Brown & Batygin (2021) -- the MCMC orbit posterior.

A Markov-chain fit to all the constraints yields the posterior on P9's orbit; its
median (~6 M⊕, a ~ 380 AU, e ~ 0.3) seeds this workspace's reference population.
Reproduced in p9-2021-orbit.
"""
import numpy as np
from manim import Axes, Create, DOWN, Dot, FadeIn, FadeOut, Scene, Text, UP, UR, VGroup, Write
import p9_manim as P
from p9_manim import dataio, layout, paper

CRATE = "p9-2021-orbit"


class Orbit2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The posterior that anchors it all")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[200, 600, 100], y_range=[0, 0.6, 0.2], x_length=9.0, y_length=4.0,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.4)
        self.play(Create(ax),
                  FadeIn(Text("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("eccentricity", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        orbits = dataio.section("nominal_orbits")
        nom = None
        if orbits:
            nom = next((o for o in orbits if "MCMC" in o["name"] or "2021" in o["name"]), None)
        med_a = float(nom["a"]) if nom else 380.0
        med_e = float(nom["e"]) if nom else 0.3
        med_m = float(nom["mass_earth"]) if nom else 6.2

        rng = np.random.default_rng(2021)
        pts = VGroup(*[Dot(ax.c2p(np.clip(med_a + rng.normal(0, 45), 200, 600),
                                  np.clip(med_e + rng.normal(0, 0.06), 0.02, 0.58)),
                           radius=0.03, color=P.TEAL).set_opacity(0.5) for _ in range(160)])
        self.play(Create(pts, lag_ratio=0.003), run_time=2.0)
        med = Dot(ax.c2p(np.clip(med_a, 200, 600), np.clip(med_e, 0.02, 0.58)),
                  radius=0.08, color=P.GREEN)
        self.play(FadeIn(med), FadeIn(Text(
            f"median: {med_m:.1f} M⊕, {med_a:.0f} AU, e={med_e:.2g}",
            font_size=15, color=P.GREEN).next_to(med, UP, buff=0.15)))

        eq = layout.explain_equation(
            self,
            [r"\hat\theta_{\rm MCMC}", r"=", r"(6.2\,M_\oplus,\ 380\,\mathrm{AU},\ 0.3)"],
            [
                (0, "the best-fit orbit from MCMC sampling"),
                (2, "≈6 Earth masses, 380 AU, eccentricity 0.3"),
            ],
            scale=0.8,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Every survey-reach figure in this project samples from this posterior.")
