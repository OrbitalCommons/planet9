"""Brown & Batygin (2021) -- the MCMC orbit posterior.

A Markov-chain fit to all the constraints yields the posterior on P9's orbit; its
median (~6 M⊕, a ~ 380 AU, e ~ 0.3) seeds this workspace's reference population.
Reproduced in p9-2021-orbit.
"""
import numpy as np
from manim import Axes, Create, DOWN, Dot, FadeIn, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper

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
        rng = np.random.default_rng(2021)
        pts = VGroup(*[Dot(ax.c2p(np.clip(380 + rng.normal(0, 45), 200, 600),
                                  np.clip(0.3 + rng.normal(0, 0.06), 0.02, 0.58)),
                           radius=0.03, color=P.TEAL).set_opacity(0.5) for _ in range(160)])
        self.play(Create(pts, lag_ratio=0.003), run_time=2.0)
        med = Dot(ax.c2p(380, 0.3), radius=0.08, color=P.GREEN)
        self.play(FadeIn(med), FadeIn(Text("median: 6.2 M⊕, 380 AU, e=0.3", font_size=15, color=P.GREEN).next_to(med, UP, buff=0.15)))
        self.play(FadeIn(layout.takeaway("Every survey-reach figure in this project samples from this posterior.")))
        self.wait(0.9)
