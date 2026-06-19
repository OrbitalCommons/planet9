"""Pichierri & Batygin (2025) -- measuring the degree of clustering and diffusion.

Quantifies both how clustered the ETNOs are AND how fast that clustering diffuses
away, turning the qualitative picture into a measurable, time-dependent statistic.
Reproduced in p9-2025-clustering.
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, Scene, Text, UP, Write, rate_functions, ValueTracker, always_redraw, Dot
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2025-clustering"


class Clustering2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "How fast does alignment fade?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 4, 1], y_range=[0, 1.0, 0.5], x_length=9.0, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax),
                  FadeIn(Text("time (Gyr)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("clustering strength  R̄", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        decay = ax.plot(lambda t: 0.85 * np.exp(-t / 2.2), x_range=[0, 4, 0.03], color=P.GREEN)
        self.play(Create(decay), run_time=1.5)
        t = ValueTracker(0.0)
        dot = always_redraw(lambda: Dot(ax.c2p(t.get_value(), 0.85 * np.exp(-t.get_value() / 2.2)), color=P.TEAL, radius=0.07))
        self.add(dot)
        self.play(t.animate.set_value(4.0), run_time=3.0, rate_func=rate_functions.linear)
        self.play(FadeIn(layout.takeaway("Clustering diffuses with time -- so seeing it today demands a maintaining force.")))
        self.wait(0.9)
