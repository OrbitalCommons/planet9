"""Pichierri & Batygin (2025) -- measuring the degree of clustering and diffusion.

Quantifies both how clustered the ETNOs are AND how fast that clustering diffuses
away, turning the qualitative picture into a measurable, time-dependent statistic.
Reproduced in p9-2025-clustering.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, LEFT, Scene, UP, Write, rate_functions, ValueTracker, always_redraw, Dot
import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2025-clustering"


class Clustering2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "How fast does alignment fade?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = widgets.axes([0, 4, 1], [0, 1.0, 0.5], x_length=9.0, y_length=3.8, font_size=16, shift_down=0.5)
        self.play(Create(ax),
                  FadeIn(layout.label("time (Gyr)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(layout.label("clustering strength  R̄", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        decay = ax.plot(lambda t: 0.85 * np.exp(-t / 2.2), x_range=[0, 4, 0.03], color=P.GREEN)
        self.play(Create(decay), run_time=1.5)
        # the clustering strength relaxes exponentially with a diffusion timescale
        eq = layout.explain_equation(
            self,
            [r"\bar R(t)", r"\sim", r"\bar R_0", r"\,e^{-t/\tau}"],
            [
                (0, "clustering strength as a function of time"),
                (2, "its initial (present-day) value"),
                (3, "exponential decay over a diffusion timescale tau"),
            ],
            scale=0.9,
        )
        self.play(FadeOut(eq))
        t = ValueTracker(0.0)
        dot = always_redraw(lambda: Dot(ax.c2p(t.get_value(), 0.85 * np.exp(-t.get_value() / 2.2)), color=P.TEAL, radius=0.07))
        self.add(dot)
        self.play(t.animate.set_value(4.0), run_time=3.0, rate_func=rate_functions.linear)
        layout.show_takeaway(self, "Clustering diffuses with time -- so seeing it today demands a maintaining force.")
