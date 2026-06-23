"""(2019) OSSOS scattering -- producing detached objects with Planet Nine.

P9 lifts scattering-disk objects above the Neptune-scattering line into stable
detached orbits, populating the high-q region the data show. Reproduced in
p9-2019-ossos-scattering (boundary, population).
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, UR, VGroup, Write, Dot
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2019-ossos-scattering"


class OssosScattering2019(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Lifting orbits out of reach")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        ax = Axes(x_range=[150, 700, 100], y_range=[30, 90, 20], x_length=9.3, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.4)
        self.play(Create(ax),
                  FadeIn(Text("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)),
                  FadeIn(Text("perihelion q (AU)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)))
        line = ax.plot(lambda a: 38.0, x_range=[150, 700, 5], color=P.RED)
        self.play(Create(line))
        self.add(Text("Neptune-scattering line", font_size=14, color=P.RED).next_to(ax.c2p(560, 38), UP, buff=0.05))
        rng = np.random.default_rng(2019)
        detached = VGroup(*[Dot(ax.c2p(rng.uniform(250, 650), rng.uniform(45, 85)), radius=0.05, color=P.GREEN)
                            for _ in range(20)])
        self.play(Create(detached, lag_ratio=0.03))
        self.add(Text("detached (lifted by P9)", font_size=15, color=P.GREEN).move_to(ax.c2p(400, 78)))

        eq = layout.explain_equation(
            self,
            [r"q", r"=", r"a(1-e)", r"\gtrsim", r"40\ \mathrm{AU}", r"\Rightarrow", r"\text{detached}"],
            [
                (0, "perihelion — closest approach to the Sun"),
                (2, "set by the orbit's size and stretch"),
                (4, "lifted well beyond Neptune"),
                (6, "no longer scattered by Neptune"),
            ],
            scale=0.75,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "P9 ferries objects above the scattering line into the detached belt.")
