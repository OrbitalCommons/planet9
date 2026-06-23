"""(2021) -- the DES distant-object catalogue.

The Dark Energy Survey's catalogue of trans-Neptunian detections: positions,
magnitudes, and the orbits used to test clustering and exclusion. Reproduced in
p9-2021-des-catalog.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, UR, VGroup, Write, Dot
import p9_manim as P
from p9_manim import layout, paper, widgets

CRATE = "p9-2021-des-catalog"


class DesCatalog2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "What DES actually catalogued")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        # the DES footprint patch with catalogued detections
        patch = widgets.footprint_rect(6.0, 3.2, fill_opacity=0.06, stroke_width=2).shift(DOWN * 0.2)
        self.play(Create(patch), FadeIn(Text("DES footprint", font_size=16, color=P.PURPLE).next_to(patch, UP, buff=0.1)))
        rng = np.random.default_rng(2021)
        dets = VGroup(*[Dot([rng.uniform(-2.8, 2.8), rng.uniform(-1.4, 1.4), 0], radius=0.05, color=P.GREEN)
                        for _ in range(24)])
        self.play(Create(dets, lag_ratio=0.03))
        self.add(Text("catalogued distant TNOs (g,r,i,z,Y)", font_size=15, color=P.GREEN).to_edge(DOWN, buff=1.5))

        eq = layout.explain_equation(
            self,
            [r"\text{catalogued}", r"\iff", r"m", r"<", r"m_{\rm lim}"],
            [
                (0, "an object makes it into the catalogue"),
                (2, "only if it's bright enough"),
                (4, "i.e. above the survey's depth limit"),
            ],
            scale=0.8,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "A deep, well-characterised catalogue -- the raw material for the exclusion tests.")
