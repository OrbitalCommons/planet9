"""(2024) -- Neptune-crossing TNOs and Planet Nine.

Some distant objects dip to perihelia that cross Neptune's orbit; their survival
and supply rate test the P9 + scattering picture. Reproduced in
p9-2024-neptune-crossing (observed TNOs, simulation).
"""
import numpy as np
from manim import Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, UR, Write
import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2024-neptune-crossing"


class NeptuneCrossing2024(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Brushing past Neptune")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        from manim import Circle
        sun = orbits.sun()
        neptune = Circle(radius=1.0, color=P.BLUE, stroke_width=2).set_stroke(opacity=0.7)
        self.play(FadeIn(sun), Create(neptune))
        self.add(Text("Neptune", font_size=14, color=P.BLUE).next_to(neptune, DOWN, buff=0.05))
        # a distant eccentric orbit whose perihelion dips inside Neptune's orbit
        crosser = orbits.ellipse_orbit(3.2, 0.78, color=P.GREEN, varpi=np.deg2rad(15))
        self.play(Create(crosser))
        peri = orbits.body_dot(3.2, 0.78, 0.0, color=P.RED, radius=0.08)
        self.play(FadeIn(peri))
        self.add(Text("perihelion crosses Neptune", font_size=15, color=P.RED).to_edge(DOWN, buff=1.5))

        eq = layout.explain_equation(
            self,
            [r"q", r"=", r"a(1-e)", r"<", r"a_N\ (30\ \mathrm{AU})"],
            [
                (0, "perihelion — closest approach to the Sun"),
                (2, "set by the orbit's size and stretch"),
                (4, "dips inside Neptune's orbit, so they cross"),
            ],
            scale=0.8,
            where=UP * 2.4)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Neptune-crossers are short-lived -- their supply constrains the distant disk.")
