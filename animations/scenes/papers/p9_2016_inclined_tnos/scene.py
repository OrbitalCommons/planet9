"""Batygin & Brown (2016) -- generation of highly inclined TNOs by Planet Nine.

Beyond confining apsides, P9 can pump some objects to very high inclinations --
even retrograde -- explaining the puzzling population of large-i TNOs. Reproduced
in p9-2016-inclined-tnos.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, Scene, UP, VGroup, Write, rate_functions,
    ValueTracker, always_redraw,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2016-inclined-tnos"


class InclinedTnos2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Kicking orbits out of the plane")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)
        # a few orbits start in-plane, then tilt up to high inclination
        n = 5
        varpis = np.linspace(0.2, 2 * np.pi, n, endpoint=False)
        tilt = ValueTracker(0.95)  # 1 = flat, ->0 squashed (high i seen edge-on)

        def orbset():
            g = VGroup()
            for v in varpis:
                e = orbits.ellipse_orbit(2.4, 0.65, color=P.GREEN, varpi=v, stroke_width=1.8, opacity=0.85)
                e.apply_matrix(np.array([[1, 0, 0], [0, max(tilt.get_value(), 0.1), 0], [0, 0, 1]]))
                g.add(e)
            return g

        live = always_redraw(orbset)
        self.add(live)
        msg = layout.label("initially near the ecliptic", font_size=18, color=P.MUTED).to_edge(UP, buff=1.5)
        self.play(FadeIn(msg))
        timing.hold_to_read(self, msg, settle=0.4)
        self.play(tilt.animate.set_value(0.12), run_time=3.5, rate_func=rate_functions.smooth)

        eq = layout.explain_equation(
            self,
            [r"\sqrt{1-e^2}", r"\,\cos i", r"\approx", r"\text{const}"],
            [
                (0, "a factor that shrinks as eccentricity grows"),
                (1, "times the cosine of the inclination"),
                (3, "this combination (Kozai) is roughly conserved"),
            ],
            scale=1.0,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "P9 can also lift orbits to steep -- even retrograde -- inclinations.")
