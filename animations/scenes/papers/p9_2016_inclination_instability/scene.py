"""Madigan & McCourt (2016) -- the inclination instability.

A flattened disk of eccentric orbits is dynamically unstable: self-gravity drives
the orbits to spontaneously tilt and cone out to a common inclination, even with
NO external planet. Reproduced in p9-2016-inclination-instability (growth rates).
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, Scene, Text, UP, VGroup, Write, rate_functions,
    ValueTracker, always_redraw,
)

import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2016-inclination-instability"


class InclinationInstability2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A disk that tips itself over")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)
        n = 8
        varpis = np.linspace(0, 2 * np.pi, n, endpoint=False)
        tilt = ValueTracker(0.02)

        def disk():
            g = VGroup()
            for v in varpis:
                e = orbits.ellipse_orbit(2.6, 0.6, color=P.GREEN, varpi=v, stroke_width=1.6, opacity=0.8)
                # squash vertically to fake an out-of-plane tilt growing with `tilt`
                s = tilt.get_value()
                e.apply_matrix(np.array([[1, 0, 0], [0, max(1 - s, 0.15), 0], [0, 0, 1]]))
                g.add(e)
            return g

        live = always_redraw(disk)
        self.add(live)
        self.play(FadeIn(Text("flat, eccentric disk -- no planet", font_size=18, color=P.MUTED).to_edge(UP, buff=1.5)))
        self.play(tilt.animate.set_value(0.85), run_time=4.0, rate_func=rate_functions.smooth)
        self.play(FadeIn(layout.takeaway(
            "Self-gravity alone can tilt a disk -- a caution: not every signal needs Planet Nine.")))
        self.wait(0.9)
