"""Inner Oort cloud origin (Batygin & Brown 2021 context).

Planet Nine may itself be an ice giant scattered out during the Sun's birth
cluster, parked in the inner Oort cloud (a ~ 2000-20000 AU) -- which is also a
source of detached objects. Reproduced in p9-2021-oort-cloud.
"""
import numpy as np
from manim import Create, DOWN, FadeIn, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2021-oort-cloud"


class OortCloud2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Born in a crowded nursery")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun(radius=0.12)
        rng = np.random.default_rng(2021)
        # a sparse inner-Oort population of large, eccentric, isotropic orbits
        swarm = orbits.etno_swarm([(rng.uniform(2.0, 3.4), rng.uniform(0.6, 0.85), rng.uniform(0, 2 * np.pi))
                                   for _ in range(9)], color=P.GREEN, opacity=0.6)
        self.play(FadeIn(sun), Create(swarm, lag_ratio=0.06), run_time=2.0)
        self.add(Text("inner Oort cloud: a ≈ 2000–20000 AU", font_size=17, color=P.GREEN).to_edge(UP, buff=1.5))
        p9 = orbits.ellipse_orbit(3.0, 0.4, color=P.BLUE, varpi=np.deg2rad(40), stroke_width=3)
        self.play(Create(p9))
        self.add(Text("P9: a scattered ice giant?", font_size=16, color=P.BLUE).to_edge(DOWN, buff=1.5))
        self.play(FadeIn(layout.takeaway("The Sun's birth cluster could have made -- and parked -- Planet Nine.")))
        self.wait(0.9)
