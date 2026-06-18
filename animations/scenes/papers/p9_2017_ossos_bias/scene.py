"""Shankman et al. (2017) -- OSSOS VI: striking biases in detecting large-a TNOs.

Surveys see distant objects only near perihelion and only where they point, so
an intrinsically isotropic population can appear clustered. The null must be
modelled carefully. Reproduced in p9-2017-ossos-bias (isotropic generator + bias).
"""
import numpy as np
from manim import (
    AnnularSector,
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    Scene,
    Text,
    UP,
    UR,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, orbits, paper

CRATE = "p9-2017-ossos-bias"


class OssosBias2017(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Is the clustering even real?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)

        rng = np.random.default_rng(2017)
        varpis = rng.uniform(0, 2 * np.pi, 12)  # truly ISOTROPIC parent
        swarm = orbits.etno_swarm([(2.5, 0.72, v) for v in varpis], color=P.GREEN, opacity=0.45)
        self.play(Create(swarm, lag_ratio=0.06), run_time=2.0)
        r0 = orbits.mean_resultant_length(varpis)
        lab0 = Text(f"intrinsic sample: isotropic (R̄ ≈ {r0:.2f})", font_size=18, color=P.GREEN).to_corner(UR, buff=0.4)
        self.play(Write(lab0))
        self.wait(0.3)

        # the survey can only DETECT objects whose perihelion falls in a sky wedge
        wedge = AnnularSector(inner_radius=0.2, outer_radius=3.4, angle=np.deg2rad(70),
                              start_angle=np.deg2rad(20), color=P.PURPLE)
        wedge.set_fill(P.PURPLE, opacity=0.12).set_stroke(P.PURPLE, width=2)
        wlab = Text("where/when the survey can actually see them", font_size=15, color=P.PURPLE)
        wlab.to_edge(DOWN, buff=1.5)
        self.play(FadeIn(wedge), FadeIn(wlab))

        # highlight only the detectable subset -> looks clustered
        detected = [v for v in varpis if np.deg2rad(20) <= (v % (2 * np.pi)) <= np.deg2rad(90)]
        keep = orbits.etno_swarm([(2.5, 0.72, v) for v in detected], color=P.TEAL, opacity=1.0, stroke_width=2.5)
        self.play(swarm.animate.set_stroke(opacity=0.12), Create(keep))
        rd = orbits.mean_resultant_length(detected) if len(detected) > 1 else 1.0
        lab1 = Text(f"detected subset looks clustered (R̄ ≈ {rd:.2f})", font_size=18, color=P.TEAL).to_corner(UR, buff=0.4)
        self.play(FadeOut(lab0), Write(lab1))
        self.play(FadeIn(layout.takeaway(
            "Selection can manufacture clustering -- the debate the field still argues.")))
        self.wait(0.9)
