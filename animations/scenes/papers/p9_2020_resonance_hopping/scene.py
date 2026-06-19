"""Resonance hopping & classification (2020).

Classifies the distant population into resonant, hopping, and scattering classes
and their fractions -- testing whether P9 leaves a resonant fingerprint.
Reproduced in p9-2020-resonance-hopping (classification).
"""
import numpy as np
from manim import AnnularSector, Create, FadeIn, Scene, Text, UP, VGroup, Write
import p9_manim as P
from p9_manim import layout, paper

CRATE = "p9-2020-resonance-hopping"


class ResonanceHopping2020(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Sorting the distant population")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))
        fracs = [("resonant", 0.45, P.PURPLE), ("hopping", 0.30, P.TEAL), ("scattering", 0.25, P.GREEN)]
        start = 0.0
        wedges = VGroup()
        for name, f, col in fracs:
            ang = f * 2 * np.pi
            w = AnnularSector(inner_radius=0.0, outer_radius=2.0, start_angle=start, angle=ang,
                              color=col, fill_opacity=0.5, stroke_width=1)
            lbl = Text(f"{name}\n{f*100:.0f}%", font_size=16, color=col).move_to(
                2.7 * np.array([np.cos(start + ang / 2), np.sin(start + ang / 2), 0]))
            wedges.add(VGroup(w, lbl))
            start += ang
        for wg in wedges:
            self.play(Create(wg[0]), FadeIn(wg[1]), run_time=0.7)
        self.play(FadeIn(layout.takeaway("A resonant sub-population would be a smoking gun for Planet Nine.")))
        self.wait(0.9)
