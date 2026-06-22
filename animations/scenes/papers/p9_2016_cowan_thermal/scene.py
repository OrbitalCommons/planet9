"""Cowan et al. (2016) -- thermal SED / detectability of Planet Nine.

Models the spectral energy distribution of a cold P9 and where it peaks, to
predict the best wavebands for a thermal detection. Reproduced in
p9-2016-cowan-thermal (SED).
"""
import numpy as np
from manim import Axes, Create, DOWN, FadeIn, Line, Scene, Text, UP, Write
import p9_manim as P
from p9_manim import layout, paper, dataio

CRATE = "p9-2016-cowan-thermal"


class CowanThermal2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Where does a cold world shine?")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 3, 1], y_range=[0, 1.1, 0.5], x_length=9.3, y_length=3.8,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        self.play(Create(ax), FadeIn(Text("wavelength (µm, log)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)))
        for lx, t in [(0, "1"), (1, "10"), (2, "100"), (3, "1000")]:
            self.add(Text(t, font_size=13, color=P.MUTED).next_to(ax.c2p(lx, 0), DOWN, buff=0.1))

        data = dataio.section("thermal")
        if data:
            # real normalised Planck spectrum of a 40 K body (log wavelength)
            wl_um = np.asarray(data["wavelength_um"], dtype=float)
            pn = np.asarray(data["planck_norm_40k"], dtype=float)
            lx = np.log10(wl_um)
            mask = (lx >= 0) & (lx <= 3)  # axis covers 1..1000 µm
            curve = ax.plot_line_graph(lx[mask], pn[mask], add_vertex_dots=False, line_color=P.YELLOW)
            self.play(Create(curve), run_time=1.4)
            wien = float(data["wien_peak_um"])
            self.add(Line(ax.c2p(np.log10(wien), 0), ax.c2p(np.log10(wien), 1.05),
                          color=P.TEAL, stroke_width=2).set_stroke(opacity=0.7))
            self.add(Text(f"Wien peak {wien:.0f} µm", font_size=13, color=P.TEAL)
                     .next_to(ax.c2p(np.log10(wien), 1.05), UP, buff=0.05))
        else:
            def sed(lx):
                wl = 10 ** lx
                x = 60.0 / wl  # ~peak near 60-90 um for ~40 K
                return float(np.clip((1.0 / wl ** 4) / (np.exp(min(x, 50)) - 1 + 1e-9) * 8e6, 0, 1.05))
            curve = ax.plot(sed, x_range=[0, 3, 0.03], color=P.YELLOW)
            self.play(Create(curve), run_time=1.4)
        for wl, name, col in [(60, "IRAS", P.RED), (90, "AKARI", P.ORANGE), (1000, "mm", P.PURPLE)]:
            self.add(Text(name, font_size=14, color=col).move_to(ax.c2p(np.log10(wl), 0.9)))
        self.play(FadeIn(layout.takeaway("A ~40 K body's SED peaks near 60–90 µm -- the far-IR is the sweet spot.")))
        self.wait(0.9)
