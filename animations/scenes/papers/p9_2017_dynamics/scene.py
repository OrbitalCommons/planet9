"""Batygin & Morbidelli (2017) -- Dynamical Evolution Induced by Planet Nine.

The secular (orbit-averaged) phase space has a stable island in which a distant
KBO's apsidal angle Δϖ relative to P9 librates -- i.e. stays confined -- rather
than circulating freely. Reproduced in p9-2017-dynamics (phase-space/Hamiltonian).
"""
import numpy as np
from manim import (
    Axes,
    Create,
    DOWN,
    Dot,
    FadeIn,
    LEFT,
    ParametricFunction,
    Scene,
    Text,
    UP,
    UR,
    VGroup,
    Write,
    rate_functions,
    ValueTracker,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2017-dynamics"


class Dynamics2017(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A trap in phase space")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[0, 360, 90], y_range=[0, 0.8, 0.2], x_length=9.5, y_length=4.2,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.4)
        xl = layout.equation(r"\Delta\varpi\ \text{(deg, relative to P9)}", color=P.FG).scale(0.6).next_to(ax, DOWN, buff=0.2)
        yl = Text("eccentricity", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        data = dataio.section("phase_portrait")
        if data:
            # Real secular Hamiltonian grid H(e, Delta-varpi). Wrap Delta-varpi
            # into [0, 360) so the libration island around 180 is centred.
            e = np.asarray(data["e"])
            dv = np.mod(np.asarray(data["dvarpi_deg"]), 360.0)
            h = np.asarray(data["h"])
            hmin, hmax = float(np.nanmin(h)), float(np.nanmax(h))
            span = (hmax - hmin) or 1.0

            # Scatter every grid cell, coloured by H (blue = low, teal = high).
            dots = VGroup()
            for i, ev in enumerate(e):
                if ev > 0.8:
                    continue
                for j, dvv in enumerate(dv):
                    frac = (h[i, j] - hmin) / span
                    col = P.BLUE if frac < 0.5 else P.TEAL
                    dots.add(Dot(ax.c2p(float(dvv), float(ev)), radius=0.022,
                                 color=col).set_opacity(0.18 + 0.55 * frac))
            self.play(FadeIn(dots), run_time=1.4)

            # The libration island = the H extremum near Delta-varpi = 180.
            mask = (dv > 120) & (dv < 240)
            sub = h[:, mask]
            ii = int(np.nanargmax(sub) // sub.shape[1])
            island_e = float(e[ii])
            center = Dot(ax.c2p(180, island_e), color=P.GREEN, radius=0.07)
            self.play(FadeIn(center))

            # Mark cells closest to the island extremum (iso-H island shading).
            thresh = hmax - 0.12 * span
            island = VGroup(*[
                Dot(ax.c2p(float(dv[j]), float(e[i])), radius=0.03, color=P.TEAL)
                for i in range(len(e)) for j in range(len(dv))
                if e[i] <= 0.8 and 120 < dv[j] < 240 and h[i, j] >= thresh])
            self.play(Create(island, lag_ratio=0.01), run_time=1.2)

            # A test particle librating around the real island centre.
            t = ValueTracker(0.0)
            amp_e = min(island_e, 0.8 - island_e, 0.22)
            body = always_redraw(lambda: Dot(
                ax.c2p(180 + 55 * np.cos(t.get_value()),
                       island_e + amp_e * np.sin(t.get_value())),
                color=P.GREEN, radius=0.07))
            self.add(body)
            self.play(t.animate.set_value(2.2 * np.pi), run_time=4.0,
                      rate_func=rate_functions.linear)
        else:
            def loop(scale):
                return ParametricFunction(
                    lambda t: ax.c2p(180 + 70 * scale * np.cos(t), 0.42 + 0.26 * scale * np.sin(t)),
                    t_range=[0, 2 * np.pi], color=P.TEAL)
            loops = VGroup(*[loop(s) for s in (0.35, 0.65, 1.0)])
            self.play(Create(loops, lag_ratio=0.2), run_time=1.8)
            circ = VGroup(*[ax.plot(lambda x, c=c: c, x_range=[0, 360, 5], color=P.MUTED).set_stroke(opacity=0.5)
                            for c in (0.66, 0.74)])
            self.play(Create(circ))
            center = Dot(ax.c2p(180, 0.42), color=P.BLUE, radius=0.06)
            self.play(FadeIn(center))
            t = ValueTracker(0.0)
            body = always_redraw(lambda: Dot(
                ax.c2p(180 + 60 * np.cos(t.get_value()), 0.42 + 0.22 * np.sin(t.get_value())),
                color=P.GREEN, radius=0.07))
            self.add(body)
            self.play(t.animate.set_value(2.2 * np.pi), run_time=4.0, rate_func=rate_functions.linear)

        lbl = Text("libration: Δϖ stays confined", font_size=18, color=P.TEAL).to_corner(UR, buff=0.5)
        self.play(Write(lbl))
        timing.hold_to_read(self, lbl, settle=0.5)

        eq = layout.equation_card(
            r"\langle H\rangle(e,\,\Delta\varpi):\ \text{libration about}\ \Delta\varpi = 180^\circ"
        ).scale(0.62).to_corner(UP + LEFT, buff=0.5).shift(DOWN * 0.7)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "Inside the island, orbits stay anti-aligned with P9 -- that is the clustering.")
