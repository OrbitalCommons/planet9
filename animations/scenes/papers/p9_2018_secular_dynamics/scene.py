"""Li, Hadden, Payne & Holman (2018) -- secular dynamics of TNOs and Planet Nine.

In the eccentricity-vector plane (k, h) = e(cosΔϖ, sinΔϖ), P9's secular forcing
shifts the equilibrium off-centre (a forced eccentricity) and orbits circulate /
librate around it. Reproduced in p9-2018-secular-dynamics (free precession rate).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, Dot, FadeIn, LEFT, ParametricFunction, Scene, Text, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2018-secular-dynamics"


class SecularDynamics2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The eccentricity-vector picture")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[-0.8, 0.8, 0.4], y_range=[-0.8, 0.8, 0.4], x_length=5.6, y_length=5.6,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.3)
        xl = layout.equation(r"k = e\cos\Delta\varpi", color=P.FG).scale(0.5).next_to(ax, DOWN, buff=0.15)
        yl = layout.equation(r"h = e\sin\Delta\varpi", color=P.FG).scale(0.5).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        data = dataio.section("phase_portrait")
        if data:
            # Map the real H(e, Delta-varpi) grid into the (k, h) =
            # e(cos Delta-varpi, sin Delta-varpi) plane, coloured by H.
            e = np.asarray(data["e"])
            dv = np.deg2rad(np.asarray(data["dvarpi_deg"]))
            H = np.asarray(data["h"])
            hmin, hmax = float(np.nanmin(H)), float(np.nanmax(H))
            span = (hmax - hmin) or 1.0

            dots = VGroup()
            for i, ev in enumerate(e):
                if ev > 0.78:
                    continue
                for j, dvv in enumerate(dv):
                    frac = (H[i, j] - hmin) / span
                    k = float(ev * np.cos(dvv))
                    hh = float(ev * np.sin(dvv))
                    col = P.BLUE if frac < 0.5 else P.TEAL
                    dots.add(Dot(ax.c2p(k, hh), radius=0.022, color=col)
                             .set_opacity(0.2 + 0.55 * frac))
            self.play(FadeIn(dots), run_time=1.6)

            # The forced eccentricity = location of the H extremum (the libration
            # centre near Delta-varpi = 180), shown as the real offset.
            mask = np.abs(np.mod(np.rad2deg(dv), 360.0) - 180.0) < 60.0
            sub = H[:, mask]
            ii = int(np.nanargmax(sub) // sub.shape[1])
            ef = float(e[ii])
            kf = -ef  # anti-aligned (Delta-varpi = 180)
            forced = Dot(ax.c2p(kf, 0), color=P.GREEN, radius=0.08)
            flbl = Text(f"forced e ≈ {ef:.2f}", font_size=15, color=P.GREEN).next_to(forced, DOWN, buff=0.1)
            loops = VGroup(*[
                ParametricFunction(lambda t, r=r: ax.c2p(kf + r * np.cos(t), r * np.sin(t)),
                                   t_range=[0, 2 * np.pi], color=P.TEAL).set_stroke(opacity=0.7)
                for r in (0.1, 0.2, 0.3)])
            self.play(FadeIn(forced), FadeIn(flbl), Create(loops, lag_ratio=0.2), run_time=1.6)
        else:
            kf = -0.32
            forced = Dot(ax.c2p(kf, 0), color=P.BLUE, radius=0.07)
            flbl = Text("forced eccentricity", font_size=15, color=P.BLUE).next_to(forced, DOWN, buff=0.1)
            self.play(FadeIn(forced), FadeIn(flbl))
            loops = VGroup(*[
                ParametricFunction(lambda t, r=r: ax.c2p(kf + r * np.cos(t), r * np.sin(t)),
                                   t_range=[0, 2 * np.pi], color=P.TEAL).set_stroke(opacity=0.8)
                for r in (0.12, 0.24, 0.36)])
            self.play(Create(loops, lag_ratio=0.2), run_time=1.8)

        eq = layout.equation_card(
            r"(k,h) = e\,(\cos\Delta\varpi,\ \sin\Delta\varpi)"
            r"\qquad \rightarrow\ e_{\rm forced}"
        ).scale(0.78).to_edge(DOWN, buff=1.7)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "P9's torque sets a forced eccentricity; orbits cycle around it, anti-aligned.")
