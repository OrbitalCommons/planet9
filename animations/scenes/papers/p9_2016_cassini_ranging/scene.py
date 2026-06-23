"""Fienga et al. (2016) -- Cassini ranging constraint.

Planet Nine tugs on Saturn, changing the precisely-measured Earth-Saturn range.
The residual depends on where P9 sits on its orbit; the data most favour a true
anomaly near ~108 deg. Reproduced in p9-2016-cassini-ranging.
"""
import numpy as np
from manim import (
    Axes,
    Create,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    MathTex,
    RIGHT,
    Scene,
    Text,
    UP,
    Write,
    always_redraw,
    rate_functions,
    ValueTracker,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, dataio, timing

CRATE = "p9-2016-cassini-ranging"


class CassiniRanging2016(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Weighing a planet through Saturn")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # residual amplitude vs P9 true anomaly, minimum near 108 deg
        ax = Axes(x_range=[0, 360, 90], y_range=[0, 1.1, 0.5], x_length=9.5, y_length=3.6,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.6)
        xl = Text("Planet Nine true anomaly ν (deg)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("Earth–Saturn range residual", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, RIGHT * -1, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        data = dataio.section("cassini")
        if data:
            # real Earth–Saturn range residual vs P9 true anomaly
            nu = np.asarray(data["nu_deg"], dtype=float)
            amp = np.asarray(data["amplitude_km"], dtype=float)
            ymax = float(amp.max()) or 1.0
            yn = amp / ymax  # normalise to own max for display
            curve = ax.plot_line_graph(nu, yn, add_vertex_dots=False, line_color=P.ORANGE)
            self.play(Create(curve), run_time=1.5)

            fav = float(data["favored_nu_deg"])
            fav_y = float(np.interp(fav, nu, yn))
            favored = ax.get_vertical_line(ax.c2p(fav, fav_y), color=P.TEAL, stroke_width=2)
            dot = Dot(ax.c2p(fav, fav_y), color=P.TEAL, radius=0.07)
            tag = Text(f"favoured ν ≈ {fav:.0f}°", font_size=18, color=P.TEAL).next_to(dot, UP, buff=0.15)
            self.play(Create(favored), FadeIn(dot), FadeIn(tag))
        else:
            # fallback: hand-coded residual curve with a minimum near nu=108 deg
            def resid(nu):
                return 0.25 + 0.75 * (np.sin(np.deg2rad((nu - 108) / 2.0))) ** 2

            curve = ax.plot(resid, x_range=[0, 360, 2], color=P.ORANGE)
            self.play(Create(curve), run_time=1.5)

            favored = ax.get_vertical_line(ax.c2p(108, resid(108)), color=P.TEAL, stroke_width=2)
            dot = Dot(ax.c2p(108, resid(108)), color=P.TEAL, radius=0.07)
            tag = Text("favoured ν ≈ 108°", font_size=18, color=P.TEAL).next_to(dot, UP, buff=0.15)
            self.play(Create(favored), FadeIn(dot), FadeIn(tag))
        self.wait(0.4)

        ro = paper.result_readout("least disturbance to Cassini at", "ν ≈ 108°–129°", color=P.ORANGE)
        ro.scale(0.8).to_corner(UP + RIGHT, buff=0.4)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro)
        self.play(FadeOut(ro))

        eq = layout.explain_equation(
            self,
            [r"\Delta\rho", r"\propto", r"GM_9",
             r"\left[\dfrac{\mathbf r_9-\mathbf r_S}{|\mathbf r_9-\mathbf r_S|^3} - \dfrac{\mathbf r_9}{r_9^3}\right]"],
            [
                (0, "shift in the measured Earth–Saturn range"),
                (2, "P9's gravitational pull (its mass)"),
                (3, "tidal tug on Saturn minus the tug on the Sun"),
            ],
            scale=0.82,
            where=DOWN * 0.6)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self,
            "Spacecraft ranging both constrains the mass and points to where P9 hides.")
