"""Batygin, Mardling & Nesvorny (2021) -- the stability boundary of the disk.

Distant scattered-disk orbits survive the age of the solar system only outside a
perihelion/semimajor-axis boundary; inside it they are cleared. This boundary
shapes which P9 orbits are dynamically allowed. Reproduced in p9-2021-stability.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, LEFT, Line, Scene, Text, UP, VGroup, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2021-stability"


class Stability2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The edge of survival")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        data = dataio.section("stability")
        if data:
            a = np.asarray(data["a_au"])
            q = np.asarray(data["q_crit_au"])
            amin, amax = float(a.min()), float(a.max())
            qmax = float(np.nanmax(q))
            ytop = float(np.ceil(qmax / 20.0) * 20.0)
            ax = widgets.axes([amin, amax, (amax - amin) / 4], [0, ytop, ytop / 3], x_length=9.5, y_length=4.0, font_size=18, shift_down=0.4)
            xl = layout.label("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
            yl = layout.label("perihelion q (AU)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, UP, buff=0.05)
            self.play(Create(ax), FadeIn(xl), FadeIn(yl))

            pts = [ax.c2p(float(av), float(qv)) for av, qv in zip(a, q)]
            boundary = VGroup()
            for k in range(len(pts) - 1):
                boundary.add(Line(pts[k], pts[k + 1], color=P.RED, stroke_width=4))
            self.play(Create(boundary), run_time=1.6)

            # Place region labels relative to the real boundary.
            a_lo = amin + 0.25 * (amax - amin)
            q_lo = float(np.interp(a_lo, a, q))
            a_hi = amin + 0.65 * (amax - amin)
            q_hi = float(np.interp(a_hi, a, q))
            unstable = layout.label("cleared (unstable)", font_size=18, color=P.RED).move_to(
                ax.c2p(a_lo, max(q_lo * 0.45, ytop * 0.12)))
            stable = Text("survives", font_size=20, color=P.TEAL, weight="BOLD").move_to(
                ax.c2p(a_hi, min(q_hi + 0.45 * (ytop - q_hi), ytop * 0.92)))
        else:
            ax = widgets.axes([200, 1000, 200], [30, 120, 30], x_length=9.5, y_length=4.0, font_size=18, shift_down=0.4)
            xl = layout.label("semi-major axis (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.2)
            yl = layout.label("perihelion q (AU)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, UP, buff=0.05)
            self.play(Create(ax), FadeIn(xl), FadeIn(yl))
            boundary = ax.plot(lambda a: 40 + 0.045 * (a - 200), x_range=[200, 1000, 5], color=P.RED)
            self.play(Create(boundary), run_time=1.4)
            unstable = layout.label("cleared (unstable)", font_size=18, color=P.RED).move_to(ax.c2p(420, 45))
            stable = Text("survives", font_size=20, color=P.TEAL, weight="BOLD").move_to(ax.c2p(700, 100))
        self.play(FadeIn(unstable), FadeIn(stable))
        timing.hold_to_read(self, unstable, stable, settle=0.5)

        eq = layout.explain_equation(
            self,
            [
                r"q_{\rm crit}(a)", "=", "a_N",
                r"\sqrt{\ln\!\left[\tfrac{576}{5}\tfrac{m_N}{M_\odot}\!\left(\tfrac{a}{a_N}\right)^{5/2}\right]}",
            ],
            [
                (0, "the smallest perihelion that survives, vs a"),
                (2, "scaled by Planet Nine's semimajor axis"),
                (3, "grows slowly with P9's mass and the orbit size"),
            ],
            scale=0.8,
        )
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Only detached-enough orbits last 4.5 Gyr -- which bounds the allowed P9 family.")
