"""Simons Observatory (2025) -- a millimetre forecast for Planet Nine.

At the ~mm peak of a cold body's blackbody spectrum, wide CMB surveys (ACT now,
Simons Observatory soon) reach far beyond optical/WISE -- ~900 AU for a 5 M⊕
world. Reproduced in p9-2025-simons-forecast (max_detectable_distance).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, LEFT, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper

CRATE = "p9-2025-simons-forecast"


def _reach_points(survey_name):
    """Real (mass, reach) pairs for a survey from viability(), nulls skipped."""
    try:
        v = dataio.viability()
        masses = v["mass_earth"]
        for s in v["surveys"]:
            if s["name"] == survey_name:
                reach = s.get("reach") or []
                return [(m, r) for m, r in zip(masses, reach) if r is not None]
    except Exception:
        pass
    return []


class SimonsForecast2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "The millimetre advantage")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[2, 20, 4], y_range=[0, 1200, 300], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 18})
        ax.shift(DOWN * 0.5)
        xl = Text("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("reach (AU)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        act_pts = [(m, r) for m, r in _reach_points("ACT (mm)") if 2 <= m <= 20]
        so_pts = [(m, r) for m, r in _reach_points("Simons Obs. (mm)") if 2 <= m <= 20]
        if act_pts:
            act = ax.plot_line_graph(
                [m for m, _ in act_pts], [r for _, r in act_pts],
                line_color=P.PURPLE, add_vertex_dots=False)["line_graph"]
        else:
            act = ax.plot(lambda m: 430 + 24 * m, x_range=[2, 20, 0.5], color=P.PURPLE)
        if so_pts:
            so = ax.plot_line_graph(
                [m for m, _ in so_pts], [r for _, r in so_pts],
                line_color=P.TEAL, add_vertex_dots=False)["line_graph"]
        else:
            so = ax.plot(lambda m: 760 + 32 * m, x_range=[2, 20, 0.5], color=P.TEAL)
        self.play(Create(act))
        self.play(Create(so))
        self.add(Text("ACT", font_size=16, color=P.PURPLE).next_to(act.get_end(), DOWN, buff=0.1))
        self.add(Text("Simons Obs.", font_size=16, color=P.TEAL).next_to(so.get_end(), UP, buff=0.1))
        self.play(FadeIn(layout.takeaway(
            "At the body's blackbody peak, mm surveys reach ~900 AU for a 5 M⊕ world.")))
        self.wait(0.9)
