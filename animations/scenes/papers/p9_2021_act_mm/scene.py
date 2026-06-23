"""(2021) -- ACT millimetre search for Planet Nine.

The Atacama Cosmology Telescope's mm maps can be searched for a moving
point-source; non-detection constrains warm/near solutions. Reproduced in
p9-2021-act-mm (mm detectability).
"""
import numpy as np
from manim import (
    Axes, Create, DOWN, FadeIn, FadeOut, Scene, Text, UP, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2021-act-mm"


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


class ActMm2021(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Hunting in the CMB maps")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = Axes(x_range=[2, 20, 4], y_range=[0, 800, 200], x_length=9.0, y_length=3.9,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
        ax.shift(DOWN * 0.5)
        xl = Text("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = Text("ACT reach (AU)", font_size=15, color=P.FG).next_to(ax, UP, buff=0.05)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        pts = [(m, r) for m, r in _reach_points("ACT (mm)") if 2 <= m <= 20]
        if pts:
            reach = ax.plot_line_graph(
                [m for m, _ in pts], [r for _, r in pts],
                line_color=P.PURPLE, add_vertex_dots=False)["line_graph"]
        else:
            reach = ax.plot(lambda m: 430 + 24 * m, x_range=[2, 20, 0.5], color=P.PURPLE)
        self.play(Create(reach), run_time=1.3)
        note = Text("229 GHz, ~8 mJy point-source limit", font_size=16, color=P.PURPLE).next_to(ax, UP, buff=0.4)
        self.play(FadeIn(note))
        timing.hold_to_read(self, note, settle=0.6)

        eq = layout.explain_equation(
            self,
            [r"F_\nu", r"\approx", r"\dfrac{2k T \nu^2}{c^2}", r"\dfrac{\pi R^2}{d^2}",
             r"\ (\text{Rayleigh--Jeans})"],
            [(0, "millimetre flux we receive"),
             (2, "long-wavelength limit of the blackbody"),
             (3, "size R over distance d, squared")],
            scale=0.78, where=DOWN * 1.6)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "CMB telescopes double as outer-solar-system surveys at millimetre wavelengths.")
