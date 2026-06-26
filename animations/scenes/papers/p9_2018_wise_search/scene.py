"""Meisner et al. / WISE W1 shift-and-stack search.

WISE scanned the whole sky in the thermal infrared, but at 3.4 µm (W1) a ~40 K
body is deep on the Wien tail -- so WISE constrains only relatively nearby /
massive Planet Nine solutions. Reproduced in p9-2018-wise-search.
"""
import numpy as np
from manim import (
    Create, DOWN, FadeIn, FadeOut, LEFT, Scene, UP, Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2018-wise-search"


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


class WiseSearch2018(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "All-sky, but only so deep")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        ax = widgets.axes([2, 20, 4], [0, 400, 100], x_length=9.0, y_length=3.9, font_size=18, shift_down=0.5)
        xl = layout.label("mass (Earth masses)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.25)
        yl = layout.label("WISE W1 reach (AU)", font_size=16, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))

        # shallow reach ~ grows slowly with mass (cold body faint at 3.4 um)
        pts = [(m, r) for m, r in _reach_points("WISE W1") if 2 <= m <= 20]
        if pts:
            reach = ax.plot_line_graph(
                [m for m, _ in pts], [r for _, r in pts],
                line_color=P.ORANGE, add_vertex_dots=False)["line_graph"]
        else:
            reach = ax.plot(lambda m: 120 + 9.0 * m, x_range=[2, 20, 0.5], color=P.ORANGE)
        self.play(Create(reach), run_time=1.4)
        note = layout.label("a 40 K world barely emits at 3.4 µm", font_size=16, color=P.ORANGE).next_to(ax, UP, buff=0.1)
        self.play(FadeIn(note))
        timing.hold_to_read(self, note, settle=0.8)

        eq = layout.explain_equation(
            self,
            [r"F_\nu", "=", r"\pi", r"B_\nu(T)", r"\left(\tfrac{R}{d}\right)^2"],
            [(0, "thermal flux we'd receive from the planet"),
             (3, "blackbody brightness at temperature T"),
             (4, "shrinks with size R over distance d, squared")],
            scale=0.85, where=DOWN * 1.6)
        self.play(FadeOut(eq))

        layout.show_takeaway(
            self, "Whole-sky thermal coverage -- but shallow for the coldest, most distant cases.")
