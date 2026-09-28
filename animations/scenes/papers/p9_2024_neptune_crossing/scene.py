"""Batygin, Morbidelli, Brown & Nesvorny (2024) -- low-inclination Neptune crossers.

Long-period objects whose perihelia dip inside Neptune's orbit do not last:
Neptune scatters them away. Yet 17 are known, on nearly flat orbits. Planet
Nine keeps supplying them by walking distant perihelia inward; without it the
supply dries up. The test compares each object's perihelion with the model's
prediction at the distance it was found. Reproduced in p9-2024-neptune-crossing
at reduced scale: the sample, the simulated footprints, the per-object
percentiles and the null distribution of the statistic are the crate's own
(anim.json -> papers -> p9-2024-neptune-crossing); the paper's two statistics
are placed on that null and labelled.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    UP,
    Axes,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2024-neptune-crossing"


def log_axes(centre, x_length, y_length):
    ax = widgets.axes([2.0, 4.0, 0.5], [0, 40, 10], x_length=x_length, y_length=y_length,
                      shift_down=0)
    ax.move_to(centre)
    ax.get_y_axis().add_numbers(font_size=16)
    xt = VGroup(*[layout.label(str(v), font_size=15).next_to(ax.c2p(np.log10(v), 0), DOWN,
                                                            buff=0.12)
                  for v in (100, 300, 1000, 3000, 10000)])
    xl = layout.label("semi-major axis a (AU, log scale)", font_size=17).next_to(
        xt, DOWN, buff=0.1).set_x(ax.get_center()[0])
    yl = layout.label("perihelion q (AU)", font_size=15).rotate(np.pi / 2).next_to(
        ax, LEFT, buff=0.4)
    return ax, VGroup(xt, xl, yl)


class EdgeAxes(Axes):
    """Axes that meet at the lower-left corner even when a range spans zero."""

    @staticmethod
    def _origin_shift(axis_range):
        return axis_range[0]


def edge_axes(x_range, y_range, x_label, y_label, x_length, y_length):
    ax = EdgeAxes(x_range=x_range, y_range=y_range, x_length=x_length, y_length=y_length,
                  axis_config={"color": P.MUTED, "include_tip": False, "font_size": 16})
    ax.add_coordinates()
    xl = layout.label(x_label, font_size=18).next_to(ax, DOWN, buff=0.12)
    yl = layout.label(y_label, font_size=15).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.12)
    return ax, VGroup(xl, yl)


class NeptuneCrossing2024(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        cuts = d["cuts"]
        wp, wo = d["with_p9"], d["without_p9"]
        run = d["run"]
        self.add(paper.scene_header(CRATE))

        # 1. the sample: perihelia inside Neptune's orbit
        ax, furn = log_axes([-1.9, 0.35, 0], 7.8, 4.4)
        nep = DashedLine(ax.c2p(2.0, cuts["q_max_au"]), ax.c2p(4.0, cuts["q_max_au"]),
                         color=P.ORANGE, stroke_width=2.5)
        nep_l = layout.label(f"Neptune's orbit, {cuts['q_max_au']:.0f} AU", font_size=16,
                             color=P.ORANGE).next_to(ax.c2p(3.5, cuts["q_max_au"]), UP,
                                                     buff=0.12)
        dots = VGroup(*[Dot(ax.c2p(np.log10(o["a_au"]), o["q_au"]), radius=0.075,
                            color=P.GREEN).set_z_index(5) for o in objs])
        self.play(Create(ax), FadeIn(furn), Create(nep), FadeIn(nep_l))
        cap = layout.caption(f"{len(objs)} known objects: a > {cuts['a_min_au']:.0f} AU, "
                             f"i < {cuts['i_max_deg']:.0f}°, perihelion inside Neptune's orbit",
                             font_size=22)
        self.play(FadeIn(dots, lag_ratio=0.08), FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.4)
        cap2 = layout.caption("Neptune soon scatters such orbits away: "
                              "something must keep replacing them", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        timing.hold_to_read(self, cap2, settle=0.6)

        # 2. the reduced simulations: with and without Planet Nine
        foot = VGroup(*[Dot(ax.c2p(np.log10(min(f["a_au"], 9999.0)), f["q_au"]), radius=0.035,
                            color=P.TEAL).set_opacity(0.6) for f in wp["footprints"]])
        tally = VGroup(
            layout.label("with Planet Nine", font_size=18, color=P.BLUE, weight="BOLD"),
            layout.label(f"{100 * wp['crossing_fraction']:.0f}% of snapshots cross",
                         font_size=17, color=P.TEAL),
            layout.label("without Planet Nine", font_size=18, color=P.FG, weight="BOLD"),
            layout.label(f"{wo['n_crossing']} of {wo['n_footprints']:,} cross", font_size=17,
                         color=P.RED),
        ).arrange(DOWN, buff=0.12, aligned_edge=LEFT).move_to([4.6, 1.6, 0])
        tally[2].shift(DOWN * 0.2)
        tally[3].shift(DOWN * 0.2)
        note = layout.label(f"{run['n_particles']} particles, {run['t_myr']:.0f} Myr, "
                            f"Planet Nine {run['p9_mass_earth']:.0f} M⊕\nto speed up its pull",
                            font_size=14, color=P.MUTED).next_to(tally, DOWN, buff=0.35,
                                                                 aligned_edge=LEFT)
        cap3 = layout.caption("Simulate the distant disk: Planet Nine walks perihelia "
                              "inward, across Neptune's orbit", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3), run_time=0.6)
        self.play(FadeIn(foot, lag_ratio=0.01), FadeIn(tally[:2]), FadeIn(note), run_time=2.0)
        timing.hold_to_read(self, cap3, settle=0.3)
        cap4 = layout.caption("Run the same disk without it: nothing crosses",
                              font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        self.play(FadeIn(tally[2:]))
        timing.hold_to_read(self, cap4, settle=0.6)
        self.play(FadeOut(VGroup(ax, furn, nep, nep_l, dots, foot, tally, note, cap4)))

        # 3. the test: where does each object sit in the model's prediction?
        x0, x1, y = -5.0, 5.0, 1.5
        line = Line([x0, y, 0], [x1, y, 0], color=P.MUTED, stroke_width=2)
        ends = VGroup(
            layout.label("0: lower q than the model allows", font_size=15, color=P.MUTED)
            .next_to([x0, y, 0], DOWN, buff=0.2, aligned_edge=LEFT),
            layout.label("1", font_size=15, color=P.MUTED).next_to([x1, y, 0], DOWN, buff=0.2),
        )
        head = layout.label("each object's perihelion as a percentile of the model's "
                            "prediction (at its discovery distance)", font_size=17)
        head.move_to([0, 2.5, 0])
        ticks = VGroup(*[Line([x0 + (x1 - x0) * v, y - 0.25, 0], [x0 + (x1 - x0) * v, y + 0.25, 0],
                              color=P.TEAL, stroke_width=4) for v in wp["xi_sorted"]])
        ks = layout.label(f"with Planet Nine: spread evenly, KS p = {wp['ks_p']:.2f} "
                          f"(paper {d['paper']['ks_p_p9']:.2f})", font_size=17, color=P.TEAL)
        ks.next_to(line, DOWN, buff=0.7)
        ks_free = layout.label(f"without: all {len(objs)} pile at 0, KS p = {wo['ks_p']:.0e}",
                               font_size=17, color=P.RED).next_to(ks, DOWN, buff=0.15)
        cap5 = layout.caption("If the model is right, the percentiles fall evenly "
                              "between 0 and 1", font_size=22)
        self.play(FadeIn(head), Create(line), FadeIn(ends), FadeIn(cap5))
        self.play(Create(ticks, lag_ratio=0.1), run_time=1.5)
        self.play(FadeIn(ks))
        self.play(FadeIn(ks_free))
        timing.hold_to_read(self, cap5, ks, ks_free, settle=0.6)
        self.play(FadeOut(VGroup(line, ends, head, ticks, ks, ks_free, cap5)))

        # 4. the paper's statistic on the null distribution
        edges = np.array(d["null"]["edges"])
        frac = np.array(d["null"]["fraction"])
        keep = edges[:-1] >= -20.0
        frac = frac[keep]
        edges = edges[np.append(keep, True)]
        ax2, lab2 = edge_axes(
            [-20, 0, 5], [0, 0.12, 0.04], "ζ = Σ log₁₀(percentile)  (lower = worse fit)",
            "fraction of random draws", x_length=9.0, y_length=3.8)
        VGroup(ax2, lab2).move_to([0, 0.3, 0])
        bars = widgets.histogram(ax2, edges, frac, color=P.MUTED, opacity=0.6)
        pp = d["paper"]
        m_p9 = widgets.marker_line(ax2, pp["zeta_p9"], (0, 0.12),
                                   f"with P9: {pp['zeta_p9']:.1f}", color=P.BLUE, side=UP)
        m_free = widgets.marker_line(ax2, pp["zeta_free"], (0, 0.09),
                                     f"without: {pp['zeta_free']:.1f}", color=P.RED, side=UP)
        cap6 = layout.caption(f"{len(objs)} random percentiles, {d['null']['n_draws']:,} times: "
                              "the spread a correct model gives", font_size=22)
        self.play(Create(ax2), FadeIn(lab2), FadeIn(cap6))
        self.play(FadeIn(bars, lag_ratio=0.03), run_time=1.2)
        timing.hold_to_read(self, cap6, settle=0.3)
        sig = d["paper_zeta_free_sigma"]
        cap7 = layout.caption(f"Paper's full runs: with Planet Nine typical; without, "
                              f"{sig:.1f}σ out in the tail", font_size=22)
        self.play(FadeOut(cap6), FadeIn(cap7), run_time=0.6)
        self.play(Create(m_p9))
        self.play(Create(m_free))
        timing.hold_to_read(self, cap7, settle=0.8)
        self.play(FadeOut(cap7))

        layout.show_takeaway(
            self, f"Neptune crossers reject a Planet-Nine-free solar system at ~{sig:.0f}σ.")
