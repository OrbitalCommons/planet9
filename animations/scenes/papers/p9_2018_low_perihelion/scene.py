"""Cáceres & Gomes (2018) -- the case for a low-perihelion Planet Nine.

The favoured Planet Nine kept its perihelion q9 = a9 (1 - e9) at 200-400 AU,
far outside the clustered TNOs. Cáceres & Gomes lowered it to 60-100 AU and
found tighter angular clustering. The reproduction keeps a9 fixed and raises
e9, so only q9 changes, and measures how widely a test TNO's perihelion
direction can swing inside Planet Nine's secular well (the half-width of the
libration island). The width is a spread, not a direction: this single-planet
secular model does not by itself fix the observed anti-alignment. Reproduced
in p9-2018-low-perihelion; the sweeps, the spreads and the real TNO orbits are
the crate's own (anim.json -> papers -> p9-2018-low-perihelion).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    AnnularSector,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2018-low-perihelion"

SUN_POS = np.array([-4.55, -0.35, 0.0])
AU_SCALE = 0.0026          # scene units per AU in the top view
PLOT_CENTRE = np.array([2.75, 0.2, 0.0])
DOWNWARD = -np.pi / 2      # P9's perihelion points down the column


def track(rows):
    """(q9, half-width) arrays sorted by increasing q9."""
    q = np.array([r["q9_au"] for r in rows])
    w = np.array([r["half_width_deg"] for r in rows])
    order = np.argsort(q)
    return q[order], w[order]


class LowPerihelion2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        a9 = d["a9_au"]
        e_lo, e_hi = d["e9_range"]
        q700, w700 = track(d["track_700"])
        q1500, w1500 = track(d["track_1500"])
        self.add(paper.scene_header(CRATE))

        # 1. the real clustered TNOs and the canonical Planet Nine, seen from above
        turn = DOWNWARD - np.radians(d["p9_varpi_deg"])
        s = AU_SCALE
        etnos = VGroup(*[
            orbits.ellipse_orbit(o["a"] * s, o["e"], color=P.GREEN,
                                 varpi=np.radians(o["varpi_deg"]) + turn,
                                 stroke_width=1.6, opacity=0.85).shift(SUN_POS)
            for o in d["etnos"]])
        nep = Circle(radius=30 * s, color=P.FG, stroke_width=1.2).move_to(SUN_POS)
        sun = orbits.sun(radius=0.045).move_to(SUN_POS)
        etno_lab = layout.label(f"{len(d['etnos'])} clustered TNOs\n(real orbits)",
                                font_size=14, color=P.GREEN)
        etno_lab.to_edge(LEFT, buff=0.3).set_y(-2.6)
        cap0 = layout.caption("Seen from above: the real clustered TNOs. Neptune's whole orbit "
                              "is the tiny ring at the Sun", font_size=20)
        self.play(FadeIn(sun), Create(nep), FadeIn(cap0), run_time=0.8)
        self.play(Create(etnos, lag_ratio=0.1), FadeIn(etno_lab), run_time=1.6)
        timing.hold_to_read(self, cap0, settle=0.3)

        e9 = ValueTracker(e_lo)

        def q9():
            return a9 * (1.0 - e9.get_value())

        p9_orbit = always_redraw(lambda: orbits.ellipse_orbit(
            a9 * s, e9.get_value(), color=P.BLUE, stroke_width=3,
            varpi=DOWNWARD).shift(SUN_POS))

        def peri_pos():
            return SUN_POS + np.array([0, -q9() * s, 0])

        peri = always_redraw(lambda: Dot(peri_pos(), radius=0.07, color=P.BLUE))

        def q_tag():
            lab = layout.label(f"q9 = {q9():.0f} AU", font_size=16, color=P.BLUE,
                               weight="BOLD")
            lab.move_to(peri_pos() + np.array([2.05, 0, 0]))
            lead = DashedLine(peri_pos() + np.array([0.1, 0, 0]),
                              lab.get_left() + np.array([-0.08, 0, 0]),
                              color=P.BLUE, stroke_width=1.5)
            return VGroup(lead, lab)

        q_lab = always_redraw(q_tag)
        p9_lab = layout.label(f"Planet Nine\na9 = {a9:.0f} AU", font_size=14, color=P.BLUE)
        p9_lab.move_to([-6.35, 2.75, 0])
        cap = layout.caption(f"The favoured Planet Nine stays far out: perihelion "
                             f"{q9():.0f} AU, outside every clustered TNO", font_size=20)
        self.play(Create(p9_orbit), FadeIn(peri), FadeIn(q_lab), FadeIn(p9_lab), FadeOut(cap0),
                  FadeIn(cap),
                  run_time=1.4)
        timing.hold_to_read(self, cap, settle=0.6)

        # 2. the measure: how far a TNO's perihelion direction can swing
        ax, labels = widgets.labeled_axes(
            [50, 300, 50], [0, 110, 20], x_label="Planet Nine perihelion q9 (AU)",
            y_label="TNO apse swing, ± degrees", y_rotate=True, numbers=True,
            x_length=6.6, y_length=4.4, shift_down=0, font_size=18)
        VGroup(ax, labels).move_to(PLOT_CENTRE)
        curve700 = widgets.curve(ax, q700, w700, color=P.BLUE)
        lab700 = layout.label(f"a9 = {a9:.0f} AU", font_size=14, color=P.BLUE)
        lab700.next_to(ax.c2p(q700[-1], w700[-1]), DOWN, buff=0.15)

        def fan():
            w = float(np.interp(q9(), q700, w700))
            apex = ax.c2p(62, 70)
            wedge = AnnularSector(inner_radius=0, outer_radius=0.85, angle=np.radians(2 * w),
                                  start_angle=np.radians(90 - w), color=P.GREEN,
                                  fill_opacity=0.25, stroke_width=0).move_arc_center_to(apex)
            spine = Line(apex, apex + 0.85 * UP, color=P.GREEN, stroke_width=2)
            txt = layout.label(f"±{w:.0f}°", font_size=16, color=P.GREEN, weight="BOLD")
            txt.next_to(apex, DOWN, buff=0.08)
            return VGroup(wedge, spine, txt)

        fan_m = always_redraw(fan)
        mark = always_redraw(lambda: Dot(ax.c2p(q9(), float(np.interp(q9(), q700, w700))),
                                         radius=0.08, color=P.BLUE))
        cap2 = layout.caption("Its pull traps each TNO's perihelion inside a wedge: "
                              "the narrower, the tighter the cluster", font_size=20)
        self.play(Create(ax), FadeIn(labels), FadeIn(mark), FadeIn(fan_m), FadeOut(cap),
                  FadeIn(cap2), run_time=1.2)
        timing.hold_to_read(self, cap2, settle=0.6)

        # 3. lower the perihelion: the wedge closes
        cap3 = layout.caption(f"Same orbit size, more eccentric: q9 falls from {q700[-1]:.0f} "
                              f"to {q700[0]:.0f} AU and the wedge closes", font_size=20)
        self.play(FadeOut(cap2), FadeIn(cap3))
        self.play(Create(curve700), e9.animate.set_value(e_hi), run_time=6.0)
        self.play(FadeIn(lab700))
        timing.hold_to_read(self, cap3, settle=0.8)

        # 4. the paper's wider planet, and where the Kuiper belt says stop
        curve1500 = widgets.curve(ax, q1500, w1500, color=P.TEAL)
        lab1500 = layout.label(f"a9 = {d['wide_a9_au']:.0f} AU", font_size=14, color=P.TEAL)
        lab1500.next_to(ax.c2p(q1500[-1], w1500[-1]), DOWN, buff=0.15).shift(LEFT * 0.3)
        wreck = Rectangle(width=ax.c2p(70, 0)[0] - ax.c2p(50, 0)[0],
                          height=ax.c2p(0, 110)[1] - ax.c2p(0, 0)[1],
                          stroke_width=0).set_fill(P.RED, opacity=0.18)
        wreck.move_to(ax.c2p(60, 55))
        wreck_lab = layout.label("paper: q9 = 60 AU\nempties the\nKuiper belt", font_size=13,
                                 color=P.RED)
        wreck_lab.next_to(ax.c2p(70, 104), RIGHT, buff=0.1)
        w90 = float(np.interp(90.0, q1500, w1500))
        pick = Dot(ax.c2p(90, w90), radius=0.08, color=P.TEAL)
        pick_line = DashedLine(ax.c2p(90, 0), ax.c2p(90, w90), color=P.TEAL, stroke_width=2)
        pick_lab = layout.label("paper's pick\nq9 = 90 AU", font_size=13, color=P.TEAL)
        pick_lab.next_to(ax.c2p(90, w90), DOWN + RIGHT, buff=0.12)
        cap4 = layout.caption("For their wider 1500 AU planet too, lower q9 means tighter "
                              "clustering, down to where the Kuiper belt breaks", font_size=19)
        self.play(FadeOut(fan_m), Create(curve1500), FadeIn(lab1500), FadeOut(cap3),
                  FadeIn(cap4), run_time=1.6)
        self.play(FadeIn(wreck), FadeIn(wreck_lab), Create(pick_line), FadeIn(pick),
                  FadeIn(pick_lab), run_time=1.2)
        timing.hold_to_read(self, cap4, settle=1.0)

        tally = paper.result_readout(
            f"apse swing, q9 {q700[-1]:.0f} → {q700[0]:.0f} AU",
            f"±{d['half_width_canonical_deg']:.0f}° → ±{d['half_width_low_deg']:.0f}°",
            color=P.GREEN).scale(0.75)
        tally.move_to(ax.c2p(235, 22))
        self.play(FadeOut(cap4), FadeIn(tally))
        timing.hold_to_read(self, tally, settle=0.8)

        layout.show_takeaway(
            self, "A closer-swinging Planet Nine clusters the TNOs more tightly.")
