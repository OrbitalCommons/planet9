"""Benakli (2026) -- a Solar-System window for hidden stellar companions.

How massive could a dark companion of the Sun be, at a given distance, without
the planetary ephemerides having noticed? The tidal pull on the planets scales
as M/d^3, so the allowed mass grows as the cube of distance. The crate computes
that envelope (anchored at 5 Earth masses, 500 AU) and the mass of smooth dark
matter enclosed within the same radius, which is tens of thousands of times too
small. Everything is from anim.json -> papers -> p9-2026-stellar-companions.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    Dot,
    DoubleArrow,
    FadeIn,
    FadeOut,
    Polygon,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2026-stellar-companions"

X_LO, X_HI = 2.0, 3.6          # log10 distance (AU)
Y_LO, Y_HI = -6.0, 4.0         # log10 mass (Earth masses)


def rows(items, font_size=15, buff=0.14):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=LEFT)


def inside(xs, ys):
    """Points of a curve that fall inside the plot, as log10 pairs."""
    out = []
    for x, y in zip(np.log10(xs), np.log10(ys)):
        if X_LO <= x <= X_HI and Y_LO <= y <= Y_HI:
            out.append((x, y))
    return out


class StellarCompanions2026(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        p9 = d["planet_nine"]

        self.add(paper.scene_header(CRATE))

        # Axes cross at their origin, so plot offsets from the lower-left corner
        # (log10 distance - X_LO, log10 mass - Y_LO) and label ticks by hand.
        ax = widgets.axes([0, X_HI - X_LO, 10], [0, Y_HI - Y_LO, 10], x_length=7.4,
                          y_length=4.2, shift_down=-0.55)

        def q(x, y):
            return ax.c2p(x - X_LO, y - Y_LO)

        marks = VGroup()
        for au in (100, 300, 1000, 3000):
            marks.add(layout.label(f"{au:,}", font_size=14, color=P.MUTED)
                      .next_to(q(np.log10(au), Y_LO), DOWN, buff=0.12))
        xl = layout.label("distance from the Sun (AU)", font_size=17, color=P.FG)
        xl.next_to(ax, DOWN, buff=0.45)
        yl = layout.label("mass  (logarithmic)", font_size=15, color=P.FG).rotate(np.pi / 2)
        yl.next_to(ax, LEFT, buff=1.05)
        bodies = [("Pluto", d["pluto_earth"]), ("Earth", 1.0), ("Saturn", d["saturn_earth"]),
                  ("Jupiter", d["jupiter_earth"])]
        guides = VGroup()
        for name, m in bodies:
            y = np.log10(m)
            guides.add(DashedLine(q(X_LO, y), q(X_HI, y), color=P.MUTED,
                                  stroke_width=1.0).set_stroke(opacity=0.5))
            marks.add(layout.label(name, font_size=14, color=P.MUTED)
                      .next_to(q(X_LO, y), LEFT, buff=0.12))
        frame = VGroup(ax, marks, xl, yl, guides)
        frame.shift(LEFT * 1.9)

        # 1. the envelope the ephemerides allow, swept out along distance
        env = inside(d["distance_au"], d["envelope_earth"])
        ex = np.array([p[0] for p in env])
        ey = np.array([p[1] for p in env])
        line = widgets.curve(ax, ex - X_LO, ey - Y_LO, color=P.TEAL)
        allowed = Polygon(*[q(x, y) for x, y in env], q(env[-1][0], Y_LO),
                          q(env[0][0], Y_LO), stroke_width=0).set_fill(P.TEAL, opacity=0.1)
        allowed_lab = layout.label("allowed: too weak a tug to notice", font_size=14,
                                   color=P.TEAL).move_to(q(3.05, -1.2))
        cap = layout.caption(
            "A distant mass tugs the planets as M/d³: ten times farther can hide "
            "a thousand times more", font_size=21)
        self.play(FadeIn(frame), run_time=1.0)
        self.play(Create(line), FadeIn(allowed), FadeIn(cap), run_time=1.5)
        self.play(FadeIn(allowed_lab))

        logd = ValueTracker(np.log10(300.0))

        def reading():
            x = logd.get_value()
            y = float(np.interp(x, ex, ey))
            dot = Dot(q(x, y), radius=0.08, color=P.TEAL)
            txt = layout.label(f"{10 ** x:,.0f} AU: up to {10 ** y:,.1f} M⊕" if y < 1.5
                               else f"{10 ** x:,.0f} AU: up to {10 ** y:,.0f} M⊕",
                               font_size=15, color=P.TEAL, weight="BOLD")
            txt.next_to(dot, UP + LEFT, buff=0.08)
            return VGroup(dot, txt)

        probe = always_redraw(reading)
        self.add(probe)
        self.play(logd.animate.set_value(np.log10(d["distance_for_jupiter"])), run_time=3.5)
        self.wait(0.6)
        self.remove(probe)
        timing.hold_to_read(self, cap, settle=0.2)

        nine = Dot(q(np.log10(p9["a_au"]), np.log10(p9["mass_earth"])), radius=0.08,
                   color=P.BLUE)
        nine_lab = layout.label("Planet Nine", font_size=14, color=P.BLUE)
        nine_lab.next_to(nine, UP + LEFT, buff=0.06)
        paper_pts = VGroup(*[
            Dot(q(np.log10(r["distance_au"]), np.log10(r["published_mass_earth"])),
                radius=0.045, color=P.FG)
            for r in d["table"]])
        key = rows([
            ("heaviest companion the", P.TEAL),
            ("planets' motions allow", P.TEAL),
            (f"{d['envelope_300']:.1f} M⊕ at 300 AU", P.TEAL),
            (f"{d['envelope_1000']:.0f} M⊕ at 1,000 AU", P.TEAL),
            (f"Saturn at {d['distance_for_saturn']:,.0f} AU", P.TEAL),
            (f"Jupiter at {d['distance_for_jupiter']:,.0f} AU", P.TEAL),
            ("white dots: the paper's table", P.FG),
            (f"{pub['envelope_1000']:.0f} M⊕ at 1,000 AU", P.FG),
        ], font_size=15)
        key[2:].shift(DOWN * 0.15)
        key[6:].shift(DOWN * 0.2)
        key.move_to([5.0, 1.45, 0])
        cap1b = layout.caption(
            f"Anchored where Planet Nine sits: {d['anchor']['mass_earth']:.0f} M⊕ at "
            f"{d['anchor']['distance_au']:.0f} AU is just allowed", font_size=21)
        self.play(FadeIn(nine), FadeIn(nine_lab), FadeIn(paper_pts), FadeIn(key),
                  FadeOut(cap), FadeIn(cap1b), run_time=1.2)
        timing.hold_to_read(self, cap1b, key, settle=0.8)

        # 2. smooth dark matter cannot fill it
        halo = inside(d["distance_au"], d["halo_earth"])
        halo_line = widgets.curve(ax, [p[0] - X_LO for p in halo],
                                  [p[1] - Y_LO for p in halo], color=P.RED)
        x_k = 3.0
        gap = DoubleArrow(q(x_k, np.log10(d["halo_1000_earth"])),
                          q(x_k, np.log10(d["envelope_1000"])), buff=0.08,
                          color=P.FG, stroke_width=2.5, tip_length=0.18)
        gap_lab = layout.label(f"× {d['halo_shortfall_1000']:,.0f}", font_size=16,
                               color=P.FG, weight="BOLD")
        gap_lab.next_to(gap, RIGHT, buff=0.1)
        key2 = rows([
            ("smooth dark matter", P.RED),
            ("inside the same radius:", P.RED),
            (f"{d['halo_1000_pluto']:.2f} Pluto masses", P.RED),
            ("within 1,000 AU", P.RED),
        ], font_size=15)
        key2.next_to(key, DOWN, buff=0.45, aligned_edge=LEFT)
        cap2 = layout.caption(
            f"The galaxy's dark matter, {d['rho_dm_msun_pc3']:.2f} solar masses per cubic "
            "parsec, adds up to almost nothing here", font_size=21)
        self.play(Create(halo_line), FadeIn(key2), FadeOut(allowed_lab), FadeOut(cap1b),
                  FadeIn(cap2), run_time=1.4)
        self.play(Create(gap), FadeIn(gap_lab))
        timing.hold_to_read(self, cap2, key2, settle=1.0)
        self.play(FadeOut(cap2))

        layout.show_takeaway(
            self, "There is room for a dark companion, but only a compact, bound one.")
