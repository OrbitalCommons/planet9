"""Benakli (2026) -- A Solar-System window for hidden stellar companions.

A tidal mass-distance envelope, calibrated to Planet-Nine-like constraints,
caps the mass an unseen companion can have at a given distance: M_max grows as
distance cubed. That window admits Earth-to-sub-Saturn masses at hundreds of AU,
rising to ~1 Jupiter mass near 2000 AU -- far more than any smooth dark-matter
halo could supply. Reproduced in p9-2026-stellar-companions.
"""
import numpy as np
from manim import (
    Create,
    DOWN,
    Dot,
    FadeIn,
    FadeOut,
    LEFT,
    Line,
    Polygon,
    RIGHT,
    Scene,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2026-stellar-companions"

# Tidal envelope, anchored to a Planet-Nine-like point: M_max(d) = M_* (d/d_*)^3
D_STAR = 300.0     # AU
M_STAR = 1.0       # Earth masses


def m_max(d):
    return M_STAR * (d / D_STAR) ** 3


class StellarCompanions2026(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "A window for hidden companions", cite="Benakli (2026)")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        # mass (Y, log) vs distance (X, log). Axes drawn over log10 coordinates.
        # X: 100..3000 AU -> log10 2..3.48 ; Y: 0.1..1000 Mearth -> log10 -1..3
        x_lo, x_hi = 2.0, np.log10(3000.0)
        y_lo, y_hi = -1.0, 3.0
        ax = widgets.axes([x_lo, x_hi, 1.0], [y_lo, y_hi, 1.0], x_length=9.6, y_length=4.4,
                          font_size=18, shift_down=0.45)
        xlab = layout.label("distance from Sun (AU)", font_size=18, color=P.FG).next_to(ax, DOWN, buff=0.55)
        ylab = layout.label("mass (Earth masses)", font_size=18, color=P.FG).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.55)
        self.play(Create(ax), FadeIn(xlab), FadeIn(ylab))

        # tick relabels in real units (decades)
        tick_grp = VGroup()
        for au in (100, 300, 1000, 3000):
            t = layout.label(str(au), font_size=14, color=P.MUTED)
            t.next_to(ax.c2p(np.log10(au), y_lo), DOWN, buff=0.14)
            tick_grp.add(t)
        for me, txt in [(-1, "0.1"), (0, "1"), (1, "10"), (2, "100"), (3, "1000")]:
            t = layout.label(txt, font_size=14, color=P.MUTED)
            t.next_to(ax.c2p(x_lo, me), LEFT, buff=0.14)
            tick_grp.add(t)
        self.add(tick_grp)

        # the M_max(d) ∝ d³ envelope line (in log-log it is a straight line)
        ds = np.linspace(100.0, 3000.0, 80)
        lx = np.log10(ds)
        ly = np.clip(np.log10(m_max(ds)), y_lo, y_hi)
        env_line = ax.plot_line_graph(lx, ly, line_color=P.TEAL, add_vertex_dots=False, stroke_width=3)

        # shade the ALLOWED window: below the envelope, above the floor
        pts = [ax.c2p(x, y) for x, y in zip(lx, ly)]
        pts += [ax.c2p(lx[-1], y_lo)]
        pts += [ax.c2p(lx[0], y_lo)]
        allowed = Polygon(*pts, color=P.TEAL, fill_opacity=0.16, stroke_width=0)
        self.play(FadeIn(allowed), Create(env_line))

        env_lbl = layout.label("M_max ∝ d³  (tidal limit)", font_size=16, color=P.TEAL, weight="BOLD")
        env_lbl.move_to(ax.c2p(np.log10(165), np.log10(400)))
        win_lbl = layout.label("allowed window", font_size=18, color=P.TEAL, weight="BOLD")
        win_lbl.move_to(ax.c2p(np.log10(1500), np.log10(0.18)))
        self.play(FadeIn(env_lbl), FadeIn(win_lbl))

        # anchor points along the envelope
        anchors = [
            (300, 1, "300 AU → ~1 M⊕"),
            (1000, 40, "1000 AU → ~40 M⊕ (sub-Saturn)"),
            (2000, 320, "2000 AU → ~320 M⊕ (≈1 Jupiter)"),
        ]
        ups = [DOWN, UP, UP]
        anchor_grp = VGroup()
        for (d, m, name), updir in zip(anchors, ups):
            dot = Dot(ax.c2p(np.log10(d), np.log10(m)), radius=0.07, color=P.GREEN)
            lbl = layout.label(name, font_size=14, color=P.GREEN)
            if updir is DOWN:
                lbl.next_to(dot, DOWN + RIGHT, buff=0.08)
            else:
                lbl.next_to(dot, updir, buff=0.1)
            anchor_grp.add(dot, lbl)
            self.play(FadeIn(dot), FadeIn(lbl), run_time=0.5)

        timing.hold_to_read(self, win_lbl, settle=0.8)

        eq = layout.explain_equation(
            self,
            [r"M_{\max}(d)", "=", r"M_\star", r"\left(\dfrac{d}{d_\star}\right)^{3}"],
            [
                (0, "the heaviest hidden companion still allowed at distance d"),
                (2, "anchored to a Planet-Nine-like reference mass and distance"),
                (3, "a tidal-constraint envelope: the cap grows as distance cubed"),
            ],
            scale=0.8,
            where=DOWN * 1.55,
        )
        self.play(FadeOut(eq), FadeOut(allowed), FadeOut(env_line), FadeOut(env_lbl),
                  FadeOut(win_lbl), FadeOut(anchor_grp), FadeOut(ax), FadeOut(xlab), FadeOut(ylab),
                  FadeOut(tick_grp))

        # second beat: a smooth dark-matter halo cannot supply such a thing
        readout = paper.result_readout(
            "dark-matter halo mass within 1000 AU", "< Pluto", color=P.PURPLE)
        readout.move_to(DOWN * 0.3)
        note = layout.label(
            "the local dark-matter density integrates to only a sub-Pluto mass —\n"
            "a smooth halo cannot be the source",
            font_size=20, color=P.FG)
        note.next_to(readout, DOWN, buff=0.45)
        self.play(FadeIn(readout, shift=UP * 0.15))
        self.play(FadeIn(note))
        timing.hold_to_read(self, note, settle=1.2)
        self.play(FadeOut(readout), FadeOut(note))

        layout.show_takeaway(
            self, "A hidden companion must be a bound, gravity-only object — not dark matter.")
