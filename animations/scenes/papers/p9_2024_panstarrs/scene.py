"""Brown, Holman & Batygin (2024) -- Pan-STARRS1 search & the combined exclusion.

Adding a PS1 3pi catalogue-linking search to ZTF and DES, the union of searches
now rules out ~78.5% of the Planet Nine prior. Reproduced in
p9-2024-panstarrs (`combined_exclusion`, which pins 0.564 / +0.050 / +0.171).
"""
import numpy as np
from manim import (
    Create,
    DOWN,
    FadeIn,
    LEFT,
    Rectangle,
    Scene,
    Text,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2024-panstarrs"


def _exclusion():
    """Real exclusion fractions from dataio.section('exclusion'); {} on failure."""
    try:
        return dataio.section("exclusion") or {}
    except Exception:
        return {}


class PanStarrs2024(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Stacking the surveys")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        W = 9.0
        track = Rectangle(width=W, height=1.0, color=P.MUTED, stroke_width=2)
        self.play(Create(track))

        ex = _exclusion()
        segs = [
            ("ZTF", ex.get("ztf", 0.564), P.RED),
            ("DES", ex.get("des_unique", 0.050), P.ORANGE),
            ("PS1", ex.get("ps1_unique", 0.171), P.PURPLE),
        ]
        x = -W / 2
        cum = 0.0
        for name, frac, col in segs:
            r = Rectangle(width=W * frac, height=1.0, color=col, stroke_width=0).set_fill(col, opacity=0.55)
            r.move_to([x + W * frac / 2, 0, 0])
            lbl = Text(f"{name}\n+{frac*100:.1f}%", font_size=16, color=P.FG).move_to(r.get_center())
            self.play(Create(r), FadeIn(lbl), run_time=0.8)
            x += W * frac
            cum += frac

        # remaining viable slice (real value if available, else 1 - cumulative)
        remaining = ex.get("remaining", 1 - cum)
        rem = Rectangle(width=W * remaining, height=1.0, color=P.TEAL, stroke_width=0).set_fill(P.TEAL, opacity=0.25)
        rem.move_to([x + W * remaining / 2, 0, 0])
        rem_lbl = Text(f"viable\n{remaining*100:.1f}%", font_size=16, color=P.TEAL).move_to(rem.get_center())
        self.play(Create(rem), FadeIn(rem_lbl))

        combined = ex.get("combined", cum)
        head = Text(f"combined: {combined*100:.1f}% excluded", font_size=24, color=P.RED).next_to(track, UP, buff=0.5)
        self.play(Write(head))
        timing.hold_to_read(self, head, settle=0.6)

        eq = layout.equation(
            r"f_{\rm comb} = 1 - \prod_s (1 - P_s)", color=P.RED, scale=0.85)
        eq.next_to(head, UP, buff=0.25)
        layout.show_equation(self, eq)

        layout.show_takeaway(
            self, "About one-fifth of the parameter space survives -- and that is where Rubin looks.")
