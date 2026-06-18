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
from p9_manim import layout, paper

CRATE = "p9-2024-panstarrs"


class PanStarrs2024(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "Stacking the surveys")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        W = 9.0
        track = Rectangle(width=W, height=1.0, color=P.MUTED, stroke_width=2)
        self.play(Create(track))

        segs = [("ZTF", 0.564, P.RED), ("DES", 0.050, P.ORANGE), ("PS1", 0.171, P.PURPLE)]
        x = -W / 2
        cum = 0.0
        for name, frac, col in segs:
            r = Rectangle(width=W * frac, height=1.0, color=col, stroke_width=0).set_fill(col, opacity=0.55)
            r.move_to([x + W * frac / 2, 0, 0])
            lbl = Text(f"{name}\n+{frac*100:.1f}%", font_size=16, color=P.FG).move_to(r.get_center())
            self.play(Create(r), FadeIn(lbl), run_time=0.8)
            x += W * frac
            cum += frac

        # remaining viable slice
        rem = Rectangle(width=W * (1 - cum), height=1.0, color=P.TEAL, stroke_width=0).set_fill(P.TEAL, opacity=0.25)
        rem.move_to([x + W * (1 - cum) / 2, 0, 0])
        rem_lbl = Text(f"viable\n{(1-cum)*100:.1f}%", font_size=16, color=P.TEAL).move_to(rem.get_center())
        self.play(Create(rem), FadeIn(rem_lbl))

        head = Text("combined: 78.5% excluded", font_size=24, color=P.RED).next_to(track, UP, buff=0.5)
        self.play(Write(head))
        self.play(FadeIn(layout.takeaway(
            "About one-fifth of the parameter space survives -- and that is where Rubin looks.")))
        self.wait(0.9)
