"""Russell & White (2025) -- what Planet Nine would look like.

Given a 6.6 M-Earth body, the most-likely composition is a mini-Neptune: a
2.0-2.6 R-Earth ball under an H-He envelope, geometric albedo 0.33-0.47. That
fixes an absolute magnitude H = -6.1 to -5.2 and an apparent magnitude of about
+21.9 to +22.7 near 625 AU. Reproduced in p9-2025-russell-albedo.
"""
from manim import (
    Annulus,
    Circle,
    Create,
    DOWN,
    FadeIn,
    FadeOut,
    LEFT,
    RIGHT,
    Scene,
    UP,
    VGroup,
    Write,
)

import p9_manim as P
from p9_manim import layout, paper, timing, widgets

CRATE = "p9-2025-russell-albedo"


class RussellAlbedo2025(Scene):
    def construct(self):
        tb = paper.title_block(CRATE, "What would Planet Nine look like?",
                               cite="Russell & White (2025)")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        bg = widgets.star_field(n=70, seed=2507, x=(-6, 6), y=(-2.2, 1.8),
                                radius=0.02, opacity=0.45)
        self.play(FadeIn(bg, lag_ratio=0.005))

        # The mini-Neptune disk with a faint reflected-light glow.
        glow = Annulus(inner_radius=0.62, outer_radius=1.1, color=P.BLUE,
                       fill_opacity=0.10, stroke_width=0).move_to(LEFT * 3.4)
        disk = Circle(radius=0.6, color=P.BLUE, fill_opacity=0.9,
                      stroke_color=P.TEAL, stroke_width=2).move_to(LEFT * 3.4)
        dlbl = layout.label("mini-Neptune\nH-He envelope 0.6-3.5%",
                            font_size=16, color=P.TEAL).next_to(disk, DOWN, buff=0.3)
        self.play(FadeIn(glow), Create(disk), FadeIn(dlbl))
        timing.hold_to_read(self, dlbl, settle=0.5)

        # mass -> radius -> albedo chain
        chain = VGroup(
            layout.label("mass  6.6 M_E", font_size=20, color=P.FG),
            layout.label("radius  2.0-2.6 R_E", font_size=20, color=P.BLUE),
            layout.label("albedo (V)  0.33-0.47", font_size=20, color=P.GREEN),
        ).arrange(DOWN, buff=0.45, aligned_edge=LEFT).move_to(RIGHT * 2.6 + UP * 0.4)
        self.play(FadeIn(chain, lag_ratio=0.3))
        timing.hold_to_read(self, chain, settle=0.8)
        self.play(FadeOut(chain), FadeOut(glow), FadeOut(disk), FadeOut(dlbl))

        # Reflected-light / absolute-magnitude relation, walked term by term.
        eq = layout.explain_equation(
            self,
            [r"m", "=", r"H", "+", r"5\log_{10}\!\big(r\,\Delta\big)"],
            [(0, "apparent magnitude -- how faint it looks from Earth"),
             (2, "absolute magnitude: brightness from reflected size and albedo"),
             (4, "distance dimming -- heliocentric r times Earth-distance Delta")],
            color=P.TEAL, scale=0.95, where=UP * 0.6)
        self.play(FadeOut(eq))

        ro = paper.result_readout("apparent magnitude near 625 AU",
                                  "m ~ 21.9-22.7", color=P.ORANGE).scale(0.72)
        ro.move_to(UP * 0.7)
        hbeat = layout.label("H = -6.1 to -5.2   (bright: 2.6 R_E / 0.47   faint: 2.0 R_E / 0.33)",
                             font_size=17, color=P.PURPLE).next_to(ro, DOWN, buff=0.4)
        self.play(FadeIn(ro), FadeIn(hbeat))
        timing.hold_to_read(self, ro, hbeat, settle=1.0)

        layout.show_takeaway(
            self, "A faint reflected dot near +22 mag -- within reach of world-class survey telescopes.")
