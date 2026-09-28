"""Batygin, Mardling & Nesvorny (2021) -- the stability boundary of the scattered disk.

A long, eccentric orbit feels Neptune only near perihelion. Each passage is a
kick, and the kicks organise into an infinite chain of 2:j resonances with
Neptune. When the perihelion is low the resonances are wide enough to overlap
and the orbit wanders chaotically; when it is high they separate and the orbit
is frozen. The overlap condition gives a closed-form critical perihelion.
Reproduced in p9-2021-stability: the resonance chains, the Sun + Neptune
N-body tracks and the boundary curve are the crate's own (anim.json -> papers
-> p9-2021-stability).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    MathTex,
    RoundedRectangle,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2021-stability"


class Stability2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        a_ref = d["a_ref_au"]
        self.add(paper.scene_header(CRATE))

        # 1. Neptune's 2:j resonances near a = 500 AU, at two perihelia
        lo, hi = 486.0, 516.0
        x0, x1 = -4.6, 5.6

        def xa(a):
            return x0 + (x1 - x0) * (a - lo) / (hi - lo)

        axis = VGroup(Line([x0, -1.9, 0], [x1, -1.9, 0], color=P.MUTED, stroke_width=2))
        for a in range(490, 516, 5):
            axis.add(Line([xa(a), -1.97, 0], [xa(a), -1.83, 0], color=P.MUTED, stroke_width=2))
            axis.add(layout.label(str(a), font_size=15, color=P.MUTED).move_to([xa(a), -2.2, 0]))
        axis.add(layout.label("semi-major axis a (AU)", font_size=17).move_to([0.5, -2.6, 0]))
        rows = VGroup()
        tags = VGroup()
        chains = sorted(d["chains"], key=lambda c: -c["q_au"])
        for k, ch in enumerate(chains):
            y = 1.3 - 2.0 * k
            chaotic = ch["q_au"] < d["q_crit_headline_au"]
            col = P.ORANGE if chaotic else P.GREEN
            row = VGroup()
            for r in ch["resonances"]:
                a_l = max(r["a_au"] - r["half_width_au"], lo)
                a_r = min(r["a_au"] + r["half_width_au"], hi)
                if a_r <= a_l:
                    continue
                box = RoundedRectangle(width=xa(a_r) - xa(a_l), height=0.55, corner_radius=0.12,
                                       color=col, stroke_width=1.5)
                box.set_fill(col, opacity=0.25).move_to([(xa(a_l) + xa(a_r)) / 2, y, 0])
                row.add(box)
            rows.add(row)
            tags.add(VGroup(
                layout.label(f"perihelion {ch['q_au']:.0f} AU", font_size=18, color=col),
                layout.label("islands overlap: chaos" if chaotic else "islands apart: stable",
                             font_size=15, color=col),
            ).arrange(DOWN, buff=0.08, aligned_edge=LEFT).next_to([x0, y, 0], LEFT, buff=0.2))
        self.play(FadeIn(axis))
        cap = layout.caption(f"Neptune kicks a distant orbit at each perihelion: its 2:j "
                             f"resonances sit {d['spacing_au']:.1f} AU apart", font_size=22)
        self.play(FadeIn(cap), Create(rows[0], lag_ratio=0.1), FadeIn(tags[0]), run_time=1.5)
        timing.hold_to_read(self, cap, settle=0.4)
        cap2 = layout.caption("Bring the perihelion closer to Neptune and every island widens",
                              font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), run_time=0.6)
        self.play(Create(rows[1], lag_ratio=0.1), FadeIn(tags[1]), run_time=1.5)
        timing.hold_to_read(self, cap2, settle=0.8)
        self.play(FadeOut(VGroup(axis, rows, tags, cap2)))

        # 2. the same two perihelia under Sun + Neptune, integrated
        nb = sorted(d["nbody"], key=lambda n: -n["q_au"])
        n_pts = len(nb[0]["a_au"][0])
        t_myr = np.arange(n_pts) * nb[0]["dt_myr"]
        prog = ValueTracker(0.0)
        panels = VGroup()
        live = []
        for k, run in enumerate(nb):
            chaotic = run["q_au"] < d["q_crit_headline_au"]
            col = P.ORANGE if chaotic else P.GREEN
            ax = widgets.axes([0, 1.0, 0.25], [380, 680, 100], x_length=5.0, y_length=3.7,
                              shift_down=0)
            ax.move_to([-3.2 + 6.6 * k, 0.35, 0])
            ax.add_coordinates()
            head = layout.label(f"perihelion {run['q_au']:.0f} AU", font_size=19, color=col,
                                weight="BOLD").next_to(ax, UP, buff=0.15)
            xl = layout.label("time (Myr)", font_size=15).next_to(ax, DOWN, buff=0.1)
            yl = layout.label("a (AU)", font_size=15).rotate(np.pi / 2).next_to(ax, LEFT,
                                                                                buff=0.1)
            panels.add(VGroup(ax, head, xl, yl))

            def tracks(ax=ax, run=run, col=col):
                m = max(2, int(prog.get_value() * (n_pts - 1)) + 1)
                return VGroup(*[widgets.curve(ax, t_myr[:m], np.clip(tr[:m], 380, 680),
                                              color=col, stroke_width=1.8)
                                for tr in run["a_au"]])
            live.append(always_redraw(tracks))
        cap3 = layout.caption("Integrate 8 orbits for 1 Myr at a = 500 AU "
                              "with the Sun and Neptune only", font_size=22)
        self.play(FadeIn(panels), FadeIn(cap3))
        self.add(*live)
        self.play(prog.animate.set_value(1.0), run_time=7.0, rate_func=linear)
        for m in live:
            m.clear_updaters()
        ly = d["lyapunov_time_yr"]
        cap4 = layout.caption(f"Below the boundary a random-walks by tens of AU; chaos sets in "
                              f"within ~{ly / 1000:.0f},000 yr, one orbit", font_size=22)
        self.play(FadeOut(cap3), FadeIn(cap4), run_time=0.6)
        timing.hold_to_read(self, cap4, settle=0.8)
        self.play(FadeOut(VGroup(panels, *live, cap4)))

        # 3. the boundary, in closed form
        eq = MathTex(r"q_{\rm crit} = a_{\rm N}\,\sqrt{\ln\!\Big(\tfrac{24^2}{5}\,"
                     r"\tfrac{m_{\rm N}}{M_\odot}\,\big(\tfrac{a}{a_{\rm N}}\big)^{5/2}\Big)}",
                     color=P.FG).scale(0.8).move_to([0, 2.45, 0])
        b = d["boundary"]
        ax3, lab3 = widgets.labeled_axes(
            [200, 1600, 200], [0, 90, 30], x_label="semi-major axis a (AU)",
            y_label="perihelion q (AU)", y_rotate=True, numbers=True,
            x_length=8.4, y_length=3.2, shift_down=0)
        VGroup(ax3, lab3).move_to([-1.0, -0.3, 0])
        crit = widgets.curve(ax3, b["a_au"], b["q_crit_au"], color=P.ORANGE, stroke_width=3.5)
        above = layout.label("detached: frozen", font_size=17, color=P.GREEN)
        above.move_to(ax3.c2p(1250, 80))
        below = layout.label("scattering: chaotic", font_size=17, color=P.ORANGE)
        below.move_to(ax3.c2p(1250, 25))
        self.play(FadeIn(eq))
        cap5 = layout.caption("Overlap sets in below a critical perihelion "
                              "that rises slowly with a", font_size=22)
        self.play(FadeIn(cap5), Create(ax3), FadeIn(lab3))
        self.play(Create(crit), FadeIn(above), FadeIn(below), run_time=1.5)
        timing.hold_to_read(self, cap5, eq, settle=0.4)
        dots = VGroup()
        for o in d["objects"]:
            col = P.ORANGE if o["chaotic"] else P.GREEN
            dots.add(Dot(ax3.c2p(o["a_au"], o["q_au"]), radius=0.07, color=col).set_z_index(4))
        qc = d["q_crit_headline_au"]
        mark = Dot(ax3.c2p(a_ref, qc), radius=0.09, color=P.ORANGE).set_z_index(5)
        mark_l = layout.label(f"{qc:.1f} AU at a = {a_ref:.0f} AU", font_size=17,
                              color=P.ORANGE).next_to(mark, DOWN + RIGHT, buff=0.1)
        cap6 = layout.caption(f"The {d['n_objects']} clustered objects: {d['n_chaotic']} sits "
                              "in the chaotic zone, the rest are frozen in place", font_size=22)
        self.play(FadeOut(cap5), FadeIn(cap6), run_time=0.6)
        self.play(FadeIn(dots, lag_ratio=0.1), FadeIn(mark), FadeIn(mark_l))
        timing.hold_to_read(self, cap6, settle=0.8)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, f"Overlapping 2:j resonances end Neptune's chaos at q ≈ {qc:.0f} AU "
                  f"(a = {a_ref:.0f} AU).")
