"""Hadden, Li, Payne & Holman (2017) -- chaotic dynamics from Planet Nine.

Mean-motion resonances with Planet Nine crowd together as an orbit approaches
the planet, and each grows wider as the orbit grows more eccentric. Where
neighbours overlap (Chirikov's criterion) the motion turns chaotic: the object
wanders from resonance to resonance instead of following a smooth secular path.
The scene sweeps the eccentricity to grow the overlapped zone, maps the chaotic
and regular parts of the (a, e) plane for a 500 AU Planet Nine at three masses,
and places the known distant objects on the map. Every number comes from
p9-2018-chaotic-dynamics (anim.json -> papers -> p9-2018-chaotic-dynamics).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Circle,
    Create,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    NumberLine,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2018-chaotic-dynamics"


def _boundary(ax, a_line, e_crit, a9, circular_zone, e_top=0.95):
    """The K = 1 curve: e_crit(a) where it exists, dropping to e = 0 at the
    inner edge of the circular overlap zone."""
    pts = [(a, e) for a, e in zip(a_line, e_crit) if e is not None and e <= e_top]
    pts.append((a9 - circular_zone, 0.0))
    vm = VMobject()
    vm.set_points_as_corners([ax.c2p(a, e) for a, e in pts])
    return vm


class ChaoticDynamics2018(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        a9, m9 = d["a9"], d["m9_earth"]
        e_grid = np.asarray(d["e_grid"])
        zone = np.asarray(d["overlap_zone_au"])

        self.add(paper.scene_header(CRATE))

        # ---- 1. resonances crowd toward the planet and merge ----------------
        line = NumberLine(x_range=[150, 520, 50], length=11.6, color=P.MUTED,
                          include_numbers=True, font_size=20, include_tip=False)
        line.move_to([0, -0.9, 0])
        xl = layout.label("semi-major axis a (AU)", font_size=18, color=P.FG)
        xl.next_to(line, DOWN, buff=0.45)
        p9 = Dot(line.n2p(a9), radius=0.13, color=P.BLUE).set_z_index(4)
        p9_lab = layout.label(f"Planet Nine, {a9:.0f} AU", font_size=16, color=P.BLUE)
        p9_lab.next_to(p9, DOWN, buff=0.5)
        self.play(Create(line), FadeIn(xl), FadeIn(p9), FadeIn(p9_lab), run_time=1.0)

        ticks = VGroup()
        names = VGroup()
        for r in d["first_order"]:
            if r["a"] > a9 - 4:
                continue
            x = line.n2p(r["a"])
            ticks.add(Line(x, x + UP * 1.5, color=P.ORANGE, stroke_width=2.0))
            if r["n"] <= 3:
                names.add(layout.label(f"{r['n'] + 1}:{r['n']}", font_size=15, color=P.ORANGE)
                          .next_to(x + UP * 1.5, UP, buff=0.08))
        for r in d["j_one"]:
            if r["j"] < 3:
                continue
            x = line.n2p(r["a"])
            ticks.add(Line(x, x + UP * 1.0, color=P.ORANGE, stroke_width=2.0, stroke_opacity=0.7))
            names.add(layout.label(f"{r['j']}:1", font_size=15, color=P.ORANGE)
                      .next_to(x + UP * 1.0, UP, buff=0.08))
        cap = layout.caption("Each tick is a resonance with Planet Nine: they crowd together "
                             "near the planet", font_size=21)
        self.play(Create(ticks, lag_ratio=0.03), FadeIn(names), FadeIn(cap), run_time=1.8)
        timing.hold_to_read(self, cap, settle=0.5)

        ecc = ValueTracker(0.0)

        def width_now():
            return float(np.interp(ecc.get_value(), e_grid, zone))

        def band():
            x0, x1 = line.n2p(a9 - width_now()), line.n2p(a9)
            r = Rectangle(width=x1[0] - x0[0], height=1.9, stroke_width=0)
            return r.set_fill(P.ORANGE, opacity=0.28).move_to([(x0[0] + x1[0]) / 2, x0[1] + 0.95, 0])

        sea = always_redraw(band)
        read = always_redraw(lambda: VGroup(
            layout.label(f"orbit eccentricity  e = {ecc.get_value():.2f}", font_size=20, color=P.FG),
            layout.label(f"overlapped, chaotic zone: {width_now():.0f} AU wide", font_size=20,
                         color=P.ORANGE),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT).move_to([-3.3, 2.3, 0]))
        cap2 = layout.caption("More eccentric orbits feel wider resonances; where neighbours "
                              "overlap, motion turns chaotic", font_size=20)
        self.play(FadeIn(sea), FadeIn(read), FadeOut(cap), FadeIn(cap2))
        self.play(ecc.animate.set_value(0.95), run_time=4.5, rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap2, settle=0.6)
        self.play(FadeOut(VGroup(line, xl, p9, p9_lab, ticks, names, sea, read, cap2)),
                  run_time=0.8)

        # ---- 2. the chaos map over (a, e) ----------------------------------
        ax, labels = widgets.labeled_axes(
            [150, 500, 50], [0, 0.95, 0.2], x_label="semi-major axis a (AU)",
            y_label="eccentricity e", y_rotate=True, numbers=True,
            x_length=8.3, y_length=4.3, shift_down=-0.5)
        ax.shift(LEFT * 1.9)
        labels.shift(LEFT * 1.9)
        self.play(Create(ax), FadeIn(labels), run_time=0.9)

        mp = d["map"]
        av, ev = np.asarray(mp["a"]), np.asarray(mp["e"])
        cells = np.asarray(mp["cells"])
        da, de = av[1] - av[0], ev[1] - ev[0]
        shade = {1: VGroup(), 2: VGroup()}
        for je in range(len(ev)):
            ia = 0
            while ia < len(av):
                c = cells[ia, je]
                if c == 0:
                    ia += 1
                    continue
                ib = ia
                while ib + 1 < len(av) and cells[ib + 1, je] == c:
                    ib += 1
                p0 = ax.c2p(av[ia] - da / 2, ev[je] - de / 2)
                p1 = ax.c2p(av[ib] + da / 2, ev[je] + de / 2)
                r = Rectangle(width=p1[0] - p0[0], height=p1[1] - p0[1], stroke_width=0)
                shade[c].add(r.set_fill(P.ORANGE, opacity=0.42 if c == 1 else 0.22)
                             .move_to((p0 + p1) / 2))
                ia = ib + 1

        b10 = next(b for b in d["boundaries"] if b["mass_earth"] == m9)
        edge = _boundary(ax, d["a_line"], b10["e_crit"], a9, b10["circular_zone_au"])
        edge.set_stroke(P.ORANGE, width=3)
        t_sea = layout.label("chaotic sea", font_size=17, color=P.ORANGE, weight="BOLD")
        t_sea.move_to(ax.c2p(462, 0.55))
        t_reg = layout.label("regular: smooth, predictable orbits", font_size=17, color=P.FG)
        t_reg.move_to(ax.c2p(285, 0.35))
        t_nep = layout.label("Neptune's resonances (q ≲ 40 AU)", font_size=14, color=P.ORANGE)
        t_nep.next_to(ax.c2p(245, 0.95), UP, buff=0.08)
        cap3 = layout.caption(f"Where motion is chaotic for a {m9:.0f} M⊕ Planet Nine at "
                              f"{a9:.0f} AU", font_size=21)
        self.play(FadeIn(shade[1]), FadeIn(shade[2]), Create(edge), FadeIn(cap3), run_time=1.4)
        self.play(FadeIn(t_sea), FadeIn(t_reg), FadeIn(t_nep))
        timing.hold_to_read(self, cap3, settle=0.4)

        # mass sweep: the boundary for 5, 10 and 20 Earth masses
        side = VGroup()
        curves = VGroup()
        for b in d["boundaries"]:
            cv = _boundary(ax, d["a_line"], b["e_crit"], a9, b["circular_zone_au"])
            op = 1.0 if b["mass_earth"] == m9 else 0.55
            cv.set_stroke(P.ORANGE, width=2.2, opacity=op)
            top = min((a for a, e in zip(d["a_line"], b["e_crit"]) if e is not None and e <= 0.95))
            tag = layout.label(f"{b['mass_earth']:.0f}", font_size=15, color=P.ORANGE)
            tag.next_to(ax.c2p(top, 0.95), UP, buff=0.08)
            curves.add(VGroup(cv, tag))
            side.add(layout.label(
                f"{b['mass_earth']:.0f} M⊕:  {100 * b['chaotic_fraction']:.0f}% chaotic",
                font_size=17, color=P.ORANGE if b["mass_earth"] == m9 else P.FG))
        head = layout.label("Planet Nine mass (M⊕)", font_size=17, color=P.FG, weight="BOLD")
        panel = VGroup(head, *side).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
        panel.next_to(ax, RIGHT, buff=0.5).shift(UP * 1.2)
        cap4 = layout.caption("The chaotic sea grows with Planet Nine's mass", font_size=21)
        self.play(FadeOut(cap3), FadeIn(cap4), FadeIn(head))
        for cv, lab in zip(curves, side):
            self.play(Create(cv), FadeIn(lab), run_time=0.8)
        timing.hold_to_read(self, cap4, panel, settle=0.4)

        # ---- 3. the known distant objects on the map ----------------------
        dots = VGroup()
        ring = VGroup()
        for o in d["etnos"]:
            chaotic = o["chaotic_p9"] or o["chaotic_neptune"]
            dt = Dot(ax.c2p(o["a"], o["e"]), radius=0.07, color=P.GREEN).set_z_index(5)
            dots.add(dt)
            if chaotic:
                ring.add(Circle(radius=0.17, color=P.GREEN, stroke_width=2.5).move_to(dt))
                ring.add(layout.label(o["name"], font_size=15, color=P.GREEN)
                         .next_to(dt, RIGHT, buff=0.22))
        tally = paper.result_readout(
            "known objects on regular orbits",
            f"{d['n_etno_regular']} of {d['n_etno']}", color=P.GREEN).scale(0.7)
        tally.next_to(panel, DOWN, buff=0.5).to_edge(RIGHT, buff=0.35)
        cap5 = layout.caption("The known distant objects: most sit in the regular zone, "
                              "one wanders in the chaotic sea", font_size=20)
        self.play(FadeIn(dots, lag_ratio=0.1), FadeOut(cap4), FadeIn(cap5), run_time=1.2)
        self.play(Create(ring), FadeIn(tally))
        timing.hold_to_read(self, cap5, tally, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "Close to Planet Nine, overlapping resonances make distant orbits chaotic.")
