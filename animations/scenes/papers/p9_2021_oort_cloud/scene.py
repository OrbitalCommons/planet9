"""Batygin & Brown (2021) -- injection of inner Oort cloud objects by Planet Nine.

Inner Oort cloud orbits have perihelia far beyond Neptune and would stay frozen
forever; Planet Nine's slow secular pull walks some perihelia down into the
distant Kuiper belt, so the observed sample is a mix of two populations. The
injected ones are herded into anti-alignment less strongly. Reproduced in
p9-2021-oort-cloud at reduced scale (mass-boosted secular Planet Nine, five
seeds); every particle track and every fraction shown is the crate's own
(anim.json -> papers -> p9-2021-oort-cloud).
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    AnnularSector,
    Arrow,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Rectangle,
    Scene,
    ValueTracker,
    VGroup,
    always_redraw,
    smooth,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2021-oort-cloud"


def rose(dvarpi_deg, centre, radius, color, n_bins=12):
    """Circular histogram of Δϖ: Planet Nine's perihelion points right (0°),
    anti-aligned orbits fill the left half."""
    v = np.asarray(dvarpi_deg) % 360.0
    counts, edges = np.histogram(v, bins=n_bins, range=(0, 360))
    frac = counts / max(counts.max(), 1)
    g = VGroup()
    for k, f in enumerate(frac):
        if f <= 0:
            continue
        g.add(AnnularSector(inner_radius=0, outer_radius=radius * f,
                            start_angle=np.deg2rad(edges[k]), angle=np.deg2rad(360 / n_bins),
                            color=color, fill_opacity=0.7, stroke_width=1,
                            stroke_color=P.BG).shift(centre))
    ring = Circle(radius=radius, color=P.MUTED, stroke_width=1).move_to(centre)
    half = AnnularSector(inner_radius=0, outer_radius=radius * 1.08, start_angle=np.pi / 2,
                         angle=np.pi, color=P.MUTED, fill_opacity=0.10,
                         stroke_width=0).shift(centre)
    p9 = Arrow(centre, centre + RIGHT * radius * 1.35, buff=0, color=P.BLUE, stroke_width=4,
               max_tip_length_to_length_ratio=0.15)
    return VGroup(half, ring, g, p9)


class OortCloud2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        ps = d["particles"]
        q_cut = d["q_injection_au"]
        self.add(paper.scene_header(CRATE))

        # 1. perihelia walked down into the Kuiper belt
        ax, labels = widgets.labeled_axes(
            [800, 2600, 300], [0, 300, 50], x_label="semi-major axis a (AU)",
            y_label="perihelion q (AU)", y_rotate=True, numbers=True,
            x_length=8.6, y_length=4.5, shift_down=0)
        VGroup(ax, labels).move_to([-1.6, 0.25, 0])
        belt = Rectangle(width=ax.c2p(2600, 0)[0] - ax.c2p(800, 0)[0],
                         height=ax.c2p(0, q_cut)[1] - ax.c2p(0, 0)[1], stroke_width=0)
        belt.set_fill(P.GREEN, opacity=0.08).move_to(
            (np.array(ax.c2p(800, 0)) + np.array(ax.c2p(2600, q_cut))) / 2)
        cut = DashedLine(ax.c2p(800, q_cut), ax.c2p(2600, q_cut), color=P.GREEN, stroke_width=2)
        belt_l = layout.label(f"distant Kuiper belt, q < {q_cut:.0f} AU", font_size=16,
                              color=P.GREEN).next_to(ax.c2p(2600, q_cut / 2), RIGHT, buff=0.15)
        cloud_l = layout.label("inner Oort cloud:\nfar beyond\nNeptune's reach", font_size=16,
                               color=P.MUTED).next_to(ax.c2p(2600, 220), RIGHT, buff=0.15)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(belt), Create(cut), FadeIn(belt_l))

        s = ValueTracker(0.0)
        a = np.array([p["a_au"] for p in ps])
        q0 = np.array([p["q0_au"] for p in ps])
        qm = np.array([p["q_min_au"] for p in ps])
        shown = [k for k, p in enumerate(ps) if p["injectable"]]

        def cloud():
            g = VGroup()
            f = s.get_value()
            for k in shown:
                q = q0[k] + f * (qm[k] - q0[k])
                col = P.GREEN if q < q_cut else P.MUTED
                g.add(Dot(ax.c2p(a[k], q), radius=0.045, color=col))
            return g

        dots = always_redraw(cloud)
        cap = layout.caption(f"{len(shown)} inner Oort cloud orbits, "
                             f"perihelia well beyond {q_cut:.0f} AU", font_size=22)
        self.play(FadeIn(dots), FadeIn(cloud_l), FadeIn(cap))
        timing.hold_to_read(self, cap, settle=0.4)
        cap2 = layout.caption(f"Add Planet Nine and run {d['equivalent_gyr']:.0f} Gyr: "
                              "its slow pull walks perihelia inward", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2))
        self.play(s.animate.set_value(1.0), run_time=5.0, rate_func=smooth)
        dots.clear_updaters()
        n_in, n_ok = d["n_injected"], d["n_injectable"]
        tally = paper.result_readout("injected into the Kuiper belt",
                                     f"{n_in} of {n_ok}", color=P.GREEN).scale(0.8)
        tally.move_to([5.1, 1.7, 0])
        ctrl = layout.label(f"without Planet Nine: {d['n_control_crossed']} of {n_ok}",
                            font_size=17, color=P.RED).next_to(tally, DOWN, buff=0.25)
        self.play(FadeOut(cloud_l), run_time=0.4)
        self.play(FadeIn(tally))
        self.play(FadeIn(ctrl))
        cap3 = layout.caption("So the distant sample mixes Kuiper belt objects "
                              "with objects from the Oort cloud", font_size=22)
        self.play(FadeOut(cap2), FadeIn(cap3))
        timing.hold_to_read(self, cap3, tally, ctrl, settle=0.8)
        self.play(FadeOut(VGroup(ax, labels, belt, cut, belt_l, dots, tally, ctrl, cap3)))

        # 2. how well each population is herded
        r = 1.45
        left, right = np.array([-3.4, 0.35, 0]), np.array([2.4, 0.35, 0])
        rose_kb = rose(d["dvarpi_scattered_deg"], left, r, P.GREEN)
        rose_ioc = rose(d["dvarpi_ioc_deg"], right, r, P.GREEN)
        p9l = layout.label("Planet Nine's\nperihelion", font_size=15, color=P.BLUE)
        p9l.next_to(left + RIGHT * r * 1.35, DOWN, buff=0.12)
        anti = layout.label("anti-aligned", font_size=15, color=P.MUTED)
        anti.next_to(left + LEFT * r * 1.08, LEFT, buff=0.1)
        pub = d["published"]

        def block(title, n, f, f_pub, centre):
            return VGroup(
                layout.label(title, font_size=19, weight="BOLD"),
                layout.label(f"{n} orbits", font_size=15, color=P.MUTED),
            ).arrange(DOWN, buff=0.08).next_to(centre + UP * r, UP, buff=0.3), VGroup(
                layout.label(f"anti-aligned {100 * f:.0f}%", font_size=20, color=P.GREEN,
                             weight="BOLD"),
                layout.label(f"paper {100 * f_pub:.0f}%", font_size=16, color=P.FG),
            ).arrange(DOWN, buff=0.08).next_to(centre + DOWN * r, DOWN, buff=0.3)

        t_kb, v_kb = block("from the Kuiper belt", len(d["dvarpi_scattered_deg"]),
                           d["f_varpi_scattered"], pub["f_varpi_scattered"], left)
        t_ioc, v_ioc = block("injected from the Oort cloud", len(d["dvarpi_ioc_deg"]),
                             d["f_varpi_ioc"], pub["f_varpi_ioc"], right)
        cap4 = layout.caption("Where do the perihelia end up, relative to Planet Nine's?",
                              font_size=22)
        self.play(FadeIn(rose_kb), FadeIn(t_kb), FadeIn(p9l), FadeIn(anti), FadeIn(cap4))
        self.play(FadeIn(v_kb))
        self.play(FadeIn(rose_ioc), FadeIn(t_ioc))
        self.play(FadeIn(v_ioc))
        timing.hold_to_read(self, cap4, v_kb, v_ioc, settle=0.4)
        cap5 = layout.caption(
            f"Paper: Oort cloud recruits are herded less tightly, "
            f"{100 * pub['f_varpi_ioc']:.0f}% against {100 * pub['f_varpi_scattered']:.0f}%",
            font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, settle=0.8)
        cap6 = layout.caption("This reduced model (no giant planets) keeps that order, "
                              "but both land lower", font_size=22)
        self.play(FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=1.0)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "A weaker-herded Oort component means a more eccentric Planet Nine.")
