"""Kaib et al. (2019) -- OSSOS XV: probing the distant solar system with
observed scattering TNOs.

Scattering objects graze Neptune at perihelion and get kicked around in
semi-major axis; lift the perihelion a few AU and they detach, frozen out of
Neptune's reach. A distant planet does the lifting, so the objects still
scattering today constrain it. The paper does this with N-body runs and the
inclinations of 69 OSSOS objects; the crate reproduces the perihelion side of
the argument: the scattering/detached divide from resonance overlap and the
secular lift a revised Planet Nine gives a synthetic scattering population
(anim.json -> papers -> p9-2019-ossos-scattering).
"""
from manim import (
    DOWN,
    LEFT,
    UP,
    Arrow,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Polygon,
    Rectangle,
    Scene,
    Transform,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2019-ossos-scattering"


class OssosScattering2019(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        bnd = d["boundary"]
        p9 = d["p9"]

        self.add(paper.scene_header(CRATE))

        # 1. the divide between scattering and detached orbits
        q_top = 200
        ax, labels = widgets.labeled_axes(
            [0, 1000, 100], [0, q_top, 40], x_label="semi-major axis (AU)",
            y_label="perihelion distance (AU)", y_rotate=True, numbers=True,
            x_length=10.6, y_length=4.6, shift_down=-0.35)
        pts = [(b["a"], b["q_crit"]) for b in bnd if b["q_crit"] > 0]
        edge = widgets.curve(ax, [a for a, _ in pts], [q for _, q in pts], color=P.ORANGE)
        below = Polygon(ax.c2p(pts[0][0], 0), *[ax.c2p(a, q) for a, q in pts],
                        ax.c2p(pts[-1][0], 0), stroke_width=0)
        below.set_fill(P.ORANGE, opacity=0.12)
        nep = DashedLine(ax.c2p(0, 30), ax.c2p(1000, 30), color=P.MUTED, stroke_width=1.6)
        nep_lab = layout.label("Neptune's distance, 30 AU", font_size=14, color=P.MUTED)
        nep_lab.next_to(ax.c2p(120, 30), DOWN, buff=0.08)
        pub = DashedLine(ax.c2p(0, d["published_q_boundary"]),
                         ax.c2p(1000, d["published_q_boundary"]), color=P.FG, stroke_width=1.2)
        pub.set_stroke(opacity=0.5)
        pub_lab = layout.label(f"paper: detach near q ≈ {d['published_q_boundary']:.0f} AU",
                               font_size=13, color=P.FG).set_opacity(0.7)
        pub_lab.next_to(ax.c2p(120, d["published_q_boundary"]), UP, buff=0.06)
        scat_lab = layout.label("scattering: Neptune kicks the orbit around", font_size=15,
                                color=P.ORANGE)
        scat_lab.move_to(ax.c2p(760, 18))
        det_lab = layout.label("detached: out of Neptune's reach", font_size=15, color=P.FG)
        det_lab.move_to(ax.c2p(300, 140))
        dots = VGroup(*[Dot(ax.c2p(o["a"], o["q"]), radius=0.05, color=P.FG) for o in objs])
        cap = layout.caption(
            "Whether Neptune can still grab an orbit depends on how close its perihelion comes",
            font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(Create(nep), FadeIn(nep_lab))
        self.play(Create(edge), FadeIn(below), FadeIn(scat_lab), FadeIn(det_lab),
                  FadeIn(cap), run_time=1.6)
        self.play(Create(pub), FadeIn(pub_lab))
        timing.hold_to_read(self, cap, settle=0.8)
        cap2 = layout.caption(
            f"{d['n_objects']} synthetic scattering objects, all just inside the divide",
            font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(FadeIn(dots, lag_ratio=0.03), FadeIn(cap2), run_time=1.4)
        timing.hold_to_read(self, cap2, settle=0.8)

        # 2. add a distant planet: perihelia are lifted
        lifted = VGroup(*[
            Dot(ax.c2p(o["a"], min(o["q_lifted"], q_top)), radius=0.05,
                color=P.BLUE if o["detached"] else P.FG) for o in objs])
        trails = VGroup(*[
            Arrow(ax.c2p(o["a"], o["q"]), ax.c2p(o["a"], min(o["q_lifted"], q_top)), buff=0.04,
                  color=P.BLUE, stroke_width=1.6, max_tip_length_to_length_ratio=0.08)
            .set_opacity(0.45) for o in objs])
        n_det = sum(1 for o in objs if o["detached"])
        tag = VGroup(
            layout.label(f"with a {p9['mass']:.0f} Earth-mass planet at {p9['a']:.0f} AU",
                         font_size=16, color=P.BLUE),
            layout.label(f"after {d['sculpting_myr']:.0f} million years:  {n_det} of "
                         f"{len(objs)} detached", font_size=16, color=P.BLUE, weight="BOLD"),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        tag.move_to(ax.c2p(30, 175), aligned_edge=LEFT)
        cap3 = layout.caption(
            "The planet's slow tug raises their perihelia out of Neptune's reach",
            font_size=22)
        self.play(FadeOut(cap2), FadeOut(det_lab), run_time=0.4)
        self.play(FadeIn(tag[0]), FadeIn(cap3))
        self.play(Transform(dots, lifted), LaggedStart(*[Create(t) for t in trails],
                                                      lag_ratio=0.02), run_time=3.0)
        self.play(FadeIn(tag[1]))
        timing.hold_to_read(self, cap3, tag, settle=1.2)
        self.play(FadeOut(VGroup(ax, labels, edge, below, nep, nep_lab, pub, pub_lab, scat_lab,
                                 dots, trails, tag, cap3)))

        # 3. how that depends on the planet's mass
        sweep = d["sweep"]
        bx0, bw, gap = -5.0, 1.3, 0.6
        base = -1.8
        h = 3.6
        axis = VGroup(
            DashedLine([bx0 - 0.4, base + h, 0], [bx0 + 4 * (bw + gap), base + h, 0],
                       color=P.MUTED, stroke_width=1.2),
            layout.label("all of them", font_size=13, color=P.MUTED)
            .next_to([bx0 - 0.4, base + h, 0], LEFT, buff=0.1),
            DashedLine([bx0 - 0.4, base, 0], [bx0 + 4 * (bw + gap), base, 0], color=P.MUTED,
                       stroke_width=1.2),
            layout.label("none", font_size=13, color=P.MUTED)
            .next_to([bx0 - 0.4, base, 0], LEFT, buff=0.1),
        )
        bars, texts = VGroup(), VGroup()
        for k, s in enumerate(sweep):
            x = bx0 + k * (bw + gap) + bw / 2
            hh = max(h * s["detached"], 0.02)
            r = Rectangle(width=bw, height=hh, stroke_width=0)
            r.set_fill(P.BLUE, opacity=0.75).move_to([x, base + hh / 2, 0])
            bars.add(r)
            unit = "Earth mass" if s["mass"] == 1 else "Earth masses"
            texts.add(layout.label(f"{s['mass']:g} {unit}", font_size=15, color=P.FG)
                      .move_to([x, base - 0.3, 0]))
            texts.add(layout.label(f"{100 * s['detached']:.0f}%", font_size=17, color=P.BLUE,
                                   weight="BOLD").next_to(r, UP, buff=0.1))
        side = VGroup(
            layout.label("share of the scattering objects", font_size=16, color=P.FG),
            layout.label(f"detached within {d['sculpting_myr']:.0f} Myr", font_size=16,
                         color=P.FG),
            layout.label(f"planet at a = {p9['a']:.0f} AU, e = {p9['e']:.2f}", font_size=14,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        side.move_to([bx0 + 4 * (bw + gap) + 0.4, base + h - 0.4, 0], aligned_edge=LEFT)
        cap4 = layout.caption(
            "Even a planet of a few Earth masses reshapes them: scattering objects are a sharp probe",
            font_size=22)
        self.play(FadeIn(axis), FadeIn(side))
        self.play(LaggedStart(*[FadeIn(b, shift=UP * 0.2) for b in bars], lag_ratio=0.2),
                  FadeIn(texts, lag_ratio=0.1), FadeIn(cap4), run_time=1.6)
        timing.hold_to_read(self, cap4, side, settle=1.4)
        self.play(FadeOut(cap4))

        layout.show_takeaway(
            self, "Objects still scattering today limit how hard a distant planet can pull.")
