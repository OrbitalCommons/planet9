"""Hu, Huang, Gladman & Zhu (2025) -- early stellar flybys are unlikely.

A close stellar flyby can lift TNO perihelia to make sednoids, but the LOW
inclinations (i < 30 deg) and primordial apsidal alignment of detached extreme
TNOs only allow special (coplanar / ecliptic-symmetric) geometries. Folding that
small geometric fraction into how often a close encounter (q* <= 1000 au) happens
leaves <= 5% odds. Reproduced in p9-2025-stellar-flybys.
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
    RIGHT,
    Scene,
    UP,
    VGroup,
    Write,
    rate_functions,
)

import p9_manim as P
from p9_manim import layout, orbits, paper, timing

CRATE = "p9-2025-stellar-flybys"


class StellarFlybys2025(Scene):
    def construct(self):
        tb = paper.title_block(
            CRATE, "Did a passing star make the sednoids?",
            cite="Hu, Huang, Gladman & Zhu (2025)")
        self.play(Write(tb[0]), FadeIn(tb[1]))
        self.play(tb.animate.scale(0.62).to_edge(UP, buff=0.3))

        sun = orbits.sun()
        self.add(sun)

        # the primordial scattered/detached disk
        rng = np.random.default_rng(2025)
        start = np.deg2rad(40) + rng.normal(0, np.deg2rad(20), 8)
        disk = VGroup(*[
            orbits.ellipse_orbit(2.1, 0.62, color=P.GREEN, varpi=v,
                                 stroke_width=1.6, opacity=0.8)
            for v in start
        ])
        self.play(Create(disk, lag_ratio=0.1), run_time=1.6)
        disk_lbl = layout.label("scattered/detached disk", font_size=18,
                                color=P.MUTED).to_edge(DOWN, buff=0.6)
        self.play(FadeIn(disk_lbl))
        timing.hold_to_read(self, disk_lbl, settle=0.5)
        self.play(FadeOut(disk_lbl))

        # a star flies past on a straight path, perturbing perihelia
        path = Line(LEFT * 6.2 + UP * 2.4, RIGHT * 6.2 + DOWN * 2.4)
        star = Dot(path.get_start(), radius=0.11, color=P.SUN).set_z_index(6)
        star_lbl = layout.label("passing field star", font_size=18, color=P.SUN)
        star_lbl.add_updater(lambda m: m.next_to(star, UP, buff=0.12))
        self.play(FadeIn(star), FadeIn(star_lbl))
        # the flyby tugs the disk apsides toward a common direction (alignment)
        target = np.deg2rad(35) + rng.normal(0, np.deg2rad(7), 8)
        self.play(
            star.animate.move_to(path.get_end()),
            *orbits.precess(disk, start, target),
            run_time=3.0, rate_func=rate_functions.smooth,
        )
        star_lbl.clear_updaters()
        self.play(FadeOut(star), FadeOut(star_lbl))

        nudge = layout.label("perihelia lifted + apsides aligned", font_size=18,
                             color=P.ORANGE).to_edge(DOWN, buff=0.6)
        self.play(FadeIn(nudge))
        timing.hold_to_read(self, nudge, settle=0.6)
        self.play(FadeOut(nudge), FadeOut(disk))

        # BUT: the detached eTNOs are observed at LOW inclination
        con = layout.label(
            "but detached extreme TNOs sit at LOW inclination: i < 30 deg",
            font_size=22, color=P.RED).to_edge(UP, buff=1.5)
        self.play(FadeIn(con))
        timing.hold_to_read(self, con, settle=0.6)
        geom = layout.label(
            "only special flyby geometries fit -- coplanar or ecliptic-symmetric",
            font_size=20, color=P.TEAL).next_to(con, DOWN, buff=0.3)
        self.play(FadeIn(geom))
        timing.hold_to_read(self, geom, settle=0.8)
        self.play(FadeOut(con), FadeOut(geom))

        # the focused encounter rate sets how OFTEN a close passage happens
        eq1 = layout.explain_equation(
            self,
            [r"\Gamma", "=", r"n_\star\,\sigma\,v", r"\quad",
             r"\sigma=\pi q_\star^2\!\left(1+\tfrac{2GM}{q_\star v^2}\right)"],
            [
                (0, "rate of close stellar encounters in the birth cluster"),
                (2, "star density x cross-section x relative speed"),
                (4, "gravitationally focused area for passages q* <= 1000 au"),
            ],
            scale=0.78,
        )
        self.play(FadeOut(eq1))

        # the overall probability: encounter happens x geometry fits x it works
        eq2 = layout.explain_equation(
            self,
            [r"P", r"\approx", r"P_{\rm enc}", r"\times",
             r"f_{\rm geom}", r"\times", r"f_{\rm succ}"],
            [
                (2, "chance a close encounter (q* <= 1000 au) ever occurred"),
                (4, "fraction of geometries giving i < 30 deg + alignment"),
                (6, "fraction of those that actually build the sednoids"),
            ],
            scale=0.92,
            where=UP * 1.7,
        )
        self.play(FadeOut(eq2))

        ro = paper.result_readout(
            "probability a constraint-satisfying flyby occurred",
            "<= 5%", color=P.RED).scale(0.8)
        self.play(FadeIn(ro))
        timing.hold_to_read(self, ro, settle=1.0)
        self.play(FadeOut(ro))

        layout.show_takeaway(
            self,
            "An early flyby clears every hurdle only ~5% of the time.")
