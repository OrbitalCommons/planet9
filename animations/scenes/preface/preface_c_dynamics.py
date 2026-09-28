"""Preface C -- The clue and the dynamics: clustering, precession, resonance.

Scenes: P05Clustering, P06Precession, P07Resonance.

Every number and curve comes from ``anim.json -> preface -> c_dynamics``
(crates/p9-anim-data/src/preface/c_dynamics.rs).
"""
import numpy as np
from manim import (
    Arrow,
    ORIGIN,
    Circle,
    Create,
    DOWN,
    DashedLine,
    DL,
    DashedVMobject,
    DecimalNumber,
    Dot,
    FadeIn,
    FadeOut,
    GrowArrow,
    GrowFromCenter,
    LEFT,
    LaggedStart,
    Line,
    MathTex,
    RIGHT,
    Scene,
    Text,
    Transform,
    UP,
    VGroup,
    VMobject,
    ValueTracker,
    Write,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, timing, widgets


def _data():
    return dataio.section("preface")["c_dynamics"]


def _header(scene, badge, title):
    scene.add(layout.concept_badge(badge))
    t = Text(title, color=P.FG, font_size=30, weight="BOLD").to_edge(UP, buff=0.55)
    scene.play(Write(t), run_time=1.0)
    return t


def _say(scene, text, hold=True, settle=0.6):
    cap = layout.caption(text)
    scene.play(FadeIn(cap, shift=UP * 0.1), run_time=0.5)
    if hold:
        timing.hold_to_read(scene, cap, settle=settle)
    return cap


def _unsay(scene, *caps):
    scene.play(*[FadeOut(c) for c in caps], run_time=0.4)


def _u(angle):
    return np.array([np.cos(angle), np.sin(angle), 0.0])


def _mapped_orbit(a, e, varpi, k, center, color, stroke_width=2.0, opacity=0.9, n=260):
    """A Kepler orbit drawn with radius compressed as k*sqrt(r): the true
    direction of every point (so the perihelion direction) is preserved."""
    pts = []
    for nu in np.linspace(-np.pi, np.pi, n):
        r = a * (1 - e * e) / (1 + e * np.cos(nu))
        pts.append(center + k * np.sqrt(r) * _u(nu + varpi))
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners(pts)
    return m.set_stroke(opacity=opacity)


def _chain(angles, length, start, color, stroke_width=3.0, center=None):
    """Unit vectors laid head to tail from ``start`` (or, with ``center``,
    shifted so the whole walk is centred there); returns (arrows, start, end)."""
    pts = [np.array(start, dtype=float)]
    for ang in angles:
        pts.append(pts[-1] + length * _u(ang))
    pts = np.array(pts)
    if center is not None:
        mid = 0.5 * (pts.min(axis=0) + pts.max(axis=0))
        pts = pts + (np.array(center, dtype=float) - mid)
    g = VGroup(*[
        Arrow(p, q, buff=0, color=color, stroke_width=stroke_width,
              max_tip_length_to_length_ratio=0.3, tip_length=0.13)
        for p, q in zip(pts[:-1], pts[1:])
    ])
    return g, pts[0], pts[-1]


def _plot_axes(x_range, y_range, x_len, y_len, center, xticks, yticks, xlabel, ylabel):
    """Film axes with legible hand-placed tick labels. Returns (ax, labels)."""
    ax = widgets.axes(x_range, y_range, x_length=x_len, y_length=y_len, shift_down=0)
    ax.move_to(center)
    xt = VGroup(*[layout.label(t, font_size=16, color=P.FG)
                  .next_to(ax.c2p(v, y_range[0]), DOWN, buff=0.12) for v, t in xticks])
    yt = VGroup(*[layout.label(t, font_size=16, color=P.FG)
                  .next_to(ax.c2p(x_range[0], v), LEFT, buff=0.12) for v, t in yticks])
    xl = layout.label(xlabel, font_size=18, color=P.FG).next_to(xt, DOWN, buff=0.12)
    yl = layout.label(ylabel, font_size=18, color=P.FG).rotate(np.pi / 2).next_to(yt, LEFT, buff=0.15)
    return ax, VGroup(xt, yt, xl, yl)


class P05Clustering(Scene):
    """Learning goal: the distant orbits point the same way, and ten random
    directions line up that well only ~1% of the time."""

    def construct(self):
        d = _data()["clustering"]
        objs = d["objects"]
        n = len(objs)
        _header(self, "PREFACE 05", "The clue: distant orbits point the same way")

        # ---- beat 1: the real orbits, seen from above -----------------------
        C = np.array([-3.3, -0.1, 0.0])
        k = 2.75 / np.sqrt(max(o["a"] * (1 + o["e"]) for o in objs))
        sun = orbits.sun(radius=0.09).move_to(C)
        nep = Circle(radius=k * np.sqrt(30.07), color=P.FG, stroke_width=1.6).move_to(C)
        nep.set_stroke(opacity=0.8)
        nep_lbl = layout.label("Neptune, 30 AU", font_size=16, color=P.FG)
        nep_lbl.move_to(C + np.array([1.55, -0.95, 0]))
        nep_tick = Line(nep_lbl.get_left() + LEFT * 0.05, C + k * np.sqrt(30.07) * _u(-0.6),
                        color=P.FG, stroke_width=1).set_stroke(opacity=0.6)
        scale_note = layout.label("distance from Sun drawn as √r", font_size=15, color=P.MUTED)
        scale_note.move_to(np.array([-6.0, -2.75, 0]), aligned_edge=LEFT)
        self.play(FadeIn(sun), Create(nep), FadeIn(nep_lbl), Create(nep_tick), FadeIn(scale_note))

        cap = _say(self, "The ten most distant well-measured orbits (a > 230 AU), from above", hold=False)
        orbs = VGroup(*[
            _mapped_orbit(o["a"], o["e"], np.deg2rad(o["varpi_deg"]), k, C, P.GREEN, 1.8, 0.75)
            for o in objs
        ])
        self.play(LaggedStart(*[Create(o) for o in orbs], lag_ratio=0.3), run_time=4.0)
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)

        peri = VGroup(*[
            Dot(C + k * np.sqrt(o["q"]) * _u(np.deg2rad(o["varpi_deg"])), radius=0.06, color=P.GREEN)
            for o in objs
        ])
        cap = _say(self, "Mark each perihelion (closest point to the Sun): all on one side", hold=False)
        self.play(LaggedStart(*[GrowFromCenter(p) for p in peri], lag_ratio=0.1), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.5)
        _unsay(self, cap)

        # ---- beat 2: keep only the direction ϖ -------------------------------
        R = 2.3
        varpis = [np.deg2rad(o["varpi_deg"]) for o in objs]
        dial = VGroup(*[
            Arrow(C, C + R * _u(v), buff=0, color=P.GREEN, stroke_width=3, tip_length=0.16)
            for v in varpis
        ])
        dial_ring = Circle(radius=R, color=P.MUTED, stroke_width=1).move_to(C).set_stroke(opacity=0.5)
        cap = _say(self, "Keep only the direction each orbit points: its longitude of perihelion ϖ",
                   hold=False)
        self.play(orbs.animate.set_stroke(opacity=0.15), FadeOut(peri), Create(dial_ring),
                  LaggedStart(*[GrowArrow(a) for a in dial], lag_ratio=0.08), run_time=1.8)
        timing.hold_to_read(self, cap, settle=0.4)
        _unsay(self, cap)
        self.play(FadeOut(orbs), run_time=0.5)
        left_head = layout.label("real distant orbits", font_size=18, color=P.GREEN)
        left_head.move_to(np.array([-3.3, 2.75, 0]))
        self.play(FadeIn(left_head))

        # ---- beat 3: head to tail ---------------------------------------------
        L = 0.45
        mid = np.array([3.0, 0.1, 0.0])
        chain, S, end = _chain(varpis, L, ORIGIN, P.GREEN, 3.0, center=mid)
        right_head = layout.label("add the arrows head to tail", font_size=18, color=P.FG)
        right_head.move_to(np.array([3.6, 2.75, 0]))
        cap = _say(self, "Add the arrows head to tail: aligned arrows travel far", hold=False)
        self.play(FadeIn(right_head))
        copies = [dial[i].copy() for i in range(n)]
        self.play(LaggedStart(*[Transform(copies[i], chain[i]) for i in range(n)],
                              lag_ratio=0.25), run_time=3.2)
        self.remove(*copies)
        self.add(chain)
        res = Arrow(S, end, buff=0, color=P.GREEN, stroke_width=6, tip_length=0.22)
        self.play(GrowArrow(res))
        rbar = MathTex(r"\bar R = \frac{\text{length}}{N} = " + f"{d['r_bar']:.2f}", color=P.GREEN)
        rbar.scale(0.75).move_to(np.array([5.2, -1.2, 0]))
        self.play(Write(rbar))
        timing.hold_to_read(self, cap, rbar, settle=0.6)
        _unsay(self, cap)
        scale_lbl = layout.label("0 = scattered, 1 = perfectly aligned", font_size=16, color=P.MUTED)
        scale_lbl.next_to(rbar, DOWN, buff=0.2)
        self.play(FadeIn(scale_lbl))
        timing.hold_to_read(self, scale_lbl, settle=0.3)

        # ---- beat 4: random skies ---------------------------------------------
        ghost = VGroup(chain, res)
        self.play(ghost.animate.set_opacity(0.18), FadeOut(rbar), FadeOut(scale_lbl), run_time=0.6)
        cap = _say(self, "Random directions wander like a drunkard's walk and rarely get far",
                   hold=False)
        prev = None
        for i, ex in enumerate(d["random_examples"]):
            angs = np.deg2rad(ex["varpi_deg"])
            ch, s0, en = _chain(angs, L, ORIGIN, P.RED, 2.5, center=mid)
            rr = Arrow(s0, en, buff=0, color=P.RED, stroke_width=5, tip_length=0.18)
            head = layout.label(f"random sky {i + 1}", font_size=18, color=P.RED).move_to(right_head)
            val = MathTex(rf"\bar R = {ex['r_bar']:.2f}", color=P.RED).scale(0.75)
            val.move_to(np.array([5.2, -1.2, 0]))
            grp = VGroup(ch, rr, val)
            outs = [FadeOut(prev)] if prev is not None else []
            self.play(*outs, Transform(right_head, head),
                      LaggedStart(*[GrowArrow(a) for a in ch], lag_ratio=0.15), run_time=1.6)
            self.play(GrowArrow(rr), FadeIn(val), run_time=0.6)
            self.wait(1.0)
            prev = grp
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)
        self.play(FadeOut(prev), FadeOut(ghost), FadeOut(right_head))

        # ---- beat 5: how often does chance do this well? ----------------------
        hist = np.array(d["mc_hist"], dtype=float)
        pct = 100.0 * hist / d["mc_trials"]
        edges = np.linspace(0, 1, len(hist) + 1)
        ax, ax_lbls = _plot_axes(
            [0, 1, 0.2], [0, 8, 2], 5.2, 3.3, np.array([3.75, 0.45, 0]),
            [(v, f"{v:.1f}") for v in (0, 0.2, 0.4, 0.6, 0.8, 1.0)],
            [(v, f"{v}%") for v in (0, 2, 4, 6, 8)],
            "R̄ of ten random directions", "share of random skies")
        cut = int(np.floor(d["r_bar"] * len(hist)))
        bars_lo = widgets.histogram(ax, edges[:cut + 1], pct[:cut], color=P.RED, opacity=0.35)
        bars_hi = widgets.histogram(ax, edges[cut:], pct[cut:], color=P.RED, opacity=0.95)
        mark = DashedLine(ax.c2p(d["r_bar"], 0), ax.c2p(d["r_bar"], 7.6), color=P.GREEN, stroke_width=3)
        mark_lbl = MathTex(rf"\bar R_{{\rm real}} = {d['r_bar']:.2f}", color=P.GREEN).scale(0.6)
        mark_lbl.next_to(mark.get_end(), RIGHT, buff=0.08)
        head = layout.label(f"{d['mc_trials']:,} random skies of {n} orbits", font_size=18, color=P.FG)
        head.move_to(np.array([3.6, 2.75, 0]))
        cap = _say(self, "Draw ten random directions again and again, and record R̄ each time", hold=False)
        self.play(FadeIn(head), Create(ax), FadeIn(ax_lbls))
        self.play(LaggedStart(*[GrowFromCenter(b) for b in bars_lo], lag_ratio=0.05),
                  LaggedStart(*[GrowFromCenter(b) for b in bars_hi], lag_ratio=0.05), run_time=2.2)
        self.play(Create(mark), FadeIn(mark_lbl))
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)

        frac = d["mc_frac_beat"]
        readout = layout.label(f"as aligned as the real ones: {100 * frac:.1f}%  (1 in {1 / frac:.0f})",
                               font_size=20, color=P.FG)
        readout.move_to(np.array([3.4, -2.3, 0]))
        self.play(bars_hi.animate.set_fill(P.ORANGE, opacity=1).set_stroke(P.ORANGE),
                  FadeIn(readout), run_time=0.8)
        timing.hold_to_read(self, readout, settle=0.8)
        cap = _say(self, "Caveat: this assumes the surveys looked everywhere equally (Preface 11)")
        _unsay(self, cap, scale_note)
        layout.show_takeaway(
            self, "Ten orbits sharing a direction happens by chance only ~1% of the time.")


class P06Precession(Scene):
    """Learning goal: the giant planets make every distant orbit precess, at
    a rate set by its own a and e, so an alignment left alone is erased within
    a few hundred Myr."""

    def construct(self):
        d = _data()["precession"]
        _header(self, "PREFACE 06", "Orbits turn: the giant planets make them precess")

        # ---- beat 1: giants blur into rings; a distant orbit's apse creeps ------
        C = np.array([-3.3, -0.1, 0.0])
        sed = d["objects"][0]
        a_s, e_s = sed["a"], sed["e"]
        k = 2.35 / np.sqrt(a_s * (1 + e_s))
        kz = 2.2 / 30.07  # zoomed-in, linear view of the giant planets
        sun = orbits.sun(radius=0.09).move_to(C)
        self.add(sun)
        giants = d["giants"]
        rings = VGroup(*[Circle(radius=kz * ag, color=P.FG, stroke_width=1.4).move_to(C)
                         .set_stroke(opacity=0.55) for _, ag in giants])
        tt = ValueTracker(0.0)
        dots = VGroup(*[
            always_redraw(lambda ag=ag, ph=j: Dot(
                C + kz * ag * _u(ph * 1.7 + 2 * np.pi * tt.get_value() * (5.203 / ag) ** 1.5),
                radius=0.07, color=P.FG))
            for j, (_, ag) in enumerate(giants)
        ])
        names = layout.label("Jupiter, Saturn, Uranus, Neptune", font_size=18, color=P.FG)
        names.next_to(C + kz * 30.07 * DOWN, DOWN, buff=0.15)
        cap = _say(self, "Jupiter, Saturn, Uranus and Neptune circle the Sun again and again...",
                   hold=False)
        self.play(Create(rings), FadeIn(dots), FadeIn(names))
        self.play(tt.animate.set_value(2.0), run_time=3.0, rate_func=rate_functions.linear)
        _unsay(self, cap)
        blur = VGroup(*[Circle(radius=kz * ag, color=P.ORANGE, stroke_width=9).move_to(C)
                        .set_stroke(opacity=0.45) for _, ag in giants])
        cap = _say(self, "...so over millions of years they act like rings of extra mass", hold=False)
        self.play(FadeOut(dots), Transform(rings, blur), run_time=1.2)
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        small = VGroup(*[Circle(radius=k * np.sqrt(ag), color=P.ORANGE, stroke_width=6).move_to(C)
                         .set_stroke(opacity=0.45) for _, ag in giants])
        zoom = layout.label("zoom out: distance from Sun drawn as √r", font_size=16, color=P.MUTED)
        zoom.move_to(np.array([-6.8, -2.85, 0]), aligned_edge=LEFT)
        self.play(FadeOut(names), Transform(rings, small), FadeIn(zoom), run_time=1.5)

        # A Sedna-shaped orbit whose apse creeps forward each lap (sped up).
        lap = ValueTracker(0.0)
        v0 = np.deg2rad(30.0)
        creep = np.deg2rad(22.0)  # per lap, for display only

        def body_pos():
            t = lap.get_value()
            nu = orbits.nu_from_mean_anomaly(e_s, 2 * np.pi * t)
            r = a_s * (1 - e_s ** 2) / (1 + e_s * np.cos(nu))
            return C + k * np.sqrt(r) * _u(nu + v0 + creep * t)

        live = always_redraw(lambda: _mapped_orbit(a_s, e_s, v0 + creep * lap.get_value(), k, C,
                                                   P.GREEN, 2.4, 0.95))
        body = always_redraw(lambda: Dot(body_pos(), radius=0.07, color=P.GREEN))
        apse = always_redraw(lambda: Arrow(
            C, C + k * np.sqrt(a_s * (1 - e_s)) * 2.4 * _u(v0 + creep * lap.get_value()),
            buff=0, color=P.ORANGE, stroke_width=4, tip_length=0.16))
        trail = VGroup()
        for j in range(5):
            trail.add(_mapped_orbit(a_s, e_s, v0 + creep * j, k, C, P.GREEN, 1.2, 0.25))
        cap = _say(self, "Their extra pull near perihelion makes a distant orbit's direction creep "
                         "forward", hold=False)
        self.play(Create(live), FadeIn(body), GrowArrow(apse))
        for j in range(4):
            self.add(trail[j])
            self.play(lap.animate.set_value(j + 1), run_time=1.4, rate_func=rate_functions.linear)
        self.add(trail[4])
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)
        sped = layout.label("(sped up enormously)", font_size=16, color=P.MUTED)
        sped.move_to(np.array([3.6, -0.9, 0]))
        facts = VGroup(
            layout.label("Sedna, for scale", font_size=20, color=P.GREEN, weight="BOLD"),
            layout.label(f"one orbit: {int(round(d['sedna']['orbit_yr'], -2)):,} yr", font_size=18,
                         color=P.FG),
            layout.label(f"one full turn of its apse: {d['sedna']['precession_myr'] / 1000:.1f} Gyr",
                         font_size=18, color=P.ORANGE),
            layout.label(f"≈ {d['sedna']['orbits_per_turn']:,.0f} orbits per turn", font_size=18,
                         color=P.FG),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.22).move_to(np.array([3.6, 0.6, 0]))
        cap = _say(self, "This is apsidal precession: slow, but the Solar System is 4.5 Gyr old",
                   hold=False)
        self.play(FadeIn(sped), LaggedStart(*[FadeIn(f, shift=LEFT * 0.1) for f in facts],
                                            lag_ratio=0.3), run_time=1.6)
        timing.hold_to_read(self, cap, facts, settle=0.8)
        _unsay(self, cap)
        self.play(*[FadeOut(m) for m in [live, body, apse, trail, rings, sun, zoom, facts, sped]])

        # ---- beat 2: the rate, term by term -----------------------------------
        eq = layout.explain_equation(
            self,
            [r"\dot\varpi", r"\approx", r"\tfrac{3}{2}\,n", r"\,\dfrac{J_2}{a^2}",
             r"\,\dfrac{1}{(1-e^2)^2}"],
            [
                (2, "n, the orbit's mean motion: distant orbits are slow"),
                (3, "J₂ = ½ Σ m a² of the giants' rings, weakened by distance a²"),
                (4, "eccentric orbits dip close to the rings, so they turn faster"),
            ],
            scale=1.0,
            where=UP * 0.4,
        )
        self.play(FadeOut(eq))

        # ---- beat 3: precession period vs distance, real ETNOs ----------------
        ax, ax_lbls = _plot_axes(
            [2.0, 3.0, 0.5], [1.0, 4.0, 1.0], 7.6, 3.9, np.array([0.3, 0.2, 0]),
            [(np.log10(v), f"{v:,}") for v in (100, 200, 500, 1000)],
            [(1, "10"), (2, "100"), (3, "1,000"), (4, "10,000")],
            "semi-major axis a (AU)", "one turn of the apse (Myr)")
        curve = widgets.curve(ax, np.log10(d["curve_a_au"]), np.log10(d["curve_period_myr"]),
                              color=P.ORANGE, stroke_width=3.5)
        c_lbl = layout.label(f"line: any orbit with perihelion {d['curve_q_au']:.0f} AU", font_size=18,
                             color=P.ORANGE)
        c_lbl.move_to(ax.c2p(2.98, 1.35), aligned_edge=RIGHT)
        age = DashedLine(ax.c2p(2.0, np.log10(d["age_myr"])), ax.c2p(3.0, np.log10(d["age_myr"])),
                         color=P.MUTED, stroke_width=2)
        age_lbl = layout.label("age of the Solar System", font_size=16, color=P.FG)
        age_lbl.next_to(ax.c2p(2.02, np.log10(d["age_myr"])), UP, buff=0.08, aligned_edge=LEFT)
        pts = VGroup(*[Dot(ax.c2p(np.log10(o["a"]), np.log10(o["period_myr"])), radius=0.08,
                           color=P.GREEN) for o in d["objects"]])
        pts_lbl = layout.label("dots: the ten real distant objects", font_size=18, color=P.GREEN)
        pts_lbl.move_to(ax.c2p(2.98, 1.7), aligned_edge=RIGHT)
        cap = _say(self, "One turn takes a few hundred Myr to a few Gyr, and every orbit has its own rate",
                   hold=False)
        self.play(Create(ax), FadeIn(ax_lbls))
        self.play(Create(curve), FadeIn(c_lbl), Create(age), FadeIn(age_lbl), run_time=1.6)
        self.play(LaggedStart(*[GrowFromCenter(p) for p in pts], lag_ratio=0.1), FadeIn(pts_lbl))
        timing.hold_to_read(self, cap, settle=0.8)
        _unsay(self, cap)
        self.play(*[FadeOut(m) for m in [ax, ax_lbls, curve, c_lbl, age, age_lbl, pts, pts_lbl]])

        # ---- beat 4: start them aligned, watch the alignment dissolve ---------
        C2 = np.array([-3.5, -0.3, 0.0])
        Rd = 2.2
        rates = np.deg2rad([o["rate_deg_per_myr"] for o in d["objects"]])
        v_start = np.deg2rad(60.0)
        clock = ValueTracker(0.0)
        ring = Circle(radius=Rd, color=P.MUTED, stroke_width=1).move_to(C2).set_stroke(opacity=0.5)
        sun2 = orbits.sun(radius=0.08).move_to(C2)
        arrows = always_redraw(lambda: VGroup(*[
            Arrow(C2, C2 + Rd * _u(v_start + r * clock.get_value()), buff=0, color=P.GREEN,
                  stroke_width=3, tip_length=0.15) for r in rates]))
        t_lbl = always_redraw(lambda: layout.label(
            f"t = {clock.get_value():,.0f} Myr", font_size=22, color=P.FG).move_to(C2 + np.array([0, 2.75, 0])))
        tt_ = np.array(d["decohere_t_myr"])
        rb_ = np.array(d["decohere_r_bar"])
        ax2, lbl2 = _plot_axes(
            [0, 4500, 1500], [0, 1, 0.5], 5.4, 3.4, np.array([3.6, 0.1, 0]),
            [(v, f"{v:,}") for v in (0, 1500, 3000, 4500)],
            [(v, f"{v:.1f}") for v in (0, 0.5, 1.0)],
            "time (Myr)", "alignment R̄")
        trace = always_redraw(lambda: widgets.curve(
            ax2, np.append(tt_[tt_ <= clock.get_value()], clock.get_value()),
            np.append(rb_[tt_ <= clock.get_value()], np.interp(clock.get_value(), tt_, rb_)),
            color=P.GREEN, stroke_width=3))
        cap = _say(self, "Thought experiment: start the ten real orbits perfectly aligned...", hold=False)
        self.play(FadeIn(sun2), Create(ring), FadeIn(arrows), FadeIn(t_lbl), Create(ax2), FadeIn(lbl2))
        self.add(trace)
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)
        cap = _say(self, "...and let each one precess at its own real rate", hold=False)
        self.play(clock.animate.set_value(d["t_half_myr"]), run_time=3.5, rate_func=rate_functions.linear)
        half = DashedLine(ax2.c2p(d["t_half_myr"], 0), ax2.c2p(d["t_half_myr"], 1), color=P.ORANGE,
                          stroke_width=2)
        half_lbl = layout.label(f"R̄ below 0.5 after {d['t_half_myr']:.0f} Myr", font_size=18,
                                color=P.ORANGE)
        half_lbl.next_to(ax2.c2p(d["t_half_myr"], 1.0), RIGHT, buff=0.12)
        self.play(Create(half), FadeIn(half_lbl))
        _unsay(self, cap)
        obs = DashedLine(ax2.c2p(0, d["r_bar_observed"]), ax2.c2p(4500, d["r_bar_observed"]),
                         color=P.GREEN, stroke_width=1.5).set_stroke(opacity=0.7)
        obs_lbl = layout.label(f"today's sky: {d['r_bar_observed']:.2f}", font_size=16, color=P.GREEN)
        obs_lbl.next_to(ax2.c2p(4500, d["r_bar_observed"]), UP, buff=0.08, aligned_edge=RIGHT)
        cap = _say(self, "Run on to today: the directions just keep churning", hold=False)
        self.play(Create(obs), FadeIn(obs_lbl), run_time=0.6)
        self.play(clock.animate.set_value(d["age_myr"]), run_time=5.5, rate_func=rate_functions.linear)
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)
        cap = _say(self, f"Chance brings back today's alignment only "
                         f"{100 * d['frac_time_above_observed']:.1f}% of the time")
        _unsay(self, cap)
        layout.show_takeaway(
            self, "Left alone, precession scrambles any alignment: something must hold it in place.")


class P07Resonance(Scene):
    """Learning goal: when two orbital periods form a simple ratio, the bodies
    meet at the same places over and over, so the tugs add up instead of
    cancelling; the same lock can protect an orbit (Neptune and Pluto)."""

    def construct(self):
        d = _data()["resonance"]
        _header(self, "PREFACE 07", "Resonance: tugs that add up")

        # ---- beat 1: where do the conjunctions land? --------------------------
        C = np.array([-3.4, -0.35, 0.0])
        r_in = 1.2
        ratio = ValueTracker(1.618)
        clock = ValueTracker(0.0)  # in inner orbital periods

        def r_out():
            return r_in * ratio.get_value() ** (2 / 3)

        n_conj = ValueTracker(60)

        def conj_angles(t_max=None):
            n_max = int(n_conj.get_value())
            r = ratio.get_value()
            out = []
            for kk in range(1, n_max + 1):
                t_k = kk * r / (r - 1.0)
                if t_max is not None and t_k > t_max:
                    break
                out.append(2 * np.pi * kk / (r - 1.0))
            return out

        inner = always_redraw(lambda: Circle(radius=r_in, color=P.FG, stroke_width=1.6)
                              .move_to(C).set_stroke(opacity=0.7))
        outer = always_redraw(lambda: Circle(radius=r_out(), color=P.GREEN, stroke_width=1.6)
                              .move_to(C).set_stroke(opacity=0.7))
        sun = orbits.sun(radius=0.08).move_to(C)
        p_in = always_redraw(lambda: Dot(C + r_in * _u(2 * np.pi * clock.get_value()), radius=0.09,
                                         color=P.FG))
        p_out = always_redraw(lambda: Dot(C + r_out() * _u(2 * np.pi * clock.get_value()
                                                         / ratio.get_value()),
                                          radius=0.09, color=P.GREEN))
        kicks = always_redraw(lambda: VGroup(*[
            Line(C + (r_out() + 0.08) * _u(th), C + (r_out() + 0.42) * _u(th), color=P.ORANGE,
                 stroke_width=4)
            for th in conj_angles(clock.get_value())]))
        lbl_in = layout.label("planet", font_size=16, color=P.FG).move_to(C + np.array([0, -0.35, 0]))
        lbl_out = layout.label("small body", font_size=16, color=P.GREEN)
        lbl_out.move_to(C + np.array([0, -2.55, 0]))

        # net push: kicks head to tail (as in Preface 05)
        S = np.array([2.2, -0.6, 0.0])
        Lk = 0.34

        def net_chain():
            angs = conj_angles(clock.get_value())
            if not angs:
                return VGroup()
            ch, _, _ = _chain(angs, Lk, S, P.ORANGE, 2.5)
            return ch

        net = always_redraw(net_chain)

        def net_text():
            angs = conj_angles(clock.get_value())
            tot = float(np.hypot(np.sum(np.cos(angs)), np.sum(np.sin(angs)))) if angs else 0.0
            return layout.label(f"{len(angs)} kicks, net push = {tot:.1f} kicks", font_size=18,
                                color=P.ORANGE).move_to(np.array([3.9, -1.5, 0]))

        net_lbl = always_redraw(net_text)
        net_head = layout.label("the kicks, head to tail", font_size=18, color=P.ORANGE)
        net_head.move_to(np.array([3.9, 0.35, 0]))
        rat_lbl = always_redraw(lambda: VGroup(
            MathTex(r"P_{\rm body} / P_{\rm planet} =", color=P.FG).scale(0.7),
            DecimalNumber(ratio.get_value(), num_decimal_places=3, color=P.FG).scale(0.7),
        ).arrange(RIGHT, buff=0.15).move_to(C + np.array([0, 2.75, 0])))

        self.play(FadeIn(sun), Create(inner), Create(outer), FadeIn(p_in), FadeIn(p_out),
                  FadeIn(lbl_in), FadeIn(lbl_out), FadeIn(rat_lbl), FadeIn(net_head))
        self.add(kicks, net, net_lbl)
        cap = _say(self, "Each time the two line up, the planet tugs the small body (orange tick)",
                   hold=False)
        self.play(clock.animate.set_value(10 * 1.618 / 0.618 + 0.3), run_time=7.5,
                  rate_func=rate_functions.linear)
        _unsay(self, cap)
        cap = _say(self, "No simple ratio: the meetings drift around, and the kicks cancel out")
        _unsay(self, cap)

        self.play(clock.animate.set_value(0.0), ratio.animate.set_value(1.5), run_time=1.2)
        cap = _say(self, "Now make the periods exactly 3 : 2 ...", hold=False)
        self.play(clock.animate.set_value(10 * 3 + 0.3), run_time=6.5, rate_func=rate_functions.linear)
        _unsay(self, cap)
        cap = _say(self, "...every meeting is at the same place, so the kicks pile up")
        _unsay(self, cap)

        # ---- beat 2: turn the knob --------------------------------------------
        self.play(FadeOut(p_in), FadeOut(p_out), FadeOut(net), FadeOut(net_head), FadeOut(net_lbl),
                  run_time=0.4)
        self.remove(p_in, p_out, net, net_lbl)
        n_conj.set_value(24)
        self.play(clock.animate.set_value(200.0), run_time=0.01)
        cap = _say(self, "Turn the knob: 24 meetings pile up in one place only very near 3 : 2",
                   hold=False)
        self.play(ratio.animate.set_value(1.46), run_time=2.0, rate_func=rate_functions.smooth)
        self.play(ratio.animate.set_value(1.54), run_time=7.0, rate_func=rate_functions.linear)
        self.play(ratio.animate.set_value(1.5), run_time=2.5, rate_func=rate_functions.smooth)
        timing.hold_to_read(self, cap, settle=0.3)
        _unsay(self, cap)
        self.play(*[FadeOut(m) for m in [inner, outer, kicks, sun, lbl_in, lbl_out, rat_lbl]])
        self.remove(kicks, inner, outer, rat_lbl)

        # ---- beat 3: Neptune and Pluto in Neptune's turning frame -------------
        pr = d["protected"]
        un = d["unprotected"]
        Cn = np.array([-3.0, -0.2, 0.0])
        sc = 2.35 / d["pluto"]["Q"]
        nep_orbit = DashedVMobject(Circle(radius=sc * d["neptune"]["a"], color=P.FG, stroke_width=1.5)
                                   .move_to(Cn), num_dashes=60).set_stroke(opacity=0.6)
        sun3 = orbits.sun(radius=0.08).move_to(Cn)
        nep = Dot(Cn + sc * d["neptune"]["a"] * RIGHT, radius=0.11, color=P.FG)
        nep_lbl = layout.label("Neptune", font_size=18, color=P.FG).next_to(nep, DL, buff=0.04)
        frame_lbl = layout.label("view turning with Neptune", font_size=18, color=P.MUTED)
        frame_lbl.move_to(Cn + np.array([0, 2.95, 0]))

        def path(xs, ys, color, width=2.6):
            m = VMobject(color=color, stroke_width=width)
            m.set_points_as_corners([Cn + sc * np.array([x, y, 0]) for x, y in zip(xs, ys)])
            return m

        pl_path = path(pr["x"], pr["y"], P.GREEN)
        idx = ValueTracker(0)
        npt = len(pr["x"])
        pl_dot = always_redraw(lambda: Dot(
            Cn + sc * np.array([pr["x"][int(idx.get_value())], pr["y"][int(idx.get_value())], 0]),
            radius=0.08, color=P.GREEN))
        facts = VGroup(
            layout.label("Pluto", font_size=20, color=P.GREEN, weight="BOLD"),
            layout.label(f"perihelion {d['pluto']['q']:.1f} AU: inside Neptune's "
                         f"{d['neptune']['a']:.1f} AU orbit", font_size=17, color=P.FG),
            layout.label(f"periods {d['pluto']['period_yr']:.0f} yr / {d['neptune']['period_yr']:.0f} yr"
                         f" = {d['period_ratio']:.3f}  ≈ 3 : 2", font_size=17, color=P.FG),
        ).arrange(DOWN, aligned_edge=LEFT, buff=0.2).move_to(np.array([3.55, 1.6, 0]))
        cap = _say(self, "A real case: Neptune and Pluto, watched from a view that turns with Neptune",
                   hold=False)
        self.play(FadeIn(sun3), Create(nep_orbit), FadeIn(nep), FadeIn(nep_lbl), FadeIn(frame_lbl))
        self.play(LaggedStart(*[FadeIn(f) for f in facts], lag_ratio=0.3), run_time=1.2)
        timing.hold_to_read(self, cap, settle=0.2)
        _unsay(self, cap)
        cycle = 3 * d["neptune"]["period_yr"]
        cap = _say(self, f"Over one 3 : 2 cycle ({cycle:.0f} yr) Pluto traces this loop", hold=False)
        self.add(pl_dot)
        self.play(Create(pl_path), idx.animate.set_value(npt - 1), run_time=6.0,
                  rate_func=rate_functions.linear)
        _unsay(self, cap)
        near = layout.label(f"closest to Neptune: {pr['min_dist_au']:.0f} AU", font_size=18,
                            color=P.GREEN)
        near2 = layout.label(f"at each meeting Pluto is out near {pr['r_at_conjunction_au']:.0f} AU",
                             font_size=17, color=P.FG)
        VGroup(near, near2).arrange(DOWN, aligned_edge=LEFT, buff=0.18).move_to(np.array([3.55, -0.2, 0]))
        cap = _say(self, "It crosses Neptune's orbit, yet the lock keeps every meeting far away",
                   hold=False)
        self.play(FadeIn(near), FadeIn(near2))
        timing.hold_to_read(self, cap, near, near2, settle=0.4)
        _unsay(self, cap)

        bad = path(un["x"], un["y"], P.RED, 2.0)
        bad_lbl = layout.label(f"same orbit, opposite phase: {un['min_dist_au']:.1f} AU",
                               font_size=18, color=P.RED)
        bad_lbl.move_to(np.array([3.55, -1.6, 0]))
        cap = _say(self, "Shift Pluto half a cycle and it would graze Neptune and be thrown out",
                   hold=False)
        self.play(FadeOut(pl_dot), pl_path.animate.set_stroke(opacity=0.35), Create(bad), FadeIn(bad_lbl),
                  run_time=2.5)
        timing.hold_to_read(self, cap, bad_lbl, settle=0.4)
        _unsay(self, cap)
        cap = _say(self, "Planet Nine's resonances can shelter distant bodies the same way")
        _unsay(self, cap)
        layout.show_takeaway(
            self, "Simple period ratios make tugs repeat: they can pump an orbit, or protect it.")
