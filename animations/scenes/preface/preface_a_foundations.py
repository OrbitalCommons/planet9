"""Preface A -- Foundations: the hook, Newton's ellipse, and Kepler's laws.

Scenes: P00Hook, P01Ellipse, P02KeplerLaws.

Every orbit, speed and period on screen comes from
``anim.json -> preface -> a_foundations`` (crates/p9-anim-data/src/preface/a_foundations.rs).
"""
import numpy as np
from manim import (
    Arrow,
    Axes,
    Circle,
    Create,
    DashedLine,
    DashedVMobject,
    Dot,
    DOWN,
    FadeIn,
    FadeOut,
    GrowArrow,
    LEFT,
    Line,
    Polygon,
    RIGHT,
    Scene,
    Sector,
    Transform,
    UP,
    UR,
    ValueTracker,
    VGroup,
    VMobject,
    Write,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, timing


def _data():
    return dataio.section("preface")["a_foundations"]


def _pts(xy, scale, origin):
    """Scene points from AU coordinates (x, y[, z]) about ``origin``."""
    o = np.asarray(origin, dtype=float)
    return [o + np.array([scale * p[0], scale * p[1], 0.0]) for p in xy]


def _poly(points, color, width=2.0, opacity=1.0):
    m = VMobject(stroke_color=color, stroke_width=width, stroke_opacity=opacity)
    m.set_points_as_corners(points)
    return m


def _panel(lines, anchor, font_size=20, buff=0.2):
    """A left-aligned stack of (text, colour) rows with its top-left at ``anchor``."""
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in lines])
    g.arrange(DOWN, buff=buff, aligned_edge=LEFT)
    g.move_to(anchor, aligned_edge=UP + LEFT)
    return g


def _pointers(varpis_deg, origin, color, length=1.6, width=3.5):
    """Slim arrows from the Sun toward each longitude of perihelion."""
    return VGroup(*[
        Arrow(origin, origin + length * np.array([np.cos(np.deg2rad(v)), np.sin(np.deg2rad(v)), 0.0]),
              buff=0, color=color, stroke_width=width, tip_length=0.16,
              max_tip_length_to_length_ratio=0.12)
        for v in varpis_deg
    ])


def _swap_caption(scene, old, text):
    new = layout.caption(text)
    if old is not None:
        scene.play(FadeOut(old), run_time=0.3)
    scene.play(FadeIn(new, shift=UP * 0.1), run_time=0.5)
    return new


class P00Hook(Scene):
    """Goal: the puzzle in one image -- real distant orbits all point one way,
    and a hidden planet on the opposite side could be why."""

    def construct(self):
        d = _data()["hook"]
        title = layout.title_card("Planet Nine", "the case for a world we have never seen")
        self.play(Write(title[0]), run_time=1.2)
        self.play(FadeIn(title[1], shift=UP * 0.2))
        self.wait(1.0)
        self.play(FadeOut(title[1]), title[0].animate.scale(0.6).to_edge(UP, buff=0.3))

        # true-scale map, seen from above the planets' plane
        etnos = d["etnos"]
        allx = [p[0] for o in etnos for p in o["xyz"]] + [p[0] for p in d["p9"]["xyz"]]
        ally = [p[1] for o in etnos for p in o["xyz"]] + [p[1] for p in d["p9"]["xyz"]]
        s = 5.9 / (max(ally) - min(ally))
        origin = np.array([-1.6 - s * (max(allx) + min(allx)) / 2,
                           0.1 - s * (max(ally) + min(ally)) / 2, 0.0])

        # start close in on Neptune's orbit, then pull back to the true scale
        big = 2.1 / (d["neptune_a"] * s)
        sun = orbits.sun(radius=0.06).move_to(origin)
        nep = Circle(radius=d["neptune_a"] * s, color=P.FG, stroke_width=1.8).move_to(origin)
        nep.scale(big, about_point=origin)
        nep_lbl = layout.label("Neptune's orbit: 30 AU, the edge of the known planets", font_size=22)
        nep_lbl.next_to(nep, DOWN, buff=0.2)
        self.play(Create(nep), FadeIn(sun), FadeIn(nep_lbl))
        self.wait(1.2)
        cap = _swap_caption(self, None, "Now zoom out, keeping everything to scale.")
        self.play(FadeOut(nep_lbl), nep.animate.scale(1 / big, about_point=origin), run_time=2.2)
        timing.hold_to_read(self, cap, settle=0.3)

        swarm = VGroup(*[_poly(_pts(o["xyz"], s, origin), P.GREEN, width=1.8, opacity=0.85)
                         for o in etnos])
        cap = _swap_caption(self, cap, f"Far beyond it: {len(etnos)} real orbits, drawn to scale (Brown 2017).")
        self.play(Create(swarm, lag_ratio=0.15), run_time=3.0)
        timing.hold_to_read(self, cap, settle=0.4)

        # each orbit's closest-approach direction
        arrows = _pointers([o["varpi_deg"] for o in etnos], origin, P.GREEN)
        cap = _swap_caption(self, cap, "Arrows: where each orbit comes closest to the Sun.")
        self.play(swarm.animate.set_stroke(opacity=0.3), GrowArrow(arrows[0]),
                  *[GrowArrow(a) for a in arrows[1:]], run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.4)
        cap = _swap_caption(self, cap, "They bunch up on one side. Random orbits would point every which way.")
        timing.hold_to_read(self, cap, settle=0.8)

        p9 = d["p9"]
        p9_orbit = DashedVMobject(_poly(_pts(p9["xyz"], s, origin), P.BLUE, width=3.0),
                                  num_dashes=70)
        p9_arrow = _pointers([p9["varpi_deg"]], origin, P.BLUE, width=5)
        cap = _swap_caption(self, cap, "A hidden planet, pointing the opposite way, could hold them there.")
        self.play(Create(p9_orbit), run_time=2.2)
        self.play(GrowArrow(p9_arrow[0]))

        legend = _panel([
            ("Neptune's orbit (30 AU)", P.FG),
            (f"{len(etnos)} distant objects", P.GREEN),
            (f"Planet Nine? ~{p9['mass_earth']:.0f} Earth masses", P.BLUE),
            (f"a = {p9['a']:.0f} AU, never seen", P.BLUE),
        ], anchor=np.array([2.6, 1.6, 0.0]), font_size=20)
        self.play(FadeIn(legend, shift=LEFT * 0.1))
        timing.hold_to_read(self, cap, legend, settle=1.0)

        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Distant orbits line up. Is an unseen planet herding them?")


class P01Ellipse(Scene):
    """Goal: an inverse-square pull turns every bound launch into an ellipse
    with the Sun at one focus; speed sets the size, and a, e name the shape."""

    def construct(self):
        d = _data()
        L = d["launch"]
        self.add(layout.concept_badge("PREFACE 01"))

        law = layout.equation_card(r"F = \frac{G M_\odot\, m}{r^2}").scale(0.9)
        law.move_to(np.array([-4.7, 2.3, 0.0]))
        self.play(FadeIn(law[0]), Write(law[1]))
        cap = _swap_caption(self, None, "Newton: one pull toward the Sun, weakening as 1/r squared.")
        timing.hold_to_read(self, cap, settle=0.2)

        s = 0.056
        S = np.array([1.9, -0.1, 0.0])
        sun = orbits.sun(radius=0.12).move_to(S)
        r0 = L["r0_au"]
        launch_pt = S + RIGHT * r0 * s
        pad = Dot(launch_pt, radius=0.05, color=P.FG)
        pad_lbl = layout.label("start: 30 AU out,\nNeptune's distance", font_size=16, color=P.MUTED)
        pad_lbl.next_to(pad, DOWN + RIGHT, buff=0.1)
        self.play(FadeIn(sun), FadeIn(pad), FadeIn(pad_lbl))
        cap = _swap_caption(self, cap, "Throw a body sideways from Neptune's distance. Vary only the speed.")

        corner = np.array([-6.6, 1.35, 0.0])
        panel = _panel([("launch speed", P.MUTED), ("", P.FG), ("", P.FG)], anchor=corner)
        self.add(panel[0])
        done = VGroup()
        arrow = None
        for path in L["paths"]:
            pts = _pts(path["xy"], s, S)
            v = path["speed_kms"]
            new_arrow = Arrow(launch_pt, launch_pt + UP * 0.26 * v, buff=0, color=P.ORANGE,
                              stroke_width=5, max_tip_length_to_length_ratio=0.12)
            if path["bound"]:
                shape = "circle" if path["e"] < 0.01 else f"ellipse, e = {path['e']:.2f}"
                reach = f"closest {path['q']:.0f} AU, farthest {path['big_q']:.0f} AU"
            else:
                shape = "escapes: never returns"
                reach = f"escape speed = {L['v_esc_kms']:.2f} km/s"
            rows = _panel([("launch speed", P.MUTED), (f"{v:.2f} km/s", P.ORANGE),
                           (shape, P.FG), (reach, P.FG)], anchor=corner)
            anims = [Transform(panel, rows)]
            anims.append(GrowArrow(new_arrow) if arrow is None else Transform(arrow, new_arrow))
            self.play(*anims, run_time=0.6)
            if arrow is None:
                arrow = new_arrow
            k = ValueTracker(1)
            trail = always_redraw(lambda: _poly(pts[: max(2, int(k.get_value()))], P.FG, width=3))
            body = always_redraw(lambda: Dot(pts[min(len(pts) - 1, int(k.get_value()))],
                                             radius=0.07, color=P.FG))
            self.add(trail, body)
            self.play(k.animate.set_value(len(pts) - 1), run_time=2.6 if path["bound"] else 1.8,
                      rate_func=rate_functions.linear)
            trail.clear_updaters()
            self.remove(body)
            self.play(trail.animate.set_stroke(color=P.MUTED, width=1.8), run_time=0.3)
            done.add(trail)
            if path["bound"] and path["e"] < 0.01:
                cap = _swap_caption(self, cap, f"At {v:.1f} km/s: a circle. That is just how fast Neptune moves.")
                timing.hold_to_read(self, cap, settle=0.3)
                cap = _swap_caption(self, cap, "Faster: it swings farther out. Fast enough, and it never returns.")
        cap = _swap_caption(self, cap, "Every bound path is an ellipse, with the Sun at one focus.")
        timing.hold_to_read(self, cap, settle=0.5)
        self.play(*[FadeOut(m) for m in (law, sun, pad, pad_lbl, panel, arrow, done, cap)])

        vv = layout.explain_equation(
            self,
            [r"v^2", "=", r"G M_\odot", r"\left(\frac{2}{r}-\frac{1}{a}\right)"],
            [
                (3, "at fixed r: more speed means a bigger orbit, a"),
                (-1, "at v squared = 2GM/r, a runs off to infinity: escape"),
            ],
            where=UP * 0.8,
        )
        self.play(FadeOut(vv))
        self._anatomy(d)

    def _anatomy(self, d):
        sed = d["sedna"]
        a_au = sed["a"]
        A = 2.6
        s = A / a_au
        C = np.array([-2.4, 0.1, 0.0])
        e = ValueTracker(0.0)

        def geom():
            ev = e.get_value()
            c = A * ev
            return ev, c, C + RIGHT * c

        def ell():
            ev, _, focus = geom()
            return orbits.ellipse_orbit(A, ev, color=P.FG).shift(focus)

        orbit = always_redraw(ell)
        sun = always_redraw(lambda: orbits.sun(radius=0.11).move_to(geom()[2]))
        empty = always_redraw(lambda: Circle(radius=0.07, color=P.MUTED, stroke_width=2)
                              .move_to(C - RIGHT * A * e.get_value()))
        axis = DashedLine(C + LEFT * A, C + RIGHT * A, color=P.MUTED, stroke_width=1.2)
        q_seg = always_redraw(lambda: Line(geom()[2], C + RIGHT * A, color=P.GREEN, stroke_width=6))
        Q_seg = always_redraw(lambda: Line(C + LEFT * A, geom()[2], color=P.ORANGE, stroke_width=6))
        q_lbl = always_redraw(lambda: layout.label("q", font_size=22, color=P.GREEN)
                              .next_to(C + RIGHT * A, RIGHT, buff=0.15))
        Q_lbl = always_redraw(lambda: layout.label("Q", font_size=22, color=P.ORANGE)
                              .next_to((geom()[2] + C + LEFT * A) / 2, DOWN, buff=0.12))

        def rows():
            ev = e.get_value()
            return _panel([
                (f"a = {a_au:.0f} AU  (size, held fixed)", P.TEAL),
                (f"e = {ev:.2f}  (stretch)", P.FG),
                (f"q = a(1 - e) = {a_au * (1 - ev):.0f} AU", P.GREEN),
                (f"Q = a(1 + e) = {a_au * (1 + ev):.0f} AU", P.ORANGE),
            ], anchor=np.array([1.3, 2.4, 0.0]), font_size=22, buff=0.28)

        panel = always_redraw(rows)
        cap = _swap_caption(self, None, "Keep the size a fixed and turn up the stretch e.")
        self.play(Create(orbit), FadeIn(sun), Create(axis))
        self.add(empty, q_seg, Q_seg, q_lbl, Q_lbl)
        self.play(FadeIn(panel))
        self.play(e.animate.set_value(0.5), run_time=3.0, rate_func=rate_functions.smooth)
        cap = _swap_caption(self, cap, "The Sun sits at a focus: q is the closest approach, Q the farthest.")
        timing.hold_to_read(self, cap, settle=0.3)
        self.play(e.animate.set_value(sed["e"]), run_time=2.5, rate_func=rate_functions.smooth)

        for m in (orbit, sun, empty, q_seg, Q_seg, q_lbl, Q_lbl, panel):
            m.clear_updaters()
        name = layout.label(f"This is Sedna's real shape: e = {sed['e']:.2f}", font_size=22,
                            color=P.GREEN, weight="BOLD")
        name.next_to(panel, DOWN, buff=0.45, aligned_edge=LEFT)
        nep_r = d["hook"]["neptune_a"] * s
        nep = Circle(radius=nep_r, color=P.FG, stroke_width=1.6).move_to(sun.get_center())
        nep_lbl = layout.label("Neptune's orbit", font_size=16, color=P.FG)
        nep_lbl.move_to(sun.get_center() + np.array([1.1, -1.0, 0.0]))
        nep_ptr = Line(nep_lbl.get_corner(UP + LEFT), nep.point_at_angle(-np.pi / 4), color=P.FG,
                       stroke_width=1.2, buff=0.05)
        self.play(orbit.animate.set_color(P.GREEN), FadeIn(name), Create(nep),
                  FadeIn(nep_lbl), Create(nep_ptr))
        ratio = sed["q"] / d["hook"]["neptune_a"]
        cap = _swap_caption(self, cap, f"Even at its closest, Sedna stays {ratio:.1f} times farther out than Neptune.")
        timing.hold_to_read(self, cap, settle=1.0)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Gravity makes ellipses: a sets the size, e the stretch, the Sun a focus.")


class P02KeplerLaws(Scene):
    """Goal: equal areas in equal times (fast when near, slow when far), and
    P^2 = a^3 (far orbits take millennia) -- with real bodies."""

    def construct(self):
        d = _data()
        self.add(layout.concept_badge("PREFACE 02"))
        self._equal_areas()
        self._sedna_clock(d["sedna"])
        self._third_law(d["kepler3"])

    def _equal_areas(self):
        a, e, n = 2.9, 0.7, 10
        C = np.array([-2.6, 0.0, 0.0])
        S = C + RIGHT * a * e
        orbit = orbits.ellipse_orbit(a, e, color=P.FG).shift(S)
        sun = orbits.sun(radius=0.12).move_to(S)
        self.play(Create(orbit), FadeIn(sun))

        def pt(M):
            return S + orbits.orbit_point(a, e, orbits.nu_from_mean_anomaly(e, M))

        def sector(M0, M1, color):
            steps = max(3, int(40 * (M1 - M0) / (2 * np.pi / n)))
            poly = [S] + [pt(m) for m in np.linspace(M0, M1, steps)]
            return Polygon(*poly, stroke_width=0.8, color=color).set_fill(color, opacity=0.35)

        colors = [P.GREEN, P.TEAL]
        M = ValueTracker(0.0)
        cap = _swap_caption(self, None, "Kepler's second law: in equal times, the Sun-body line sweeps equal areas.")
        panel = _panel([
            ("each wedge: 1/10 of the orbit's time", P.FG),
            ("and exactly 1/10 of its area", P.FG),
        ], anchor=np.array([1.4, 2.6, 0.0]), font_size=21)
        self.play(FadeIn(panel))

        def live():
            m = M.get_value()
            k = min(int(m / (2 * np.pi / n)), n - 1)
            m0 = k * 2 * np.pi / n
            if m - m0 < 1e-3:
                return VGroup()
            return sector(m0, m, colors[k % 2])

        wedge = always_redraw(live)
        body = always_redraw(lambda: Dot(pt(M.get_value()), radius=0.08, color=P.FG).set_z_index(4))
        spoke = always_redraw(lambda: Line(S, pt(M.get_value()), color=P.FG, stroke_width=1.5))
        self.add(wedge, spoke, body)
        done = VGroup()
        for k in range(n):
            self.play(M.animate.set_value((k + 1) * 2 * np.pi / n), run_time=1.1,
                      rate_func=rate_functions.linear)
            w = sector(k * 2 * np.pi / n, (k + 1) * 2 * np.pi / n, colors[k % 2])
            done.add(w)
            self.add(w)
        self.remove(wedge)
        ratio = (1 + e) / (1 - e)
        more = _panel([
            ("near the Sun: thin wedges, moving fast", P.FG),
            ("far out: fat wedges, crawling", P.FG),
            (f"speed at perihelion: {ratio:.1f}x aphelion", P.ORANGE),
        ], anchor=panel.get_corner(DOWN + LEFT) + DOWN * 0.45, font_size=21)
        self.play(FadeIn(more))
        timing.hold_to_read(self, more, settle=1.0)
        self.play(*[FadeOut(m) for m in (orbit, sun, done, spoke, body, panel, more, cap)])

    def _sedna_clock(self, sed):
        kyr = [t / 1000.0 for t in sed["t_yr"]]
        ax = Axes(x_range=[0, 12, 2], y_range=[0, 1000, 200], x_length=7.4, y_length=4.2,
                  axis_config={"color": P.MUTED, "include_tip": False})
        ax.move_to(np.array([-2.3, 0.35, 0.0]))
        for t in range(2, 13, 2):
            ax.add(layout.label(str(t), font_size=17).next_to(ax.c2p(t, 0), DOWN, buff=0.12))
        for r in range(200, 1001, 200):
            ax.add(layout.label(f"{r:,}", font_size=17).next_to(ax.c2p(0, r), LEFT, buff=0.12))
        xl = layout.label("time since perihelion (thousand years)", font_size=18, color=P.FG)
        xl.next_to(ax, DOWN, buff=0.15)
        yl = layout.label("distance from Sun (AU)", font_size=18, color=P.FG).rotate(np.pi / 2)
        yl.next_to(ax, LEFT, buff=0.15)
        band = Polygon(ax.c2p(0, 0), ax.c2p(12, 0), ax.c2p(12, 100), ax.c2p(0, 100),
                       stroke_width=0).set_fill(P.GREEN, opacity=0.22)
        band_lbl = layout.label("closer than 100 AU", font_size=16, color=P.GREEN)
        band_lbl.next_to(ax.c2p(6, 100), UP, buff=0.08)
        curve = VMobject(stroke_color=P.GREEN, stroke_width=3.5)
        curve.set_points_as_corners([ax.c2p(t, r) for t, r in zip(kyr, sed["r_au"])])
        cap = _swap_caption(self, None, "The same law, for a real body: Sedna's distance through one orbit.")
        self.play(Create(ax), FadeIn(xl), FadeIn(yl))
        self.play(FadeIn(band), FadeIn(band_lbl))
        self.play(Create(curve), run_time=4.0, rate_func=rate_functions.linear)
        pct_out = 100 * (1 - sed["frac_inside_100"])
        panel = _panel([
            ("Sedna", P.GREEN),
            (f"one orbit: {sed['period_yr']:,.0f} yr", P.FG),
            (f"inside 100 AU: {sed['yr_inside_100']:.0f} yr", P.FG),
            (f"beyond 100 AU: {pct_out:.0f}% of the time", P.FG),
        ], anchor=np.array([2.6, 2.3, 0.0]), font_size=21)
        self.play(FadeIn(panel))
        cap = _swap_caption(self, cap, "It races past perihelion, then crawls far out, too faint to see.")
        timing.hold_to_read(self, cap, panel, settle=1.0)
        self.play(*[FadeOut(m) for m in (ax, xl, yl, band, band_lbl, curve, panel, cap)])

    def _third_law(self, k3):
        eq = layout.explain_equation(
            self,
            [r"P^2", "=", r"a^3"],
            [
                (0, "orbital period, in years, squared"),
                (2, "semi-major axis, in AU, cubed"),
                (-1, "so the period grows as a to the 3/2: ten times farther, 32 times slower"),
            ],
            scale=1.3,
            where=UP * 0.8,
        )
        self.play(FadeOut(eq))

        ax = Axes(x_range=[0, 3, 1], y_range=[0, 4.5, 1], x_length=6.6, y_length=4.4,
                  axis_config={"color": P.MUTED, "include_tip": False})
        ax.move_to(np.array([-2.4, 0.45, 0.0]))
        ticks = VGroup()
        for kx, txt in enumerate(["1", "10", "100", "1000"]):
            ticks.add(layout.label(txt, font_size=16, color=P.FG).next_to(ax.c2p(kx, 0), DOWN, buff=0.12))
        for ky, txt in enumerate(["1", "10", "100", "1,000", "10,000"]):
            ticks.add(layout.label(txt, font_size=16, color=P.FG).next_to(ax.c2p(0, ky), LEFT, buff=0.12))
        xl = layout.label("semi-major axis a (AU)", font_size=18).next_to(ax, DOWN, buff=0.45)
        yl = layout.label("period (years)", font_size=18).rotate(np.pi / 2).next_to(ax, LEFT, buff=0.75)
        line = Line(ax.c2p(0, 0), ax.c2p(3, 4.5), color=P.MUTED, stroke_width=2)
        self.play(Create(ax), FadeIn(ticks), FadeIn(xl), FadeIn(yl))
        self.play(Create(line))

        dots = VGroup()
        labels = VGroup()
        for j, b in enumerate(k3["bodies"]):
            p = ax.c2p(np.log10(b["a"]), np.log10(b["period_yr"]))
            dots.add(Dot(p, radius=0.07, color=P.GREEN))
            lbl = layout.label(b["name"], font_size=16, color=P.GREEN)
            if j == 0:
                lbl.next_to(p, RIGHT + UP * 0.6, buff=0.12)
            elif j % 2 == 1:
                lbl.next_to(p, UP + LEFT, buff=0.06)
            else:
                lbl.next_to(p, DOWN + RIGHT, buff=0.06)
            labels.add(lbl)
        cap = _swap_caption(self, None, "Real bodies, from Earth out to Pluto, all sit on one line.")
        self.play(FadeIn(dots, lag_ratio=0.2), FadeIn(labels, lag_ratio=0.2), run_time=2.0)
        timing.hold_to_read(self, cap, settle=0.4)

        p9 = k3["p9"]
        p9_dots = VGroup(*[Dot(ax.c2p(np.log10(o["a"]), np.log10(o["period_yr"])), radius=0.08,
                               color=P.BLUE) for o in p9])
        lo = min(o["period_yr"] for o in p9)
        hi = max(o["period_yr"] for o in p9)
        p9_lbl = layout.label("Planet Nine estimates", font_size=16, color=P.BLUE)
        p9_lbl.next_to(p9_dots, UP + LEFT, buff=0.15)
        cap = _swap_caption(self, cap, f"Planet Nine's published orbits: one lap takes {lo:,.0f} to {hi:,.0f} years.")
        self.play(FadeIn(p9_dots), FadeIn(p9_lbl))
        timing.hold_to_read(self, cap, settle=0.4)

        # how far round its orbit Planet Nine has moved since Newton
        deg = k3["p9_deg_since_principia"]
        Wc = np.array([4.6, 0.2, 0.0])
        R = 1.2
        ring = Circle(radius=R, color=P.BLUE, stroke_width=2).move_to(Wc).set_stroke(opacity=0.6)
        wedge = Sector(radius=R, angle=np.deg2rad(deg), start_angle=np.pi / 2, arc_center=Wc,
                       fill_opacity=0.75, stroke_width=0, color=P.BLUE)
        head = layout.label(f"since Newton, {k3['years_since_principia']:.0f} yr ago:", font_size=18)
        head.next_to(ring, UP, buff=0.3)
        foot = _panel([
            (f"Planet Nine: {deg:.0f} degrees of one lap", P.BLUE),
            (f"Neptune: {k3['neptune_orbits_since_found']:.1f} laps since found", P.FG),
        ], anchor=ring.get_bottom() + DOWN * 0.3 + LEFT * 1.6, font_size=17)
        cap = _swap_caption(self, cap, "Since Newton's Principia, a Planet Nine would barely have moved.")
        self.play(Create(ring), FadeIn(head))
        self.play(FadeIn(wedge), FadeIn(foot))
        timing.hold_to_read(self, cap, foot, settle=1.0)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Far orbits are slow: one lap of Planet Nine takes ~10,000 years.")
