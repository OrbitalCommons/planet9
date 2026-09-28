"""Preface B -- Orienting an orbit, and the geography of the outer system.

Scenes: P03Elements, P04Geography.

Every orbit and number on screen comes from
``anim.json -> preface -> b_elements`` (crates/p9-anim-data/src/preface/b_elements.rs).
"""
import numpy as np
from manim import (
    Annulus,
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
    UP,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    rate_functions,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, timing


def _data():
    return dataio.section("preface")["b_elements"]


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


def _swap_caption(scene, old, text):
    new = layout.caption(text)
    if old is not None:
        scene.play(FadeOut(old), run_time=0.3)
    scene.play(FadeIn(new, shift=UP * 0.1), run_time=0.5)
    return new


# ---- a hand-rolled oblique view of 3-D orbit geometry -------------------------

AZ = np.deg2rad(-28.0)


def _rot_z(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _rot_x(t):
    c, s = np.cos(t), np.sin(t)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


class _View:
    """Projects ecliptic-frame 3-D points onto the screen for a camera at
    elevation ``elev`` above the ecliptic (90 deg = looking straight down)."""

    def __init__(self, origin, elev):
        self.o = np.asarray(origin, dtype=float)
        self.elev = elev

    def p(self, v):
        x, y, z = v
        xr = x * np.cos(AZ) - y * np.sin(AZ)
        yr = x * np.sin(AZ) + y * np.cos(AZ)
        return self.o + np.array([xr, yr * np.sin(self.elev) + z * np.cos(self.elev), 0.0])


def _frame(i, node, argp):
    """Rotation taking orbit-plane (perihelion along +x) coordinates to the ecliptic."""
    return _rot_z(node) @ _rot_x(i) @ _rot_z(argp)


def _orbit_xyz(a, e, i, node, argp, n=240):
    E = np.linspace(0.0, 2 * np.pi, n + 1)
    pf = np.stack([a * (np.cos(E) - e), a * np.sqrt(1 - e * e) * np.sin(E), 0.0 * E], axis=1)
    return pf @ _frame(i, node, argp).T


def _split_by_plane(view, xyz, color, width=3.0):
    """Orbit drawn solid above the ecliptic, faint below it."""
    above, below = VGroup(), VGroup()
    run, sign = [xyz[0]], xyz[0][2] >= -1e-9
    for v in xyz[1:]:
        s = v[2] >= -1e-9
        run.append(v)
        if s != sign:
            (above if sign else below).add(_poly([view.p(u) for u in run], color, width,
                                                 1.0 if sign else 0.3))
            run, sign = [v], s
    (above if sign else below).add(_poly([view.p(u) for u in run], color, width,
                                         1.0 if sign else 0.3))
    below.set_z_index(0)
    above.set_z_index(3)
    return VGroup(below, above)


def _arc_xyz(frame, r, t0, t1, n=60):
    t = np.linspace(t0, t1, n)
    return (np.stack([r * np.cos(t), r * np.sin(t), 0.0 * t], axis=1)) @ frame.T


class P03Elements(Scene):
    """Goal: six numbers pin an orbit in space -- a, e for its shape; i, the
    node and the argument of perihelion for its orientation; the longitude of
    perihelion is the combination that says which way it points."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 03"))
        d = _data()["demo"]
        R = 2.75                              # ecliptic disk radius (scene units)
        k = 3.3 / d["big_q"]                  # scene units per AU
        a, e = d["a"] * k, d["e"]
        I = ValueTracker(0.0)
        NODE = ValueTracker(0.0)
        ARGP = ValueTracker(0.0)
        ELEV = ValueTracker(np.deg2rad(30.0))
        origin = np.array([-2.0, -0.05, 0.0])

        def view():
            return _View(origin, ELEV.get_value())

        def disk():
            v = view()
            ring = [v.p((R * np.cos(t), R * np.sin(t), 0.0)) for t in np.linspace(0, 2 * np.pi, 120)]
            g = Polygon(*ring, stroke_color=P.MUTED, stroke_width=1.5)
            return g.set_fill(P.MUTED, opacity=0.16).set_z_index(1)

        def ref_arrow():
            v = view()
            return Arrow(v.p((0, 0, 0)), v.p((R + 0.35, 0, 0)), buff=0, color=P.FG,
                         stroke_width=2.5, tip_length=0.18).set_z_index(2)

        def ref_label():
            v = view()
            return layout.label("reference\ndirection", font_size=16, color=P.FG).next_to(
                v.p((R + 0.35, 0, 0)), RIGHT, buff=0.1)

        def frame():
            return _frame(I.get_value(), NODE.get_value(), ARGP.get_value())

        def orbit():
            xyz = _orbit_xyz(a, e, I.get_value(), NODE.get_value(), ARGP.get_value())
            return _split_by_plane(view(), xyz, P.GREEN)

        def peri():
            v = view()
            tip = frame() @ np.array([a * (1 - e), 0.0, 0.0])
            return VGroup(Line(v.p((0, 0, 0)), v.p(tip), color=P.ORANGE, stroke_width=3),
                          Dot(v.p(tip), radius=0.07, color=P.ORANGE)).set_z_index(4)

        sun = always_redraw(lambda: orbits.sun(radius=0.1).move_to(view().p((0, 0, 0))).set_z_index(5))
        plane = always_redraw(disk)
        ref = always_redraw(ref_arrow)
        ref_lbl = always_redraw(ref_label)
        orb = always_redraw(orbit)
        apse = always_redraw(peri)

        cap = _swap_caption(self, None, "Six numbers pin down any orbit. Here is a real one: 2012 VP113.")
        self.play(FadeIn(plane), FadeIn(sun), FadeIn(ref), FadeIn(ref_lbl))
        ecl_lbl = layout.label("ecliptic: the planets' plane", font_size=17, color=P.MUTED)
        ecl_lbl.move_to(view().p((-0.35 * R, -0.95 * R, 0.0)) + DOWN * 0.3)
        self.play(FadeIn(ecl_lbl), Create(orb), FadeIn(apse))

        def rows():
            deg = np.rad2deg
            return _panel([
                (f"{d['name']}", P.GREEN),
                (f"a = {d['a']:.0f} AU      size", P.FG),
                (f"e = {d['e']:.2f}          stretch", P.FG),
                (f"i = {deg(I.get_value()):.1f}°   tilt", P.TEAL),
                (f"Ω = {deg(NODE.get_value()):.1f}°   node", P.PURPLE),
                (f"ω = {deg(ARGP.get_value()):.1f}°   perihelion", P.ORANGE),
            ], anchor=np.array([2.6, 2.75, 0.0]), font_size=21, buff=0.24)

        panel = always_redraw(rows)
        self.play(FadeIn(panel))
        timing.hold_to_read(self, cap, settle=0.3)
        cap = _swap_caption(self, cap, "a and e give its size and shape. Now orient it.")
        timing.hold_to_read(self, cap, settle=0.2)

        # --- tilt: i ---
        def node_dir():
            n = NODE.get_value()
            return np.array([np.cos(n), np.sin(n), 0.0])

        def node_line():
            v = view()
            nd = node_dir()
            return DashedLine(v.p(-R * nd), v.p(R * nd), color=P.PURPLE, stroke_width=2,
                              dash_length=0.1).set_z_index(2)

        def asc_node():
            # where the orbit climbs through the ecliptic: true anomaly = -omega
            v = view()
            nu = -ARGP.get_value()
            r = a * (1 - e * e) / (1 + e * np.cos(nu))
            p = v.p(r * node_dir())
            lbl = layout.label("ascending\nnode", font_size=15, color=P.PURPLE)
            out = v.p(1.6 * r * node_dir()) - v.p((0, 0, 0))
            lbl.move_to(p + 0.75 * out / max(np.linalg.norm(out), 1e-6))
            return VGroup(Dot(p, radius=0.07, color=P.PURPLE), lbl).set_z_index(6)

        def tilt_arc():
            v = view()
            n = node_dir()
            perp = np.cross([0, 0, 1.0], n)
            i = I.get_value()
            c = a * (1 - e * e) / (1 + e * np.cos(-ARGP.get_value())) * n
            pts = [v.p(c + 0.55 * (np.cos(t) * perp + np.sin(t) * np.array([0, 0, 1.0])))
                   for t in np.linspace(0, i, 20)]
            g = VGroup(_poly(pts, P.TEAL, 3))
            g.add(layout.label("i", font_size=22, color=P.TEAL).move_to(
                v.p(c + 0.85 * (np.cos(i / 2) * perp + np.sin(i / 2) * np.array([0, 0, 1.0])))))
            # a tilt angle is only readable from the side: fade it as the camera rises
            fade = float(np.clip(np.cos(ELEV.get_value()) / np.cos(np.deg2rad(30.0)), 0.0, 1.0))
            g[0].set_stroke(opacity=fade)
            g[1].set_opacity(fade)
            return g.set_z_index(6)

        nline = always_redraw(node_line)
        anode = always_redraw(asc_node)
        iarc = always_redraw(tilt_arc)
        cap = _swap_caption(self, cap, "i: tip the orbit's plane out of the ecliptic, about a hinge.")
        self.add(nline, anode, iarc)
        self.play(I.animate.set_value(np.deg2rad(d["i_deg"])), run_time=2.5)
        timing.hold_to_read(self, cap, settle=0.3)

        # --- node: Omega ---
        def node_arc():
            v = view()
            n = NODE.get_value()
            pts = [v.p(p) for p in _arc_xyz(np.eye(3), 1.3, 0.0, n)]
            g = VGroup(_poly(pts, P.PURPLE, 3))
            g.add(layout.label("Ω", font_size=22, color=P.PURPLE).move_to(
                v.p((1.6 * np.cos(n / 2), 1.6 * np.sin(n / 2), 0.0))))
            return g.set_z_index(6)

        narc = always_redraw(node_arc)
        self.add(narc)
        cap = _swap_caption(self, cap, "Ω: swing the hinge round. It is measured from the reference direction.")
        self.play(NODE.animate.set_value(np.deg2rad(d["node_deg"])), run_time=3.0)
        timing.hold_to_read(self, cap, settle=0.3)

        # --- argument of perihelion: omega ---
        def argp_arc():
            v = view()
            fr = _rot_z(NODE.get_value()) @ _rot_x(I.get_value())
            w = ARGP.get_value()
            pts = [v.p(p) for p in _arc_xyz(fr, 0.8, 0.0, w)]
            g = VGroup(_poly(pts, P.ORANGE, 3))
            mid = _arc_xyz(fr, 1.05, w / 2, w / 2, n=1)[0]
            g.add(layout.label("ω", font_size=22, color=P.ORANGE).move_to(v.p(mid)))
            return g.set_z_index(6)

        warc = always_redraw(argp_arc)
        self.add(warc)
        cap = _swap_caption(self, cap, "ω: slide perihelion round within the tilted plane, starting at the node.")
        self.play(ARGP.animate.set_value(np.deg2rad(d["argp_deg"])), run_time=3.5)
        timing.hold_to_read(self, cap, settle=0.3)

        # --- look down from above: longitude of perihelion ---
        cap = _swap_caption(self, cap, "Now look straight down on the planets' plane.")
        self.play(FadeOut(ecl_lbl), ELEV.animate.set_value(np.pi / 2), run_time=3.0,
                  rate_func=rate_functions.smooth)
        vp = np.deg2rad(d["varpi_deg"])
        top = view()
        varpi_arc = _poly([top.p((2.1 * np.cos(t), 2.1 * np.sin(t), 0.0))
                           for t in np.linspace(0.0, vp, 40)], P.FG, 4).set_z_index(7)
        varpi_lbl = layout.label("ϖ", font_size=26, color=P.FG, weight="BOLD").move_to(
            top.p((2.45 * np.cos(vp / 2), 2.45 * np.sin(vp / 2), 0.0))).set_z_index(7)
        total = d["node_deg"] + d["argp_deg"]
        eq = VGroup(
            layout.equation(r"\varpi = \Omega + \omega", color=P.FG).scale(0.75),
            layout.equation(rf"= {d['node_deg']:.1f}^\circ + {d['argp_deg']:.1f}^\circ"
                            rf" = {total:.1f}^\circ", color=P.FG).scale(0.62),
            layout.equation(rf"\equiv {d['varpi_deg']:.1f}^\circ", color=P.FG).scale(0.62),
        ).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        eq.next_to(panel, DOWN, buff=0.4, aligned_edge=LEFT)
        eq_note = layout.label("longitude of perihelion", font_size=18, color=P.MUTED)
        eq_note.next_to(eq, DOWN, buff=0.15, aligned_edge=LEFT)
        cap = _swap_caption(self, cap, "ϖ = Ω + ω, a two-step angle: which way the orbit points.")
        self.play(Create(varpi_arc), FadeIn(varpi_lbl), FadeIn(eq), FadeIn(eq_note))
        timing.hold_to_read(self, cap, eq_note, settle=0.8)
        cap = _swap_caption(self, cap, "Planet Nine's evidence lives in ϖ, and in how the planes tilt.")
        timing.hold_to_read(self, cap, settle=0.5)

        # --- the sixth number: where along the orbit the body is ---
        for m in (sun, plane, ref, ref_lbl, orb, apse, nline, anode, iarc, narc, warc, panel):
            m.clear_updaters()
        M = ValueTracker(0.0)
        fr = frame()

        def body():
            nu = orbits.nu_from_mean_anomaly(e, M.get_value())
            r = a * (1 - e * e) / (1 + e * np.cos(nu))
            return Dot(top.p(fr @ np.array([r * np.cos(nu), r * np.sin(nu), 0.0])), radius=0.09,
                       color=P.GREEN).set_z_index(8)

        dot = always_redraw(body)
        m_row = layout.label(f"M: where it is now\n(one lap: {d['period_yr']:,.0f} yr)",
                             font_size=19, color=P.GREEN)
        m_row.next_to(eq_note, DOWN, buff=0.35, aligned_edge=LEFT)
        cap = _swap_caption(self, cap, "The sixth number, M, says where along the orbit the body is.")
        self.add(dot)
        self.play(FadeIn(m_row))
        self.play(M.animate.set_value(2 * np.pi), run_time=5.0, rate_func=rate_functions.linear)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Shape: a, e.  Plane: i, Ω.  Pointing: ϖ = Ω + ω.  Position: M.")


# ---- P04 ----------------------------------------------------------------------

def _nice_length(target):
    """The 1-2-5 round number of AU closest (in log) to ``target`` AU."""
    best = None
    for k in range(-1, 5):
        for m in (1, 2, 5):
            v = m * 10 ** k
            if best is None or abs(np.log(v / target)) < abs(np.log(best / target)):
                best = v
    return best


def _light(hours):
    if hours < 1.0:
        return f"{hours * 60:.0f} light-minutes"
    if hours < 48.0:
        return f"{hours:.1f} light-hours"
    return f"{hours / 24:.1f} light-days"


def _window(x, lo, hi, soft=1.35):
    """Opacity ramp: 1 inside [lo, hi], fading to 0 a factor ``soft`` outside."""
    if x <= 0:
        return 0.0
    if x < lo:
        return float(np.clip(np.log(x / (lo / soft)) / np.log(soft), 0, 1))
    if x > hi:
        return float(np.clip(1 - np.log(x / hi) / np.log(soft), 0, 1))
    return 1.0


class P04Geography(Scene):
    """Goal: a to-scale tour of where everything is -- planets, Kuiper belt,
    the distant objects, Planet Nine -- and the (a, q) map that shows the
    distant objects sit out of Neptune's reach."""

    def construct(self):
        self.add(layout.concept_badge("PREFACE 04"))
        g = _data()["geography"]
        self._zoom(g)
        self._aq_map(g)

    def _zoom(self, g):
        LZ = ValueTracker(np.log10(1.0 / 2.0))          # log10(AU per scene unit)
        lz_end = np.log10(300.0)
        o0, o1 = np.array([-1.0, 0.0, 0.0]), np.array([-1.3, 0.55, 0.0])

        def origin():
            f = np.clip((LZ.get_value() - 1.4) / (lz_end - 1.4), 0.0, 1.0)
            return o0 + (o1 - o0) * f

        def to_scene(xy, s, o):
            return [o + np.array([p[0] / s, p[1] / s, 0.0]) for p in xy]

        # label angle per planet, alternating sides so neighbours never collide
        angles = {"Earth": 55, "Jupiter": 55, "Saturn": 125, "Uranus": 55, "Neptune": 125}
        light = g["light_hr_per_au"]
        DIST = ValueTracker(0.0)

        def world():
            s = 10 ** LZ.get_value()
            o = origin()
            out = VGroup()
            # Kuiper belt
            k_in, k_out = [x / s for x in g["kuiper_belt_au"]]
            if 0.06 < k_out < 14:
                ring = Annulus(inner_radius=k_in, outer_radius=k_out, stroke_width=0,
                               fill_color=P.GREEN, fill_opacity=0.07).move_to(o)
                out.add(ring)
                w = _window(k_out, 1.2, 2.6)
                if w > 0:
                    out.add(layout.label("Kuiper belt, 30-50 AU", font_size=18, color=P.GREEN)
                            .set_opacity(w).move_to(o + np.array([0.0, -k_out - 0.25, 0.0])))
            for pl in g["planets"]:
                r = pl["a"] / s
                if 0.04 < r < 12:
                    fade = float(np.clip((r - 0.04) / 0.1, 0, 1))
                    out.add(Circle(radius=r, color=P.FG, stroke_width=1.8)
                            .move_to(o).set_stroke(opacity=fade))
                    w = _window(r, 0.9, 2.9)
                    if w > 0:
                        t = np.deg2rad(angles[pl["name"]])
                        pos = o + r * np.array([np.cos(t), np.sin(t), 0.0])
                        side = UP + (RIGHT if np.cos(t) > 0 else LEFT)
                        out.add(layout.label(pl["name"], font_size=19).set_opacity(w)
                                .next_to(pos, side, buff=0.05))
            pl = g["pluto"]
            pts = to_scene(pl["xyz"], s, o)
            rmax = max(np.linalg.norm(p - o) for p in pts)
            if 0.05 < rmax < 14:
                out.add(_poly(pts, P.GREEN, 1.8, 0.9))
                w = _window(rmax, 1.0, 2.8)
                if w > 0:
                    out.add(layout.label("Pluto", font_size=18, color=P.GREEN).set_opacity(w)
                            .next_to(o + LEFT * rmax, LEFT, buff=0.1))
            vr = g["voyager1_au"] / s
            if 0.2 < vr < 14:
                out.add(DashedVMobject(Circle(radius=vr, color=P.PURPLE, stroke_width=1.6)
                                       .move_to(o), num_dashes=60))
                w = _window(vr, 0.9, 2.2)
                if w > 0:
                    out.add(layout.label(f"Voyager 1 today: {g['voyager1_au']:.0f} AU", font_size=17,
                                         color=P.PURPLE).set_opacity(w)
                            .next_to(o + UP * vr, UP, buff=0.08))
            w = DIST.get_value()
            if w > 0.01:
                for obj in g["distant"]:
                    out.add(_poly(to_scene(obj["xyz"], s, o), P.GREEN, 1.4, 0.75 * w))
            out.add(orbits.sun(radius=0.08).move_to(o))
            return out

        def scale_bar():
            s = 10 ** LZ.get_value()
            L = _nice_length(1.4 * s)
            left = np.array([-6.6, -2.75, 0.0])
            bar = Line(left, left + RIGHT * L / s, color=P.FG, stroke_width=3)
            ends = VGroup(*[Line(p + UP * 0.08, p + DOWN * 0.08, color=P.FG, stroke_width=3)
                            for p in (bar.get_start(), bar.get_end())])
            txt = layout.label(f"{L:g} AU  =  {_light(L * light)}", font_size=18)
            txt.next_to(bar, UP, buff=0.12, aligned_edge=LEFT)
            return VGroup(bar, ends, txt)

        view = always_redraw(world)
        bar = always_redraw(scale_bar)
        self.add(view, bar)
        cap = _swap_caption(self, None, "Everything to scale. Earth: 1 AU, the Earth-Sun distance.")
        timing.hold_to_read(self, cap, settle=0.3)

        a_jup = next(p["a"] for p in g["planets"] if p["name"] == "Jupiter")
        stops = [
            (np.log10(a_jup / 1.9),
             f"Jupiter at {a_jup:.1f} AU: sunlight takes {a_jup * light * 60:.0f} minutes to get there."),
            (np.log10(g["neptune_a"] / 2.4),
             f"Neptune, the last planet: 30 AU, {g['light_hr_neptune']:.1f} light-hours from the Sun."),
            (np.log10(g["kuiper_belt_au"][1] / 2.2),
             "The Kuiper belt: icy bodies from 30 to 50 AU, Pluto among them."),
            (np.log10(g["voyager1_au"] / 1.6),
             f"Voyager 1, our farthest probe: only {g['voyager1_au']:.0f} AU out after "
             f"{g['voyager1_years_flown']:.0f} years."),
            (lz_end, "Beyond: the distant objects, on orbits hundreds of AU long."),
        ]
        for lz, text in stops:
            anims = [LZ.animate.set_value(lz)]
            if lz == lz_end:
                anims.append(DIST.animate.set_value(1.0))
            self.play(*anims, run_time=2.6, rate_func=rate_functions.smooth)
            cap = _swap_caption(self, cap, text)
            timing.hold_to_read(self, cap, settle=0.2)

        view.clear_updaters()
        bar.clear_updaters()
        s = 10 ** lz_end
        o = origin()
        ref = next(p for p in g["p9"] if p["label"].startswith("Batygin et al. 2019"))
        p9 = DashedVMobject(_poly([o + np.array([p[0] / s, p[1] / s, 0.0]) for p in ref["xyz"]],
                                  P.BLUE, 3.2), num_dashes=80).set_z_index(4)
        p9_lbl = layout.label("Planet Nine?", font_size=20, color=P.BLUE, weight="BOLD")
        top = max((pt for pt in p9.get_all_points()), key=lambda q: q[1])
        p9_lbl.next_to(top, RIGHT, buff=0.2)
        cap = _swap_caption(self, cap, f"Planet Nine's proposed orbit: a = {ref['a']:.0f} AU.")
        self.play(Create(p9), FadeIn(p9_lbl), run_time=2.2)
        panel = _panel([
            ("sunlight takes", P.MUTED),
            (f"to Neptune:  {g['light_hr_neptune']:.1f} hours", P.FG),
            (f"to Planet Nine:  {g['light_day_p9']:.1f} days", P.BLUE),
            ("sunlight there is", P.MUTED),
            (f"1/{1 / g['sunlight_p9_vs_neptune']:.0f} as bright as at Neptune", P.BLUE),
            (f"Voyager 1 would need {g['voyager1_yr_to_p9']:.0f} yr", P.PURPLE),
        ], anchor=np.array([2.7, 2.3, 0.0]), font_size=20, buff=0.2)
        self.play(FadeIn(panel, shift=LEFT * 0.1))
        timing.hold_to_read(self, cap, panel, settle=0.8)
        self.play(*[FadeOut(m) for m in (view, bar, p9, p9_lbl, panel, cap)])

    def _aq_map(self, g):
        # log10(a) in [2, 3.5] and q in [25, 85] AU, offset so the axes meet bottom-left
        ax = Axes(x_range=[0, 1.5, 0.5], y_range=[0, 60, 10], x_length=8.2, y_length=4.5,
                  axis_config={"color": P.MUTED, "include_tip": False})
        ax.move_to(np.array([-1.3, 0.3, 0.0]))

        def c2p(a_au, q_au):
            return ax.c2p(np.log10(a_au) - 2.0, q_au - 25.0)

        ticks = VGroup()
        for a_t in (100, 200, 500, 1000, 2000):
            ticks.add(layout.label(f"{a_t:,}", font_size=17).next_to(c2p(a_t, 25), DOWN, buff=0.12))
        for q_t in range(30, 86, 10):
            ticks.add(layout.label(str(q_t), font_size=17).next_to(c2p(100, q_t), LEFT, buff=0.12))
        xl = layout.label("orbit size a (AU)", font_size=19).next_to(ticks, DOWN, buff=0.15)
        xl.set_x(ax.get_center()[0])
        yl = layout.label("closest approach q (AU)", font_size=19).rotate(np.pi / 2)
        yl.next_to(ax, LEFT, buff=0.55)

        qn = g["neptune_a"]
        qd = g["detached_q_min"]
        reach = Polygon(c2p(100, 25), c2p(3162, 25), c2p(3162, qd), c2p(100, qd),
                        stroke_width=0).set_fill(P.ORANGE, opacity=0.14)
        nep = DashedLine(c2p(100, qn), c2p(3162, qn), color=P.FG, stroke_width=1.8)
        strip = (25.0 + qn) / 2
        nep_lbl = layout.label("Neptune's orbit", font_size=16).move_to(c2p(1800, strip))
        reach_lbl = layout.label(f"q < {qd:.0f} AU: Neptune still tugs on these", font_size=16,
                                 color=P.ORANGE).move_to(c2p(300, strip))
        det_lbl = layout.label("detached: out of Neptune's reach", font_size=17, color=P.GREEN)
        det_lbl.move_to(c2p(200, 57))

        cap = _swap_caption(self, None, "A map of the distant objects: orbit size against closest approach.")
        self.play(Create(ax), FadeIn(ticks), FadeIn(xl), FadeIn(yl))
        self.play(FadeIn(reach), Create(nep), FadeIn(nep_lbl), FadeIn(reach_lbl))
        dots = VGroup(*[Dot(c2p(o["a"], o["q"]), radius=0.065, color=P.GREEN) for o in g["aq"]])
        self.play(FadeIn(dots, lag_ratio=0.08), run_time=2.0)
        timing.hold_to_read(self, cap, settle=0.3)

        named = {"Sedna": RIGHT, "2012 VP113": RIGHT, "2015 TG387": UP,
                 "Ammonite (2023 KQ14)": LEFT}
        tags = VGroup()
        for o in g["aq"]:
            if o["name"] in named:
                name = "Ammonite" if o["name"].startswith("Ammonite") else o["name"]
                tags.add(layout.label(name, font_size=16, color=P.GREEN)
                         .next_to(c2p(o["a"], o["q"]), named[o["name"]], buff=0.1))
        cap = _swap_caption(self, cap, "Many never come near Neptune. Neptune alone cannot have put them there.")
        self.play(FadeIn(det_lbl), FadeIn(tags))
        timing.hold_to_read(self, cap, settle=0.5)

        q9 = min(p["q"] for p in g["p9"])
        q9_hi = max(p["q"] for p in g["p9"])
        base = c2p(500, 78)
        up = Arrow(base, base + UP * 1.0, buff=0, color=P.BLUE, stroke_width=5, tip_length=0.2)
        up_lbl = layout.label(f"Planet Nine: q = {q9:.0f}-{q9_hi:.0f} AU, far off the top",
                              font_size=17, color=P.BLUE)
        up_lbl.next_to(up, RIGHT, buff=0.15).shift(UP * 0.3)
        cap = _swap_caption(self, cap, "Something beyond them could have lifted and shaped these orbits.")
        self.play(GrowArrow(up), FadeIn(up_lbl))
        timing.hold_to_read(self, cap, settle=0.8)
        self.play(FadeOut(cap))
        layout.show_takeaway(self, "Past Neptune lie orbits hundreds of AU long: Planet Nine's realm.")
