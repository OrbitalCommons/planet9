"""Pfalzner, Wagner & Bischoff (2026) -- the flyby, 4.5 billion years on.

A 0.8 solar-mass star on a parabolic orbit passes 110 AU from the young Sun,
inclined 70 degrees to its disc. The scene flies the crate's exact perturber
trajectory past the crate's disc of tracers, then morphs each drawn tracer
from its circular pre-flyby orbit to the orbit the crate's integration left it
on, coloured by the Table 1 family it falls in. It then sets the computed
family census against Table 1, and follows the tracers left with perihelia near
Neptune through the crate's reduced-scale long-term integration.
Data: anim.json -> papers -> p9-2026-flyby-evolution.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Circle,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    LaggedStart,
    Line,
    Polygon,
    Scene,
    Transform,
    ValueTracker,
    VGroup,
    VMobject,
    always_redraw,
    linear,
)

import p9_manim as P
from p9_manim import dataio, layout, orbits, paper, timing, widgets

CRATE = "p9-2026-flyby-evolution"

FAMILY = {
    "cold KB": ("Kuiper belt", P.GREEN),
    "hot KB": ("Kuiper belt", P.GREEN),
    "detached": ("detached and Sedna-like", P.TEAL),
    "Sedna-like": ("detached and Sedna-like", P.TEAL),
    "inclined": ("steeply inclined or retrograde", P.PURPLE),
    "retrograde": ("steeply inclined or retrograde", P.PURPLE),
    "inner": ("thrown into the planets' region", P.RED),
    None: ("in no Table 1 family", P.MUTED),
}
BAR_COLOUR = {
    "cold KB": P.GREEN, "hot KB": P.GREEN, "detached": P.TEAL, "Sedna-like": P.TEAL,
    "inclined": P.PURPLE, "retrograde": P.PURPLE, "inner": P.RED,
}

# view of the disc frame: turned by AZ about the disc pole, tilted EL toward us
AZ, EL = np.radians(30.0), np.radians(28.0)
SCALE = 3.5 / 300.0  # scene units per AU
CENTRE = np.array([-1.6, 0.35, 0.0])
CLIP_AU = 280.0
# drawn post-flyby orbits are cut where they run beyond this distance
ORBIT_CLIP_AU = 230.0


def view(x, y, z):
    """Project disc-frame AU coordinates onto the screen."""
    xr = x * np.cos(AZ) - y * np.sin(AZ)
    yr = x * np.sin(AZ) + y * np.cos(AZ)
    return CENTRE + SCALE * np.array([xr, yr * np.sin(EL) + z * np.cos(EL), 0.0])


def polyline(points, color, stroke_width=1.6, opacity=0.8):
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners(points)
    return m.set_stroke(opacity=opacity)


def ring(r, n, color, stroke_width=1.4, opacity=0.6):
    t = np.linspace(0, 2 * np.pi, n)
    return polyline([view(r * np.cos(u), r * np.sin(u), 0.0) for u in t], color,
                    stroke_width, opacity)


def pct(x):
    return f"{100 * x:.1f}%" if x < 0.1 else f"{100 * x:.0f}%"


class FlybyEvolution2026(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        fb = d["flyby"]
        disc = d["disc"]

        self.add(paper.scene_header(CRATE))

        # 1. the encounter: the star's real path past the disc of tracers
        sun = orbits.sun(radius=0.07).move_to(view(0, 0, 0)).set_z_index(5)
        nep_au = dataio.body("Neptune")["a_au"]
        nep = ring(nep_au, 64, P.FG, stroke_width=1.6, opacity=0.9)
        nep_lab = layout.label("Neptune's orbit (white)", font_size=14, color=P.FG)
        nep_lab.move_to(view(0, -disc["r_max_au"], 0) + DOWN * 0.45 + LEFT * 1.2)
        nep_lab = VGroup(nep_lab, Line(nep_lab.get_right() + RIGHT * 0.05,
                                       view(nep_au * 0.5, -nep_au * 0.87, 0), color=P.FG,
                                       stroke_width=1.0).set_stroke(opacity=0.7))
        edges = VGroup(ring(disc["r_min_au"], 72, P.MUTED, 1.0, 0.5),
                       ring(disc["r_max_au"], 96, P.MUTED, 1.0, 0.5))
        drawn = d["drawn_orbits"]
        before = VGroup(*[ring(o["r0_au"], len(o["x"]), P.GREEN, 1.2, 0.55) for o in drawn])
        disc_lab = layout.label(
            f"the young disc: {disc['r_min_au']:.0f}–{disc['r_max_au']:.0f} AU, circular orbits",
            font_size=15, color=P.GREEN)
        disc_lab.next_to(view(0, disc["r_max_au"], 0), UP, buff=0.3)

        path = fb["path"]
        pts = [view(p["x"], p["y"], p["z"]) for p in path
               if np.linalg.norm([p["x"], p["y"], p["z"]]) <= CLIP_AU]
        track = polyline(pts, P.ORANGE, stroke_width=2.0, opacity=0.9)
        k = ValueTracker(0.0)
        star = always_redraw(lambda: Dot(pts[min(int(k.get_value()), len(pts) - 1)], radius=0.1,
                                         color=P.ORANGE).set_z_index(6))
        trail = always_redraw(lambda: polyline(pts[:max(2, int(k.get_value()) + 1)], P.ORANGE,
                                               stroke_width=2.0, opacity=0.9))
        peri = view(fb["periastron"]["x"], fb["periastron"]["y"], fb["periastron"]["z"])
        peri_dot = Dot(peri, radius=0.05, color=P.ORANGE)
        peri_lab = layout.label(f"closest: {fb['q_au']:.0f} AU", font_size=15, color=P.ORANGE)
        peri_lab.next_to(peri, UP + RIGHT, buff=0.08)

        facts = VGroup(
            layout.label("the passing star", font_size=17, color=P.ORANGE, weight="BOLD"),
            layout.label(f"{fb['mass_solar']:.1f} solar masses", font_size=16),
            layout.label(f"closest approach {fb['q_au']:.0f} AU", font_size=16),
            layout.label(f"path tilted {fb['inclination_deg']:.0f}° to the disc", font_size=16),
            layout.label(f"{fb['periastron_speed_kms']:.1f} km/s at closest", font_size=16),
            layout.label(f"crossing takes ~{2 * abs(path[0]['t_yr']):,.0f} yr", font_size=16,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.15, aligned_edge=LEFT).move_to([4.9, 1.2, 0])

        cap = layout.caption("Before: a flat disc of small bodies around the young Sun",
                             font_size=22)
        self.play(FadeIn(sun), FadeIn(edges), Create(nep), FadeIn(nep_lab),
                  LaggedStart(*[Create(r) for r in before], lag_ratio=0.01), FadeIn(disc_lab),
                  FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, disc_lab, settle=0.4)
        cap2 = layout.caption("A star sweeps past on its exact computed path", font_size=22)
        self.play(FadeOut(cap), FadeIn(cap2), FadeIn(facts), FadeOut(disc_lab), run_time=0.6)
        self.add(trail, star)
        i_peri = int(np.argmin([np.linalg.norm(p - peri) for p in pts]))
        self.play(k.animate.set_value(i_peri), run_time=3.0, rate_func=linear)
        self.play(FadeIn(peri_dot), FadeIn(peri_lab), run_time=0.4)

        after = VGroup()
        for o in drawn:
            col = FAMILY[o["group"]][1]
            p3 = [view(x, y, z) if np.sqrt(x * x + y * y + z * z) <= ORBIT_CLIP_AU
                  else None for x, y, z in zip(o["x"], o["y"], o["z"])]
            last = next(p for p in p3 if p is not None)
            filled = []
            for p in p3:
                last = p if p is not None else last
                filled.append(last)
            after.add(polyline(filled, col, 1.3, 0.75))
        cap3 = layout.caption("After: computed orbits, coloured by the paper's families",
                              font_size=22)
        self.play(k.animate.set_value(len(pts) - 1),
                  *[Transform(b, a) for b, a in zip(before, after)],
                  FadeOut(cap2), FadeIn(cap3), run_time=3.5, rate_func=linear)
        self.remove(star)
        seen = []
        key_rows = []
        for g in ("cold KB", "detached", "inclined", "inner", None):
            name, col = FAMILY[g]
            if name in seen:
                continue
            seen.append(name)
            key_rows.append(VGroup(Line(np.zeros(3), RIGHT * 0.35, color=col, stroke_width=3),
                                   layout.label(name, font_size=15, color=col))
                            .arrange(RIGHT, buff=0.12))
        key = VGroup(*key_rows).arrange(DOWN, buff=0.14, aligned_edge=LEFT)
        key.move_to([4.9, -1.3, 0])
        unb = VGroup(
            layout.label(f"{pct(d['unbound_fraction'])} of tracers torn away entirely",
                         font_size=15, color=P.MUTED),
            layout.label(f"orbits drawn out to {ORBIT_CLIP_AU:.0f} AU", font_size=15,
                         color=P.MUTED),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT).next_to(key, DOWN, buff=0.2,
                                                               aligned_edge=LEFT)
        self.play(FadeIn(key), FadeIn(unb), run_time=0.6)
        timing.hold_to_read(self, cap3, key, unb, settle=1.0)
        self.play(FadeOut(VGroup(sun, edges, nep, nep_lab, before, track, trail, peri_dot, peri_lab,
                                 facts, key, unb, cap3)), run_time=0.7)

        # 2. the family census against the paper's Table 1
        fams = d["families"]
        x_bar = -2.3
        per = 16.0  # scene units per unit fraction
        rows = VGroup()
        ticks = VGroup()
        for j, f in enumerate(fams):
            y = 2.25 - 0.66 * j
            col = BAR_COLOUR[f["label"]]
            lab = layout.label(f["label"], font_size=17, color=col).next_to([x_bar, y, 0], LEFT,
                                                                            buff=0.25)
            bar = Polygon([x_bar, y - 0.2, 0], [x_bar + per * f["fraction"], y - 0.2, 0],
                          [x_bar + per * f["fraction"], y + 0.2, 0], [x_bar, y + 0.2, 0],
                          stroke_width=0).set_fill(col, opacity=0.7)
            tick = Line([x_bar + per * f["published_t0"], y - 0.3, 0],
                        [x_bar + per * f["published_t0"], y + 0.3, 0], color=P.FG,
                        stroke_width=3)
            nums = VGroup(
                layout.label(pct(f["fraction"]), font_size=16, color=col),
                layout.label(pct(f["published_t0"]), font_size=16, color=P.FG),
            ).arrange(RIGHT, buff=0.5)
            nums.move_to([5.6, y, 0])
            rows.add(VGroup(lab, bar, nums))
            ticks.add(tick)
        heads = VGroup(layout.label("here", font_size=15, color=P.GREEN),
                       layout.label("paper", font_size=15, color=P.FG)).arrange(RIGHT, buff=0.55)
        heads.move_to([5.6, 2.85, 0])
        base = Line([x_bar, 2.65, 0], [x_bar, -2.0, 0], color=P.MUTED, stroke_width=1.2)
        tick_key = VGroup(Line(UP * 0.15, DOWN * 0.15, color=P.FG, stroke_width=3),
                          layout.label("white tick: the paper's Table 1, just after the flyby",
                                       font_size=15, color=P.FG)).arrange(RIGHT, buff=0.15)
        tick_key.move_to([0.4, -2.45, 0])
        cap4 = layout.caption(
            f"Computed from {d['n_tracers']:,} tracers: which family each one ends up in",
            font_size=22)
        self.play(FadeIn(base), FadeIn(heads), FadeIn(cap4), run_time=0.5)
        self.play(LaggedStart(*[FadeIn(r) for r in rows], lag_ratio=0.12), run_time=1.6)
        self.play(LaggedStart(*[Create(t) for t in ticks], lag_ratio=0.08), FadeIn(tick_key),
                  run_time=1.0)
        timing.hold_to_read(self, cap4, settle=0.6)
        cap5 = layout.caption(
            "Belt, Sedna-like and injected shares land near the paper's; too few steep orbits here",
            font_size=22)
        self.play(FadeOut(cap4), FadeIn(cap5), run_time=0.5)
        timing.hold_to_read(self, cap5, settle=1.0)
        self.play(FadeOut(VGroup(rows, ticks, heads, base, tick_key, cap5)), run_time=0.7)

        # 3. then 4.5 billion years of the planets: the near-Neptune tracers
        ev = d["evolution"]
        ax, labs = widgets.labeled_axes(
            [0, 150, 25], [0, 45, 10], x_label="semi-major axis  (AU)",
            y_label="perihelion  (AU)", y_rotate=True, numbers=True,
            x_length=6.0, y_length=4.0, shift_down=0)
        ax.move_to([-3.1, 0.4, 0])
        labs[0].next_to(ax, DOWN, buff=0.45)
        labs[1].next_to(ax, LEFT, buff=0.45)
        band = Polygon(ax.c2p(0, 30), ax.c2p(150, 30), ax.c2p(150, 40), ax.c2p(0, 40),
                       stroke_width=0).set_fill(P.ORANGE, opacity=0.14)
        band_lab = layout.label("perihelion near Neptune: 30–40 AU", font_size=14,
                                color=P.ORANGE).next_to(ax.c2p(75, 40), UP, buff=0.08)
        movers = VGroup()
        for tr in ev["tracers"]:
            s = tr["start"]
            if s["a_au"] > 150:
                continue
            dot = Dot(ax.c2p(s["a_au"], s["q_au"]), radius=0.04, color=P.GREEN)
            e = tr["end"]
            end = None if e is None else ax.c2p(min(e["a_au"], 150), max(e["q_au"], 0))
            movers.add(dot)
            dot.end = end
            dot.left_band = not tr["still_near_neptune"]
        t_ax, t_labs = widgets.labeled_axes(
            [0, ev["span_myr"], 1], [0, 0.7, 0.1], x_label="million years after the flyby",
            y_label="share that has left", y_rotate=True, numbers=True,
            x_length=4.4, y_length=4.0, shift_down=0)
        t_ax.move_to([3.9, 0.4, 0])
        t_labs[0].next_to(t_ax, DOWN, buff=0.45)
        t_labs[1].next_to(t_ax, LEFT, buff=0.45)
        tt = np.array(ev["t_myr"])
        loss = np.array(ev["loss"])
        tk = ValueTracker(0.0)
        loss_curve = always_redraw(lambda: widgets.curve(
            t_ax, tt[tt <= tk.get_value() + 1e-9], loss[tt <= tk.get_value() + 1e-9],
            color=P.RED) if tk.get_value() > tt[1] else VGroup())
        pub_lines = VGroup()
        for v, when in ((ev["published_loss_10myr"], "by 10 Myr"),
                        (ev["published_loss_4_5gyr"], "by 4.5 Gyr")):
            pub_lines.add(DashedLine(t_ax.c2p(0, v), t_ax.c2p(ev["span_myr"], v), color=P.FG,
                                     stroke_width=1.6))
            pub_lines.add(layout.label(f"paper: {pct(v)} {when}", font_size=14, color=P.FG)
                          .next_to(t_ax.c2p(0, v), UP, buff=0.06)
                          .align_to(t_ax.c2p(0, v), LEFT).shift(RIGHT * 0.1))
        cap6 = layout.caption(
            f"Follow {ev['n_followed']} tracers left with perihelia near Neptune",
            font_size=22)
        self.play(FadeIn(ax), FadeIn(labs), FadeIn(band), FadeIn(band_lab), FadeIn(movers),
                  FadeIn(t_ax), FadeIn(t_labs), FadeIn(cap6), run_time=1.0)
        timing.hold_to_read(self, cap6, settle=0.3)
        self.add(loss_curve)
        anims = []
        for dot in movers:
            if dot.end is not None:
                anims.append(dot.animate.move_to(dot.end).set_color(
                    P.RED if dot.left_band else P.GREEN))
        cap7 = layout.caption(
            f"Neptune's kicks start clearing the band: {pct(d['near_neptune_loss'])} gone in "
            f"{ev['span_myr']:.0f} Myr", font_size=22)
        self.play(*anims, tk.animate.set_value(tt[-1]), FadeOut(cap6), FadeIn(cap7),
                  run_time=4.0, rate_func=linear)
        self.play(FadeIn(pub_lines), run_time=0.5)
        timing.hold_to_read(self, cap7, pub_lines, settle=1.2)
        self.play(FadeOut(cap7), run_time=0.4)

        layout.show_takeaway(
            self, "Neptune slowly clears the 30–40 AU excess: the flyby fits even better")
