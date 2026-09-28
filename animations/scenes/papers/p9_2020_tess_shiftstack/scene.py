"""Rice & Laughlin (2020) -- a targeted TESS shift-stacking search for Planet
Nine and distant TNOs in the galactic plane.

TESS pixels are 21 arcsec wide and a sector lasts 27 days, so what a
shift-stack can find is set by how many pixels a body crosses: nearby objects
race across dozens, Planet Nine crawls across a handful. The pipeline tries
every shift vector in a grid of on-sky rates, recovers three known TNOs, and
finds that blind recovery is reliable only for V < 21 within 150 AU.
Reproduced in p9-2020-tess-shiftstack: the drift, trial-track grid and depth
curves are the crate's own (anim.json -> papers -> p9-2020-tess-shiftstack).
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
    Rectangle,
    Scene,
    SurroundingRectangle,
    ValueTracker,
    VGroup,
    always_redraw,
)

import p9_manim as P
from p9_manim.plot import Plot
from p9_manim import dataio, layout, paper, timing

CRATE = "p9-2020-tess-shiftstack"


def _readout(title, value, note, colour):
    rows = [layout.label(title, font_size=16, color=P.FG),
            layout.label(value, font_size=28, color=colour, weight="BOLD")]
    if note:
        rows.append(layout.label(note, font_size=15, color=P.MUTED))
    g = VGroup(*rows).arrange(DOWN, buff=0.1)
    box = SurroundingRectangle(g, color=colour, buff=0.18, corner_radius=0.08)
    box.set_fill(colour, opacity=0.06)
    return VGroup(box, g)


class TessShiftstack2020(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        p9 = d["p9"]
        dist = np.array(d["distance_au"])
        drift = np.array(d["drift_px_per_sector"])
        rate = np.array(d["rate_arcsec_per_day"])
        self.add(paper.scene_header(CRATE))

        # 1. how far a body crawls across TESS pixels in one sector
        n_px, w = 64, 0.19
        strip = VGroup(*[Rectangle(width=w, height=w, stroke_width=0.8, color=P.MUTED)
                         .set_fill("#16171f", 1.0) for _ in range(n_px)]).arrange(RIGHT, buff=0)
        strip.move_to(UP * 0.4)
        x0 = strip[0].get_left()[0]
        moon_px = 1865.0 / d["pixel_scale_arcsec"]  # full Moon: 31 arcmin, a scale reference
        strip_lab = layout.label(
            f"TESS pixels, {d['pixel_scale_arcsec']:.0f} arcsec each: the full Moon spans "
            f"{moon_px:.0f} of them", font_size=17, color=P.MUTED)
        strip_lab.next_to(strip, UP, buff=0.15).align_to(strip, LEFT)
        t = ValueTracker(float(dist[0]))

        def track():
            n = float(np.interp(t.get_value(), dist, drift))
            y = strip.get_center()[1]
            col = P.GREEN if t.get_value() <= d["published_blind_distance_au"] else P.BLUE
            bar = Line([x0, y, 0], [x0 + n * w, y, 0], color=col, stroke_width=9)
            return VGroup(bar, Dot([x0 + n * w, y, 0], radius=0.08, color=col))

        mover = always_redraw(track)
        tag = always_redraw(lambda: layout.label(
            f"{t.get_value():.0f} AU away: crosses "
            f"{np.interp(t.get_value(), dist, drift):.0f} pixels in a "
            f"{d['sector_days']:.0f}-day sector",
            font_size=20, color=P.FG).next_to(strip, DOWN, buff=0.4))
        cap = layout.caption("The farther the body, the slower it crawls across the image",
                             font_size=22)
        self.play(FadeIn(strip), FadeIn(strip_lab), FadeIn(mover), FadeIn(tag), FadeIn(cap))
        self.play(t.animate.set_value(p9["perihelion_au"]), run_time=3.0)
        timing.hold_to_read(self, cap, settle=0.2)
        self.play(t.animate.set_value(p9["aphelion_au"]), run_time=1.5)
        lo = np.interp(p9["aphelion_au"], dist, drift)
        hi = np.interp(p9["perihelion_au"], dist, drift)
        zone = Rectangle(width=(hi - lo) * w, height=w * 1.8, stroke_width=0).set_fill(P.BLUE, 0.3)
        zone.move_to([x0 + 0.5 * (lo + hi) * w, strip.get_center()[1], 0])
        zone_lab = layout.label(
            f"Planet Nine at {p9['perihelion_au']:.0f}-{p9['aphelion_au']:.0f} AU: "
            f"{lo:.0f}-{hi:.0f} pixels", font_size=17, color=P.BLUE)
        zone_lab.next_to(strip, DOWN, buff=1.05)
        near = np.interp(d["published_blind_distance_au"], dist, drift)
        near_mark = Line([x0 + near * w, strip.get_bottom()[1] - 0.1, 0],
                         [x0 + near * w, strip.get_top()[1] + 0.1, 0], color=P.GREEN,
                         stroke_width=3)
        near_lab = layout.label(f"{d['published_blind_distance_au']:.0f} AU: {near:.0f} pixels",
                                font_size=17, color=P.GREEN)
        near_lab.next_to(zone_lab, DOWN, buff=0.2)
        cap2 = layout.caption(
            "Planet Nine barely moves: little to separate it from the fixed stars", font_size=22)
        self.play(FadeIn(zone), FadeIn(zone_lab), Create(near_mark), FadeIn(near_lab),
                  FadeOut(cap), FadeIn(cap2))
        timing.hold_to_read(self, cap2, zone_lab, settle=0.6)
        self.remove(mover, tag)
        self.play(FadeOut(VGroup(strip, strip_lab, zone, zone_lab, near_mark, near_lab, cap2)))

        # 2. the trial tracks: a grid of on-sky velocities
        r_max = 50.0
        n_side = int(round(np.sqrt(d["trial_tracks"])))
        half = 2.35
        s = half / r_max
        centre = np.array([-2.6, 0.2, 0])
        frame = Rectangle(width=2 * half, height=2 * half, color=P.PURPLE, stroke_width=1.6)
        frame.move_to(centre)
        grid = VGroup()
        for k in range(1, n_side):
            u = -half + 2 * half * k / n_side
            grid.add(Line(centre + [u, -half, 0], centre + [u, half, 0]),
                     Line(centre + [-half, u, 0], centre + [half, u, 0]))
        grid.set_stroke(P.PURPLE, width=0.6, opacity=0.45)
        ax_x = layout.label("east-west rate", font_size=15, color=P.MUTED)
        ax_x.next_to(frame, DOWN, buff=0.12)
        ax_y = layout.label("north-south rate", font_size=15, color=P.MUTED).rotate(np.pi / 2)
        ax_y.next_to(frame, LEFT, buff=0.12)
        scale_lab = layout.label(f"±{r_max:.0f} arcsec/day", font_size=14, color=P.PURPLE)
        scale_lab.next_to(frame, UP, buff=0.1).align_to(frame, RIGHT)
        cap3 = layout.caption("Every cell is one shift vector to stack along", font_size=22)
        self.play(Create(frame), FadeIn(ax_x), FadeIn(ax_y), FadeIn(scale_lab), FadeIn(cap3))
        self.play(Create(grid, lag_ratio=0.02), run_time=1.6)

        def ring(au, col, width=2.5):
            r = float(np.interp(au, dist, rate)) * s
            return Circle(radius=r, color=col, stroke_width=width).move_to(centre)

        r70 = ring(dist[0], P.GREEN)
        r150 = ring(d["published_blind_distance_au"], P.GREEN)
        r_in, r_out = ring(p9["perihelion_au"], P.BLUE), ring(p9["aphelion_au"], P.BLUE)
        star = Dot(centre, radius=0.07, color=P.FG)
        keys = VGroup(
            layout.label(f"{dist[0]:.0f} AU: {rate[0]:.0f} arcsec/day", font_size=17,
                         color=P.GREEN),
            layout.label(f"{d['published_blind_distance_au']:.0f} AU: "
                         f"{np.interp(d['published_blind_distance_au'], dist, rate):.0f} arcsec/day",
                         font_size=17, color=P.GREEN),
            layout.label(f"Planet Nine: {np.interp(p9['aphelion_au'], dist, rate):.0f}-"
                         f"{np.interp(p9['perihelion_au'], dist, rate):.0f} arcsec/day",
                         font_size=17, color=P.BLUE),
            layout.label("stars: zero motion", font_size=17, color=P.FG),
        ).arrange(DOWN, buff=0.22, aligned_edge=LEFT)
        keys.move_to([3.0, 1.35, 0]).align_to(np.array([0.8, 0, 0]), LEFT)
        self.play(Create(r70), Create(r150), Create(r_in), Create(r_out), FadeIn(star),
                  FadeIn(keys, lag_ratio=0.2), run_time=1.6)
        tracks = _readout("trial tracks per sector", f"{d['trial_tracks']:.0f}",
                          f"paper: {d['published_trial_tracks']:.0f}", P.PURPLE)
        tracks.next_to(keys, DOWN, buff=0.45).align_to(keys, LEFT)
        cap4 = layout.caption(
            "Planet Nine's rate sits in the few cells nearest the stars", font_size=22)
        self.play(FadeIn(tracks), FadeOut(cap3), FadeIn(cap4))
        timing.hold_to_read(self, cap4, keys, settle=0.8)
        self.play(FadeOut(VGroup(frame, grid, ax_x, ax_y, scale_lab, r70, r150, r_in, r_out,
                                 star, keys, tracks, cap4)))

        # 3. bright enough, but outside the blind pipeline's range
        v = np.array(d["v_p9"])
        plot = Plot([70, 800], [23.5, 12.0], [100, 200, 300, 400, 500, 600, 700, 800],
                    [22, 20, 18, 16, 14, 12], "distance from Earth (AU)",
                    "V magnitude (brighter upward)", centre=(0.35, 0.45), height=4.3)
        ok = plot.band(70, d["published_blind_distance_au"], d["published_blind_limit_v"], 12.0,
                       P.GREEN, opacity=0.18)
        ok_lab = layout.label("blind recovery reliable (paper)", font_size=15, color=P.GREEN)
        ok_lab.next_to(plot.p(152, 12.4), RIGHT, buff=0.12)
        curve = plot.curve(dist, v, P.BLUE)
        curve_lab = layout.label(f"Planet Nine ({p9['mass_earth']:.1f} Earth masses)",
                                 font_size=16, color=P.BLUE)
        curve_lab.next_to(plot.p(475, 18.4), RIGHT, buff=0.1)
        depth = plot.hline(d["depth_one_sector"], P.PURPLE)
        depth_lab = layout.label(f"one sector reaches V = {d['depth_one_sector']:.1f}",
                                 font_size=15, color=P.PURPLE)
        depth_lab.next_to(plot.p(165, d["depth_one_sector"]), DOWN + RIGHT, buff=0.1)
        self.play(FadeIn(plot))
        self.play(Create(curve), FadeIn(curve_lab), Create(depth), FadeIn(depth_lab),
                  run_time=1.4)
        tnos = VGroup(layout.label("pipeline test, recovered:", font_size=15, color=P.GREEN),
                      *[layout.label(f"{o['name']}  V = {o['v_mag']:.2f}", font_size=15,
                                     color=P.GREEN) for o in d["tnos"]])
        tnos.arrange(DOWN, buff=0.1, aligned_edge=LEFT).move_to(plot.p(690, 13.6))
        cap5 = layout.caption("Blind, the pipeline finds V < 21 bodies only within 150 AU",
                              font_size=22)
        self.play(FadeIn(ok), FadeIn(ok_lab), FadeIn(tnos), FadeIn(cap5))
        timing.hold_to_read(self, cap5, tnos, settle=0.4)
        zone = plot.band(p9["perihelion_au"], p9["aphelion_au"], 23.5, 12.0, P.BLUE,
                         opacity=0.14)
        zone_lab2 = layout.label("its orbit", font_size=15, color=P.BLUE)
        zone_lab2.next_to(plot.p(0.5 * (p9["perihelion_au"] + p9["aphelion_au"]), 23.5), UP,
                          buff=0.12)
        cap6 = layout.caption(
            f"Planet Nine would be V = {p9['v_perihelion']:.1f}-{p9['v_aphelion']:.1f}: "
            "bright enough, but too slow for the blind search", font_size=22)
        self.play(FadeIn(zone), FadeIn(zone_lab2), FadeOut(cap5), FadeIn(cap6))
        timing.hold_to_read(self, cap6, settle=1.0)
        self.play(FadeOut(cap6))

        layout.show_takeaway(
            self, "A working TESS pipeline, but Planet Nine's distance is still beyond its reach.")
