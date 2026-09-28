"""Brown, Holman & Batygin (2024) -- a Pan-STARRS1 search for Planet Nine.

The third search scored against the Brown & Batygin (2021) reference
population. Every synthetic planet is pushed through the ZTF, DES and PS1
survey models; PS1 (3pi sky north of -30 deg, V = 21.5 at 50% completeness,
nine detections required) removes the members the two earlier searches missed.
The survivors are what defines the updated orbit: more distant and fainter.
Everything drawn is the crate's own (anim.json -> papers -> p9-2024-panstarrs).
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
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2024-panstarrs"


def dot_key(items, font_size=16):
    """A row of (label, colour) keys with round swatches."""
    row = VGroup()
    for text, col in items:
        row.add(VGroup(Dot(radius=0.06, color=col),
                       layout.label(text, font_size=font_size, color=col)).arrange(RIGHT, buff=0.1))
    return row.arrange(RIGHT, buff=0.5)


def compare_row(name, before, after, published, unit, digits=0):
    """One line: prior median -> survivor median, with the paper's value."""
    cells = [
        layout.label(name, font_size=17, color=P.MUTED),
        layout.label(f"{before:.{digits}f}", font_size=17, color=P.MUTED),
        layout.label("→", font_size=17, color=P.MUTED),
        layout.label(f"{after:.{digits}f} {unit}", font_size=19, color=P.TEAL, weight="BOLD"),
        layout.label(f"paper: {published:.{digits}f}", font_size=15, color=P.FG),
    ]
    return cells


class PanStarrs2024(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pop = d["population"]
        ps1 = d["ps1"]
        pub = d["published"]

        self.add(paper.scene_header(CRATE))

        # 1. the state of the search the day before
        m = sky.SkyMap(width=11.0, dec_range=(-75, 75), centre=(0.0, 0.42, 0.0))
        ecl, gal = m.reference_curves()
        old = [s for s in pop if s["p_before"] >= 0.5]
        new = [s for s in pop if s["p_before"] < 0.5 and s["p_after"] >= 0.5]
        left = [s for s in pop if s["p_after"] < 0.5]
        dots_old = m.dots(old, color=P.RED, opacity=0.3)
        dots_new = m.dots(new, color=P.TEAL)
        dots_left = m.dots(left, color=P.TEAL)
        des = VGroup(*[m.box(b["ra_lo"], b["ra_hi"], b["dec_lo"], b["dec_hi"],
                             color=P.MUTED, opacity=0.18, stroke_width=1.0)
                       for b in d["des_boxes"]])
        key1 = dot_key([(f"ruled out by ZTF + DES  {100 * d['before']:.1f}%", P.RED),
                        (f"not yet searched  {100 * (1 - d['before']):.1f}%", P.TEAL)])
        key1.next_to(m.frame, DOWN, buff=0.62)
        cap = layout.caption(
            f"{len(pop)} predicted Planet Nines, scored by ZTF and DES (grey boxes: DES)",
            font_size=21)
        self.play(FadeIn(m), run_time=0.7)
        self.play(Create(ecl), Create(gal), FadeIn(des), run_time=1.0)
        self.play(FadeIn(dots_old), FadeIn(dots_new), FadeIn(dots_left), FadeIn(key1),
                  FadeIn(cap), run_time=1.4)
        timing.hold_to_read(self, cap, key1, settle=0.8)

        # 2. what Pan-STARRS1 adds
        foot = m.dec_band(ps1["dec_limit_deg"], 90)
        cap2 = layout.caption(
            f"Pan-STARRS1: all sky north of {ps1['dec_limit_deg']:.0f}°, to V ≈ {ps1['depth_v']:.1f}, "
            f"linked if seen on {ps1['linking_threshold']} of {ps1['n_epochs']} epochs",
            font_size=21)
        self.play(FadeIn(foot), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        timing.hold_to_read(self, cap2, settle=0.6)

        # 3. the members only PS1 would have linked
        key3 = dot_key([(f"ZTF + DES  {100 * d['before']:.1f}%", P.MUTED),
                        (f"new with PS1  +{100 * d['ps1_unique']:.1f}%", P.RED),
                        (f"still viable  {100 * d['remaining']:.1f}%", P.TEAL)])
        key3.next_to(m.frame, DOWN, buff=0.62)
        cap3 = layout.caption(
            f"Nothing found. PS1 would have linked {100 * d['ps1_total']:.0f}% of them "
            f"(paper: {100 * pub['ps1_total']:.0f}%); red ones are new",
            font_size=21)
        self.play(dots_new.animate.set_color(P.RED).set_opacity(1.0),
                  dots_old.animate.set_color(P.MUTED).set_opacity(0.45),
                  FadeOut(key1), FadeIn(key3), FadeOut(cap2), FadeIn(cap3), run_time=1.6)
        timing.hold_to_read(self, cap3, key3, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, des, foot, dots_old, dots_new, dots_left, key3,
                                 cap3)))

        # 4. why: PS1 goes a magnitude deeper over the same sky
        v = np.array([s["v_mag"] for s in pop])
        w_old = np.array([s["p_before"] for s in pop])
        w_all = np.array([s["p_after"] for s in pop])
        edges = np.arange(17.0, 26.01, 0.5)
        h_old, _ = np.histogram(v, bins=edges, weights=w_old)
        h_all, _ = np.histogram(v, bins=edges, weights=w_all)
        tot, _ = np.histogram(v, bins=edges)
        top = max(20, int(np.ceil(tot.max() / 20.0) * 20))
        ax, labels = widgets.labeled_axes(
            [17, 26, 1], [0, top, top // 4], x_label="apparent magnitude V  (fainter →)",
            y_label="synthetic planets", y_rotate=True, numbers=True,
            x_length=7.0, y_length=4.0, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.3)
        bars_old = widgets.histogram(ax, edges, h_old, color=P.MUTED, opacity=0.6)
        bars_new = widgets.histogram(ax, edges, h_all - h_old, color=P.RED, opacity=0.8,
                                     base=h_old)
        bars_left = widgets.histogram(ax, edges, tot - h_all, color=P.TEAL, opacity=0.55,
                                      base=h_all)
        ztf_line = widgets.marker_line(ax, d["ztf_depth_v"], (0, top),
                                       f"ZTF {d['ztf_depth_v']:.1f}", color=P.MUTED, side=UP)
        ps1_line = widgets.marker_line(ax, ps1["depth_v"], (0, top),
                                       f"PS1 {ps1['depth_v']:.1f}", side=UP)
        ztf_line[1].shift(LEFT * 0.45)
        ps1_line[1].shift(RIGHT * 0.45)
        legend = VGroup(
            layout.label("ruled out by ZTF + DES", font_size=15, color=P.MUTED),
            layout.label("newly ruled out by PS1", font_size=15, color=P.RED),
            layout.label("still viable", font_size=15, color=P.TEAL),
        ).arrange(DOWN, buff=0.1, aligned_edge=LEFT)
        legend.move_to([4.3, 2.35, 0])
        cap4 = layout.caption("A magnitude deeper than ZTF over the same three quarters of the sky",
                              font_size=21)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(bars_old, lag_ratio=0.1), FadeIn(bars_left, lag_ratio=0.1),
                  Create(ztf_line), run_time=1.2)
        self.play(FadeIn(bars_new, lag_ratio=0.1), Create(ps1_line), FadeIn(legend),
                  FadeIn(cap4), run_time=1.4)
        timing.hold_to_read(self, cap4, legend, settle=0.8)

        # 5. what is left defines the updated planet
        prior, surv = d["prior"], d["survivors"]
        head = layout.label("median of what survives", font_size=17, color=P.FG, weight="BOLD")
        rows = [
            compare_row("a", prior["a_au"], surv["a_au"], pub["a_au"], "AU"),
            compare_row("mass", prior["mass_earth"], surv["mass_earth"], pub["mass_earth"], "M⊕",
                        digits=1),
            compare_row("V", prior["v_mag"], surv["v_mag"], pub["v_mag"], "mag", digits=1),
        ]
        table = VGroup()
        for r, cells in enumerate(rows):
            xs = [2.35, 3.05, 3.6, 4.55, 5.95]
            for x, c in zip(xs, cells):
                c.move_to([x, 0.55 - 0.55 * r, 0])
                table.add(c)
        head.move_to([4.2, 1.2, 0])
        cap5 = layout.caption(
            "The survivors are farther and fainter than the prediction the searches started from",
            font_size=21)
        self.play(FadeIn(head), FadeIn(table, lag_ratio=0.1), FadeOut(cap4), FadeIn(cap5))
        timing.hold_to_read(self, cap5, table, settle=1.0)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, f"Three searches, {100 * d['cumulative']:.0f}% ruled out: what is left is "
                  f"near V ≈ {surv['v_mag']:.0f}.")
