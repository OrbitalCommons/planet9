"""Brown & Batygin (2021) -- a search for Planet Nine in the ZTF archive.

The first search scored against the whole predicted population: every synthetic
Planet Nine drawn from the orbit posterior is placed on the sky with its
brightness, and the ZTF survey model decides whether three years of public
images would have linked it. Nothing was found, so every member ZTF would have
caught is ruled out -- the bright, nearby ones. Reproduced in p9-2021-ztf;
the population, detection probabilities and efficiency curve shown here are the
crate's own (anim.json -> papers -> p9-2021-ztf).
"""
import numpy as np
from manim import (
    DOWN,
    UP,
    Create,
    FadeIn,
    FadeOut,
    LaggedStart,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2021-ztf"


class Ztf2021(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pop = d["population"]
        frac = d["ztf"]
        dec_limit = d["dec_limit_deg"]

        self.add(paper.scene_header(CRATE))

        # 1. where the predicted planet could be tonight
        m = sky.SkyMap(width=12.0, dec_range=(-75, 75), centre=(0.0, 0.1, 0.0))
        ecl, gal = m.reference_curves()
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), run_time=1.2)
        key = m.legend([("ecliptic", P.ORANGE), ("galactic plane ±10°", P.PURPLE)])
        self.play(FadeIn(key))

        caught = [s for s in pop if s["p_detect"] >= 0.5]
        missed = [s for s in pop if s["p_detect"] < 0.5]
        dots_c = m.dots(caught, color=P.TEAL)
        dots_m = m.dots(missed, color=P.TEAL)
        cap = layout.caption(f"{len(pop)} synthetic Planet Nines drawn from the predicted orbits",
                             font_size=22)
        self.play(LaggedStart(*[FadeIn(x) for x in (dots_c, dots_m)], lag_ratio=0.2),
                  FadeIn(cap), FadeOut(key), run_time=1.6)
        timing.hold_to_read(self, cap, settle=0.8)

        # 2. what ZTF can see
        foot = m.dec_band(dec_limit, 90)
        cap2 = layout.caption(
            f"ZTF images everything north of {dec_limit:.0f}°, every few nights, to V ≈ 20.5",
            font_size=22)
        self.play(FadeIn(foot), FadeOut(cap), FadeIn(cap2), run_time=1.0)
        timing.hold_to_read(self, cap2, settle=0.6)

        # 3. the ones it would have linked are ruled out
        cap3 = layout.caption(
            "Nothing was found: every one ZTF would have linked is ruled out", font_size=22)
        self.play(dots_c.animate.set_color(P.RED), dots_m.animate.set_opacity(0.55),
                  FadeOut(cap2), FadeIn(cap3), run_time=1.6)
        tally = paper.result_readout("predicted orbits ruled out", f"{100 * frac:.1f}%",
                                     color=P.RED).scale(0.62)
        tally.move_to(m.frame.get_corner(UP + P.layout.RIGHT) + np.array([-1.35, -0.6, 0]))
        self.play(FadeIn(tally))
        timing.hold_to_read(self, cap3, tally, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, foot, dots_c, dots_m, tally, cap3)))

        # 4. why the survivors survive: they are faint
        v = np.array([s["v_mag"] for s in pop])
        pdet = np.array([s["p_detect"] for s in pop])
        edges = np.arange(17.0, 26.01, 0.5)
        hit, _ = np.histogram(v, bins=edges, weights=pdet)
        tot, _ = np.histogram(v, bins=edges)
        top = max(10, int(np.ceil(tot.max() / 20.0) * 20))
        ax, labels = widgets.labeled_axes(
            [17, 26, 1], [0, top, top // 4], x_label="apparent magnitude V  (fainter →)",
            y_label="synthetic planets", y_rotate=True, numbers=True,
            x_length=10.0, y_length=4.2, shift_down=-0.35)
        bars_hit = widgets.histogram(ax, edges, hit, color=P.RED, opacity=0.75)
        bars_left = widgets.histogram(ax, edges, tot - hit, color=P.TEAL, opacity=0.55, base=hit)
        depth = widgets.marker_line(ax, 20.5, (0, top), "ZTF depth  V ≈ 20.5", side=UP)
        self.play(Create(ax), FadeIn(labels))
        self.play(FadeIn(bars_hit, lag_ratio=0.1), FadeIn(bars_left, lag_ratio=0.1), run_time=1.4)
        self.play(Create(depth))
        legend = VGroup(
            layout.label("ruled out by ZTF", font_size=16, color=P.RED),
            layout.label("still hidden: too faint or outside the footprint",
                         font_size=16, color=P.TEAL),
        ).arrange(DOWN, buff=0.12, aligned_edge=P.layout.LEFT)
        # Upper right, over the thin faint tail, clear of the tall bars.
        legend.next_to(ax.c2p(26, top), DOWN + P.layout.LEFT, buff=0.1)
        self.play(FadeIn(legend))
        timing.hold_to_read(self, legend, settle=1.2)

        layout.show_takeaway(
            self, "ZTF removes the bright half of the prediction; what is left is fainter than V ≈ 21.")
