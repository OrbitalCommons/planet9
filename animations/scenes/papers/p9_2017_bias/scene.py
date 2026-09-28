"""Brown (2017) -- observational bias and the clustering of distant eccentric
Kuiper belt objects.

Distant eccentric objects are bright enough to find only near perihelion, so
where the surveys looked decides which perihelion directions can be found at
all. The paper models that bias and asks how often it alone aligns ten orbits
as well as the real ones. The sky positions, the bias curve, the Monte Carlo
distributions and the probabilities are the crate's own (anim.json -> papers ->
p9-2017-bias); the crate's bias is a simplified stand-in for the paper's
survey-by-survey model, so its numbers sit beside the published ones.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Create,
    DashedLine,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    Rectangle,
    Scene,
    VGroup,
    VMobject,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, sky, timing, widgets

CRATE = "p9-2017-bias"

# Published in the paper's abstract; drawn only as labelled comparisons.
PAPER_P_VARPI = 0.012
PAPER_P_COMBINED = 0.00025


def plated(text, color, font_size=14, fill="#16171f"):
    lab = layout.label(text, font_size=font_size, color=color)
    box = Rectangle(width=lab.width + 0.1, height=lab.height + 0.1, stroke_width=0)
    box.set_fill(fill, opacity=0.9).move_to(lab)
    return VGroup(box, lab)


def step_outline(ax, edges, counts, color, stroke_width=2.5):
    """A histogram drawn as a stepped outline, so it reads through a filled one."""
    pts = [ax.c2p(edges[0], 0)]
    for k, n in enumerate(counts):
        pts += [ax.c2p(edges[k], n), ax.c2p(edges[k + 1], n)]
    pts.append(ax.c2p(edges[-1], 0))
    m = VMobject(color=color, stroke_width=stroke_width)
    m.set_points_as_corners(pts)
    return m


def swatch_row(swatch, text, color, font_size=15):
    return VGroup(swatch, layout.label(text, font_size=font_size, color=color)).arrange(
        RIGHT, buff=0.15)


class Bias2017(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        objs = d["objects"]
        bias = d["bias"]
        chance = d["chance"]

        self.add(paper.scene_header(CRATE))

        # 1. where the ten objects were when they could be found
        m = sky.SkyMap(width=12.0, dec_range=(-60, 60), centre=(0.0, 0.45, 0.0))
        ecl, gal = m.reference_curves()
        spots = VGroup(*[Dot(m.p(o["peri_ra_deg"], o["peri_dec_deg"]), radius=0.07,
                             color=P.GREEN) for o in objs])
        key = m.legend([("ecliptic", P.ORANGE), ("galactic plane ±10°", P.PURPLE),
                        ("perihelion of a distant object", P.GREEN)])
        key.next_to(m.frame, DOWN, buff=0.75)
        cap = layout.caption(
            f"The {d['n_sample']} known orbits beyond 230 AU: each dot is where one comes to "
            f"perihelion", font_size=22)
        self.play(FadeIn(m), run_time=0.8)
        self.play(Create(ecl), Create(gal), FadeIn(key), run_time=1.0)
        self.play(FadeIn(spots, lag_ratio=0.15), FadeIn(cap), run_time=1.6)
        timing.hold_to_read(self, cap, settle=1.0)
        cap1b = layout.caption(
            "They are found only near perihelion, and only where surveys have looked",
            font_size=22)
        self.play(FadeOut(cap), run_time=0.4)
        self.play(FadeIn(cap1b), run_time=0.6)
        timing.hold_to_read(self, cap1b, settle=1.0)
        self.play(FadeOut(VGroup(m, ecl, gal, spots, key, cap1b)))

        # 2. the modelled bias against longitude of perihelion
        ax, labels = widgets.labeled_axes(
            [0, 360, 45], [0, 1.2, 0.2], x_label="longitude of perihelion (degrees)",
            y_label="relative chance of discovery", y_rotate=True, numbers=True,
            x_length=10.5, y_length=4.3, shift_down=-0.3)
        curve = widgets.curve(ax, bias["lon_deg"], bias["weight"], color=P.PURPLE)
        dips = sorted(zip(bias["weight"], bias["lon_deg"]))
        deep = dips[0][1]
        shallow = next(lon for w, lon in dips if abs(lon - deep) > 90)
        dip_labs = VGroup(
            plated("toward the galactic centre", P.PURPLE, fill=P.BG)
            .next_to(ax.c2p(deep, 0.03), UP, buff=0.35),
            plated("Milky Way, anticentre", P.PURPLE, fill=P.BG)
            .next_to(ax.c2p(shallow, 0.62), DOWN, buff=0.45),
        )
        cap2 = layout.caption(
            "The modelled bias: perihelia in the Milky Way are hard to find",
            font_size=22)
        self.play(Create(ax), FadeIn(labels))
        self.play(Create(curve), FadeIn(cap2), run_time=1.6)
        self.play(FadeIn(dip_labs))
        timing.hold_to_read(self, cap2, settle=0.8)

        marks = VGroup(*[
            Line(ax.c2p(o["varpi_deg"], 1.02), ax.c2p(o["varpi_deg"], 1.14), color=P.GREEN,
                 stroke_width=4) for o in objs])
        marks_lab = layout.label("the ten real orbits", font_size=15, color=P.GREEN)
        marks_lab.next_to(ax.c2p(d["mean_varpi_deg"], 1.16), UP, buff=0.1)
        inside = sum(1 for o in objs
                     if abs((o["varpi_deg"] - d["mean_varpi_deg"] + 180) % 360 - 180) < 90)
        cap3 = layout.caption(
            f"{inside} of the ten point within 90° of the same direction, "
            f"longitude {d['mean_varpi_deg']:.0f}°", font_size=22)
        self.play(FadeOut(cap2), run_time=0.4)
        self.play(FadeIn(marks, lag_ratio=0.1), FadeIn(marks_lab), FadeIn(cap3))
        timing.hold_to_read(self, cap3, settle=1.0)
        self.play(FadeOut(VGroup(ax, labels, curve, dip_labs, marks, marks_lab, cap3)))

        # 3. how often the bias alone aligns ten orbits this well
        edges = np.array(chance["edges"])
        flat = 100 * np.array(chance["uniform"])
        skew = 100 * np.array(chance["biased"])
        top = float(np.ceil(max(flat.max(), skew.max())))
        ax2, labels2 = widgets.labeled_axes(
            [0, 1, 0.2], [0, top, 2], x_label="how well ten orbits line up  (0 = not at all, "
                                                "1 = perfectly)",
            y_label="share of random samples (%)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=-0.3)
        plot = VGroup(ax2, labels2).shift(LEFT * 2.6)
        bars_skew = widgets.histogram(ax2, edges, skew, color=P.PURPLE, opacity=0.45)
        bars_flat = step_outline(ax2, edges, flat, color=P.FG)
        seen = DashedLine(ax2.c2p(d["r_bar"], 0), ax2.c2p(d["r_bar"], top), color=P.GREEN,
                          stroke_width=3)
        seen_lab = layout.label("the ten real orbits", font_size=15, color=P.GREEN)
        seen_lab.next_to(seen.get_end(), UP, buff=0.08)
        legend = VGroup(
            swatch_row(Line(LEFT * 0.2, RIGHT * 0.2, color=P.FG, stroke_width=2.5),
                       "no bias: any direction equally likely", P.FG),
            swatch_row(Rectangle(width=0.4, height=0.2, stroke_width=0)
                       .set_fill(P.PURPLE, opacity=0.6),
                       "with the modelled survey bias", P.PURPLE),
        ).arrange(DOWN, buff=0.16, aligned_edge=LEFT)
        legend.move_to([2.65, 2.3, 0], aligned_edge=LEFT)
        rows = VGroup(
            layout.label("chance of lining up this well", font_size=16, color=P.MUTED),
            layout.label(f"no bias:  {100 * chance['p_uniform']:.1f}%", font_size=19),
            layout.label(f"with bias:  {100 * d['p_varpi']:.1f}%", font_size=21,
                         color=P.PURPLE, weight="BOLD"),
            layout.label(f"paper, its full survey model:  {100 * PAPER_P_VARPI:.1f}%",
                         font_size=17, color=P.MUTED),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        rows.move_to([4.2, 0.6, 0], aligned_edge=LEFT).shift(LEFT * 1.55)
        cap4 = layout.caption(
            "Bias makes a chance alignment likelier, but still rare", font_size=22)
        self.play(Create(ax2), FadeIn(labels2))
        self.play(FadeIn(bars_flat, lag_ratio=0.05), FadeIn(bars_skew, lag_ratio=0.05),
                  FadeIn(legend), run_time=1.4)
        self.play(Create(seen), FadeIn(seen_lab), FadeIn(cap4))
        self.play(FadeIn(rows, lag_ratio=0.2), run_time=1.2)
        timing.hold_to_read(self, cap4, rows, settle=1.2)

        both = VGroup(
            layout.label("and a second alignment at once", font_size=16, color=P.MUTED),
            layout.label(f"here, with ω:  {100 * d['p_combined']:.3f}%", font_size=21,
                         color=P.GREEN, weight="BOLD"),
            layout.label(f"paper, with orbital poles:  {100 * PAPER_P_COMBINED:.3f}%",
                         font_size=17, color=P.MUTED),
        ).arrange(DOWN, buff=0.2, aligned_edge=LEFT)
        both.next_to(rows, DOWN, buff=0.45, aligned_edge=LEFT)
        cap5 = layout.caption(
            "Ask for a second, independent alignment as well and chance all but vanishes",
            font_size=22)
        self.play(FadeOut(cap4), run_time=0.4)
        self.play(FadeIn(both, lag_ratio=0.2), FadeIn(cap5))
        timing.hold_to_read(self, cap5, both, settle=1.2)
        self.play(FadeOut(cap5))

        layout.show_takeaway(
            self, "The bias is real, but it does not produce an alignment this strong.")
