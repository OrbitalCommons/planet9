"""Russell & White (2025) -- the radius, composition, albedo and absolute
magnitude of Planet Nine.

From cold exoplanets of measured mass and radius the paper argues a 6.6
Earth-mass Planet Nine is a mini-Neptune: 2.0-2.6 Earth radii, a thin
hydrogen-helium envelope, albedo 0.33-0.47. The crate folds those endpoints
through the reflected-light law (p9-core photometry) to get the absolute and
apparent magnitudes. Everything is from anim.json -> papers ->
p9-2025-russell-albedo.
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
    DashedVMobject,
    FadeIn,
    FadeOut,
    Polygon,
    Scene,
    VGroup,
)

import p9_manim as P
from p9_manim import dataio, layout, paper, timing, widgets

CRATE = "p9-2025-russell-albedo"


def rows(items, font_size=15, buff=0.14, align=LEFT):
    g = VGroup(*[layout.label(t, font_size=font_size, color=c) for t, c in items])
    return g.arrange(DOWN, buff=buff, aligned_edge=align)


def globe(radius_earth, albedo, scale, color, title, lines, baseline):
    """A disc drawn to scale, brightness following the albedo, with its facts."""
    disc = Circle(radius=radius_earth * scale, color=color, stroke_width=2)
    disc.set_fill(color, opacity=0.15 + 1.2 * albedo)
    head = layout.label(title, font_size=16, color=color, weight="BOLD")
    body = VGroup(*[layout.label(t, font_size=14, color=P.FG) for t in lines])
    body.arrange(DOWN, buff=0.1)
    text = VGroup(head, body).arrange(DOWN, buff=0.14)
    disc.move_to([0, baseline + radius_earth * scale, 0])
    text.move_to([0, baseline - 0.25 - text.height / 2, 0])
    return VGroup(disc, text)


class RussellAlbedo2025(Scene):
    def construct(self):
        d = (dataio.section("papers") or {})[CRATE]
        pub = d["published"]
        faint, bright, nep = d["faint"], d["bright"], d["scaled_neptune"]
        r_star = d["most_likely_distance_au"]

        self.add(paper.scene_header(CRATE))

        # 1. what a 6.6 Earth-mass planet is made of, to scale
        scale = 0.42
        base = -0.6
        earth = globe(1.0, 0.3, scale, P.MUTED, "Earth", ["1 R⊕"], base)
        small = globe(faint["radius_earth"], faint["albedo"], scale, P.BLUE, "mini-Neptune, thin",
                      [f"{faint['radius_earth']:.1f} R⊕",
                       f"{100 * faint['envelope_fraction']:.1f}% hydrogen and helium",
                       f"reflects {100 * faint['albedo']:.0f}%"], base)
        large = globe(bright["radius_earth"], bright["albedo"], scale, P.BLUE,
                      "mini-Neptune, thick",
                      [f"{bright['radius_earth']:.1f} R⊕",
                       f"{100 * bright['envelope_fraction']:.1f}% hydrogen and helium",
                       f"reflects {100 * bright['albedo']:.0f}%"], base)
        old = globe(nep["radius_earth"], nep["albedo"], scale, P.PURPLE,
                    "scaled-down Neptune",
                    [f"{nep['radius_earth']:.1f} R⊕", "what the searches assumed",
                     f"reflects {100 * nep['albedo']:.0f}%"], base)
        for g, x in ((earth, -5.4), (small, -2.2), (large, 1.2), (old, 4.8)):
            g.shift(RIGHT * x)
        cap = layout.caption(
            f"Cold exoplanets of {pub['mass_earth']:.1f} M⊕ are mini-Neptunes: "
            "rock and ice under a thin envelope", font_size=21)
        self.play(FadeIn(earth), FadeIn(cap), run_time=0.8)
        self.play(FadeIn(small), FadeIn(large), run_time=1.2)
        timing.hold_to_read(self, cap, small[1], large[1], settle=0.6)
        cap_b = layout.caption(
            f"Smaller than the {nep['radius_earth']:.1f} R⊕ scaled-down Neptune, so it "
            "reflects less sunlight", font_size=21)
        self.play(FadeIn(old), FadeOut(cap), FadeIn(cap_b), run_time=1.0)
        timing.hold_to_read(self, cap_b, old[1], settle=0.8)
        self.play(FadeOut(VGroup(earth, small, large, old, cap_b)))

        # 2. how bright that makes it
        dist = np.array(d["distance_au"])
        ax, labels = widgets.labeled_axes(
            [300, 1000, 100], [18, 25, 1], x_label="heliocentric distance (AU)",
            y_label="apparent magnitude V  (fainter ↑)", y_rotate=True, numbers=True,
            x_length=7.4, y_length=4.2, shift_down=-0.3)
        VGroup(ax, labels).shift(LEFT * 2.3)
        band = Polygon(
            *[ax.c2p(x, y) for x, y in zip(dist, faint["v_curve"])],
            *[ax.c2p(x, y) for x, y in zip(dist[::-1], bright["v_curve"][::-1])],
            stroke_width=0).set_fill(P.BLUE, opacity=0.3)
        c_faint = widgets.curve(ax, dist, faint["v_curve"], color=P.BLUE, stroke_width=2.5)
        c_bright = widgets.curve(ax, dist, bright["v_curve"], color=P.BLUE, stroke_width=2.5)
        c_old = DashedVMobject(
            widgets.curve(ax, dist, nep["v_curve"], color=P.PURPLE, stroke_width=2.5),
            num_dashes=40)
        depths = VGroup()
        for name, mag in (("Pan-STARRS1", d["depths"]["ps1"]), ("DES", d["depths"]["des"])):
            line = DashedLine(ax.c2p(300, mag), ax.c2p(1000, mag), color=P.MUTED,
                              stroke_width=1.6)
            lab = layout.label(f"{name} depth {mag:.1f}", font_size=13, color=P.MUTED)
            if name == "DES":
                lab.next_to(ax.c2p(300, mag), UP, buff=0.05).shift(RIGHT * 0.95)
            else:
                lab.next_to(ax.c2p(1000, mag), DOWN, buff=0.06).shift(LEFT * 1.1)
            depths.add(VGroup(line, lab))
        key = rows([
            ("mini-Neptune, thick to thin envelope", P.BLUE),
            ("scaled-down Neptune", P.PURPLE),
        ])
        key.move_to([4.55, 2.4, 0])
        cap2 = layout.caption("Reflected sunlight fades as the fourth power of distance",
                              font_size=21)
        self.play(Create(ax), FadeIn(labels), run_time=0.9)
        self.play(FadeIn(depths), Create(c_old), FadeIn(key[1]), FadeIn(cap2), run_time=1.0)
        self.play(FadeIn(band), Create(c_faint), Create(c_bright), FadeIn(key[0]),
                  run_time=1.2)
        timing.hold_to_read(self, cap2, key, settle=0.5)

        here = DashedLine(ax.c2p(r_star, 18), ax.c2p(r_star, 25), color=P.FG, stroke_width=1.6)
        span = VGroup(*[
            widgets.curve(ax, [r_star, r_star], [bright["v"], faint["v"]], color=P.TEAL,
                          stroke_width=6)])
        facts = rows([
            ("absolute magnitude H", P.FG),
            (f"{d['h_bright']:.1f} to {d['h_faint']:.1f}", P.TEAL),
            (f"paper: {pub['h'][0]:.1f} to {pub['h'][1]:.1f}", P.MUTED),
            (f"apparent magnitude at {r_star:.0f} AU", P.FG),
            (f"V = {d['v_bright']:.1f} to {d['v_faint']:.1f}", P.TEAL),
            (f"paper: {pub['v'][0]:.1f} to {pub['v'][1]:.1f}", P.MUTED),
            ("disk seen from Earth", P.FG),
            (f"{faint['disk_arcsec']:.2f}″ to {bright['disk_arcsec']:.2f}″ across", P.TEAL),
        ], font_size=15, buff=0.1)
        facts[3:].shift(DOWN * 0.2)
        facts[6:].shift(DOWN * 0.2)
        facts.next_to(key, DOWN, buff=0.4, aligned_edge=LEFT)
        cap3 = layout.caption(
            f"No distance is published: {r_star:.0f} AU is solved here to match the "
            "paper's magnitudes", font_size=21)
        self.play(Create(here), Create(span), FadeIn(facts, lag_ratio=0.1), FadeOut(cap2),
                  FadeIn(cap3), run_time=1.4)
        timing.hold_to_read(self, cap3, facts, settle=1.0)
        self.play(FadeOut(cap3))

        layout.show_takeaway(
            self, "A mini-Neptune Planet Nine is smaller than assumed and shines near "
                  f"V = {0.5 * (d['v_bright'] + d['v_faint']):.0f}.")
