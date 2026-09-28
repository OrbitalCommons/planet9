"""State-of-the-field gauges: compact panels that show one quantity the whole
literature is arguing about, *before* and *after* a given paper.

Four quantities run through the film:

* **clustering** -- how significant is the apsidal alignment (sigma), a see-saw
  between the discovery camp and the bias critiques;
* **orbit** -- the best-fit Planet Nine (mass, semi-major axis, eccentricity);
* **excluded** -- the fraction of the predicted orbits direct searches have
  ruled out;
* **sample** -- how many distant objects the inference rests on.

Every gauge is built from ledger payloads and returns ``(panel, animate)``:
``panel`` is the static "before" state, ``animate(scene)`` plays the change.
All gauges fit a 6.0 x 3.3 box centred on the origin; the caller positions it.
"""
import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Arrow,
    Create,
    Dot,
    FadeIn,
    GrowFromCenter,
    Line,
    Rectangle,
    Transform,
    VGroup,
)

from . import layout, ledger, orbits
from . import theme as T

W, H = 6.0, 3.3


def _frame(title):
    box = Rectangle(width=W, height=H, color=T.MUTED, stroke_width=1.2)
    box.set_stroke(opacity=0.55).set_fill("#1f2030", opacity=0.85)
    cap = layout.label(title, font_size=15, color=T.MUTED)
    cap.move_to(box.get_top() + DOWN * 0.24)
    return VGroup(box, cap)


def _local(frame, mob):
    """Map ``mob`` from panel-local coordinates (origin-centred, unscaled) into
    the panel's current placement, which the caller may have moved and scaled."""
    return mob.scale(frame[0].width / W, about_point=np.zeros(3)).shift(frame[0].get_center())


# ---- excluded ---------------------------------------------------------------

def excluded(entry):
    """Stacked bar of the searched fraction; earlier surveys muted, this one in
    colour, the remainder teal."""
    prior = ledger.earlier(entry["crate"], "excluded")
    pay = entry["excluded"]
    after = float(pay["cumulative"])
    before = float(prior[-1]["excluded"]["cumulative"]) if prior else 0.0
    after = max(after, before)

    frame = _frame("predicted orbits ruled out by direct searches")
    bar_w, bar_h = W - 0.9, 0.62
    x0 = -bar_w / 2
    track = Rectangle(width=bar_w, height=bar_h, color=T.MUTED, stroke_width=1.5)
    track.set_fill(T.TEAL, opacity=0.10)

    segs = VGroup()
    labels = VGroup()
    cum = 0.0
    for e in prior:
        c = float(e["excluded"]["cumulative"])
        if c <= cum:
            continue
        seg = Rectangle(width=bar_w * (c - cum), height=bar_h, stroke_width=0)
        seg.set_fill(T.RED, opacity=0.30)
        seg.move_to([x0 + bar_w * (cum + c) / 2, 0, 0])
        segs.add(seg)
        lab = layout.label(str(e["excluded"].get("label", "")), font_size=12, color=T.MUTED)
        lab.next_to(seg, DOWN, buff=0.1)
        if lab.width < seg.width + 0.35:
            labels.add(lab)
        cum = c
    before_txt = layout.label(f"{100 * before:.0f}%", font_size=30, color=T.FG, weight="BOLD")
    before_txt.move_to([0, 0.95, 0])
    remain = layout.label("still viable", font_size=12, color=T.TEAL)
    remain.move_to([x0 + bar_w * (1 + before) / 2, 0, 0])
    panel = VGroup(frame, track, segs, labels, before_txt, remain)

    def animate(scene):
        new = Rectangle(width=max(bar_w * (after - before), 1e-3), height=bar_h, stroke_width=0)
        new.set_fill(T.RED, opacity=0.78)
        _local(frame, new.move_to([x0 + bar_w * (before + after) / 2, 0, 0]))
        lab = layout.label(f"{pay.get('label', '')}  +{100 * (after - before):.1f}%",
                           font_size=14, color=T.RED, weight="BOLD")
        lab.next_to(new, UP, buff=0.12)
        after_txt = layout.label(f"{100 * after:.0f}%", font_size=30, color=T.RED, weight="BOLD")
        after_txt.move_to(before_txt)
        remain2 = remain.copy().move_to(
            _local(frame, Dot([x0 + bar_w * (1 + after) / 2, 0, 0])).get_center())
        if frame[0].width * (1 - after) * bar_w / W < remain.width + 0.2:
            remain2.set_opacity(0)
        scene.play(GrowFromCenter(new), FadeIn(lab, shift=UP * 0.1),
                   Transform(before_txt, after_txt), Transform(remain, remain2), run_time=1.6)
        panel.add(new, lab)

    return panel, animate


# ---- clustering -------------------------------------------------------------

SIGMA_MAX = 5.0


def clustering(entry):
    """A 0-5 sigma scale carrying every earlier significance claim as a dot
    (green = finds clustering, red = finds it consistent with bias); the new
    claim drops in with an arrow from the previous one."""
    prior = ledger.earlier(entry["crate"], "clustering")
    pay = entry["clustering"]
    frame = _frame("significance of the orbital clustering")
    ax_w = W - 1.2
    x0 = -ax_w / 2
    axis = Line([x0, -0.55, 0], [x0 + ax_w, -0.55, 0], color=T.MUTED, stroke_width=2)
    ticks = VGroup()
    for s in range(0, int(SIGMA_MAX) + 1):
        x = x0 + ax_w * s / SIGMA_MAX
        ticks.add(Line([x, -0.63, 0], [x, -0.47, 0], color=T.MUTED, stroke_width=1.5))
        ticks.add(layout.label(f"{s}σ", font_size=13, color=T.MUTED).move_to([x, -0.9, 0]))
    chance = layout.label("consistent with chance", font_size=11, color=T.MUTED)
    chance.move_to([x0 + ax_w * 0.12, -1.25, 0])
    strong = layout.label("hard to explain by chance", font_size=11, color=T.MUTED)
    strong.move_to([x0 + ax_w * 0.82, -1.25, 0])

    def pos(sig, level):
        x = x0 + ax_w * min(max(float(sig), 0.0), SIGMA_MAX) / SIGMA_MAX
        return np.array([x, -0.55 + 0.34 + 0.36 * level, 0.0])

    def colour(stance):
        return T.GREEN if stance == "pro" else T.RED

    dots = VGroup()
    recent = prior[-3:]
    for k, e in enumerate(recent):
        p = e["clustering"]
        d = Dot(pos(p["sigma"], k), radius=0.07, color=colour(p.get("stance", "pro")))
        d.set_opacity(0.55)
        lab = layout.label(str(p.get("label", "")), font_size=11, color=T.MUTED)
        lab.next_to(d, RIGHT, buff=0.08)
        dots.add(VGroup(d, lab))
    panel = VGroup(frame, axis, ticks, chance, strong, dots)

    def animate(scene):
        level = len(recent)
        col = colour(pay.get("stance", "pro"))
        d = _local(frame, Dot(pos(pay["sigma"], level), radius=0.11, color=col))
        here = d.get_center()
        lab = layout.label(f"{pay.get('label', '')}  {float(pay['sigma']):.1f}σ",
                           font_size=14, color=col, weight="BOLD")
        lab.next_to(d, RIGHT, buff=0.1)
        if lab.get_right()[0] > frame.get_right()[0] - 0.1:
            lab.next_to(d, LEFT, buff=0.1)
        anims = [GrowFromCenter(d), FadeIn(lab)]
        if recent:
            prev = dots[-1][0].get_center()
            if np.linalg.norm(here - prev) > 0.35:
                anims.append(Create(Arrow(prev, here, buff=0.12, color=col,
                                          stroke_width=2.5, max_tip_length_to_length_ratio=0.12)))
        scene.play(*anims, run_time=1.4)
        panel.add(d, lab)

    return panel, animate


# ---- orbit ------------------------------------------------------------------

def _orbit_curve(pay, scale, color, opacity=1.0, width=3.0):
    a = float(pay["a"]) * scale
    e = float(pay.get("e", 0.3))
    return orbits.ellipse_orbit(a, e, color=color, varpi=np.pi, stroke_width=width,
                                opacity=opacity)


def _orbit_text(pay):
    bits = [f"{float(pay['mass']):.1f} M⊕", f"a = {float(pay['a']):.0f} AU"]
    if pay.get("e") is not None:
        bits.append(f"e = {float(pay['e']):.2f}")
    if pay.get("i") is not None:
        bits.append(f"i = {float(pay['i']):.0f}°")
    return "   ".join(bits)


def orbit(entry):
    """Top-down view: the previous best-fit orbit (muted) against the new one,
    Neptune's orbit for scale."""
    prior = ledger.earlier(entry["crate"], "orbit")
    pay = entry["orbit"]
    prev = prior[-1]["orbit"] if prior else None
    frame = _frame("best-fit Planet Nine orbit (top view)")
    biggest = max([float(pay["a"]) * (1 + float(pay.get("e", 0.3)))]
                  + ([float(prev["a"]) * (1 + float(prev.get("e", 0.3)))] if prev else []))
    scale = (W * 0.62) / (2.0 * biggest)
    centre = np.array([0.95, 0.15, 0.0])
    sun = Dot(np.zeros(3), radius=0.05, color=T.SUN)
    nep = orbits.ellipse_orbit(30.0 * scale, 0.0, color=T.MUTED, stroke_width=1.2)
    art = VGroup(sun, nep)
    old = None
    old_txt = None
    if prev:
        old = _orbit_curve(prev, scale, T.MUTED, opacity=0.9, width=2.0)
        art.add(old)
        old_txt = layout.label("before:  " + _orbit_text(prev), font_size=13, color=T.MUTED)
    art.shift(centre)
    panel = VGroup(frame, art)
    if old_txt is not None:
        old_txt.move_to([0, -H / 2 + 0.62, 0])
        panel.add(old_txt)

    def animate(scene):
        new = _local(frame, _orbit_curve(pay, scale, T.BLUE).shift(centre))
        txt = layout.label("after:   " + _orbit_text(pay), font_size=14, color=T.BLUE,
                           weight="BOLD")
        _local(frame, txt.move_to([0, -H / 2 + 0.3, 0]))
        scene.play(Create(new), FadeIn(txt), run_time=1.8)
        panel.add(new, txt)

    return panel, animate


# ---- sample -----------------------------------------------------------------

def sample(entry):
    """Dots for the distant objects the inference rests on: before (muted) and
    the ones this paper adds (green)."""
    prior = ledger.earlier(entry["crate"], "sample")
    pay = entry["sample"]
    before = int(prior[-1]["sample"]["n"]) if prior else 0
    after = int(pay["n"])
    frame = _frame(str(pay.get("what", "distant objects in the sample")))
    per_row = 12
    n_show = max(before, after)

    def spot(k):
        row, col = divmod(k, per_row)
        return np.array([-W / 2 + 0.75 + col * 0.4, 0.55 - row * 0.4, 0.0])

    dots = VGroup(*[Dot(spot(k), radius=0.09, color=T.MUTED).set_opacity(0.8)
                    for k in range(min(before, n_show))])
    count = layout.label(f"N = {before}", font_size=26, color=T.FG, weight="BOLD")
    count.move_to([0, -H / 2 + 0.5, 0])
    panel = VGroup(frame, dots, count)

    def animate(scene):
        anims = []
        if after >= before:
            new = VGroup(*[_local(frame, Dot(spot(k), radius=0.09, color=T.GREEN))
                           for k in range(before, after)])
            if len(new):
                anims.append(FadeIn(new, lag_ratio=0.15))
                panel.add(new)
        else:
            drop = VGroup(*dots[after:])
            anims.append(drop.animate.set_opacity(0.15))
        col = T.GREEN if after >= before else T.RED
        count2 = layout.label(f"N = {after}", font_size=26, color=col, weight="BOLD")
        count2.move_to(count)
        anims.append(Transform(count, count2))
        scene.play(*anims, run_time=1.5)

    return panel, animate


BUILDERS = {
    "excluded": excluded,
    "clustering": clustering,
    "orbit": orbit,
    "sample": sample,
}


def build(entry):
    """The gauge for the entry's track, or ``(None, None)`` when it moves none."""
    track = entry.get("track")
    if track in BUILDERS and entry.get(track):
        pay = entry[track]
        if any(v is None for v in pay.values()):
            return None, None
        return BUILDERS[track](entry)
    return None, None
