"""The two beats that frame every paper scene.

* :func:`present_claim` -- before the scene: where we are in time, who is
  speaking, and what they set out to show.
* :func:`present_contribution` -- after the scene: what the field knew the day
  after that it did not know the day before, shown on the gauge the paper
  moves, with the workspace's reproduced number set against the published one.
"""
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Dot,
    FadeIn,
    FadeOut,
    Line,
    RoundedRectangle,
    SurroundingRectangle,
    VGroup,
)

from . import gauges, layout, ledger, timing
from . import theme as T

RIBBON_W = 12.4
RIBBON_Y = 3.2


def timeline(entry):
    """The publication ribbon: every ledger paper as a tick coloured by kind,
    past ticks bright, future ticks dim, this paper marked and dated."""
    y0, y1 = ledger.span()
    x0 = -RIBBON_W / 2

    def x_of(t):
        return x0 + RIBBON_W * (t - y0) / (y1 - y0)

    base = Line([x0, RIBBON_Y, 0], [x0 + RIBBON_W, RIBBON_Y, 0], color=T.MUTED, stroke_width=1.5)
    years = VGroup()
    for y in range(y0, y1 + 1, 2):
        years.add(Line([x_of(y), RIBBON_Y - 0.07, 0], [x_of(y), RIBBON_Y + 0.07, 0],
                       color=T.MUTED, stroke_width=1.2))
        years.add(layout.label(str(y), font_size=11, color=T.MUTED)
                  .move_to([x_of(y), RIBBON_Y - 0.26, 0]))
    ticks = VGroup()
    marker = None
    stack = {}
    for e in ledger.entries():
        slot = round(e["t"] * 6)
        level = stack.get(slot, 0)
        stack[slot] = level + 1
        _, col = ledger.kind_label(e)
        x = x_of(e["t"])
        yb = RIBBON_Y + 0.08 + 0.13 * level
        tick = Line([x, yb, 0], [x, yb + 0.1, 0], color=col, stroke_width=2.2)
        if e["crate"] == entry["crate"]:
            marker = (x, col)
            tick.set_stroke(width=0)
        elif e["t"] > entry["t"]:
            tick.set_stroke(opacity=0.22)
        else:
            tick.set_stroke(opacity=0.9)
        ticks.add(tick)
    g = VGroup(base, years, ticks)
    if marker:
        x, col = marker
        dot = Dot([x, RIBBON_Y, 0], radius=0.09, color=col).set_z_index(3)
        when = layout.label(ledger.date_label(entry), font_size=14, color=col, weight="BOLD")
        when.move_to([min(max(x, x0 + 0.5), x0 + RIBBON_W - 0.5), RIBBON_Y - 0.55, 0])
        g.add(dot, when)
    return g


def kind_chip(entry):
    text, col = ledger.kind_label(entry)
    lab = layout.label(text, font_size=13, color=col, weight="BOLD")
    box = RoundedRectangle(width=lab.width + 0.3, height=lab.height + 0.2,
                           corner_radius=0.08, color=col, stroke_width=1.3)
    box.set_fill(col, opacity=0.10)
    lab.move_to(box)
    return VGroup(box, lab)


def header(entry):
    """Kind chip + citation line, left-aligned under the ribbon."""
    chip = kind_chip(entry)
    cite = layout.label(f"{entry['cite']}   arXiv:{entry['arxiv']}", font_size=16, color=T.MUTED)
    row = VGroup(chip, cite).arrange(RIGHT, buff=0.25)
    row.move_to([0, 2.2, 0]).to_edge(LEFT, buff=0.9)
    return row


def _body(text, width, font_size, color=None, weight="NORMAL"):
    return layout.label(timing.wrap(text, width=width), font_size=font_size,
                        color=color or T.FG, weight=weight, line_spacing=0.9)


def result_row(entry):
    """'reproduced here X  |  paper Y' with a tick when they agree to 5%."""
    r = entry.get("result")
    if not r:
        return None
    fmt = r.get("fmt")
    rep = ledger.format_value(r.get("reproduced"), fmt)
    pub = ledger.format_value(r.get("published"), fmt) if r.get("published") is not None else None
    if rep is None and pub is None:
        return None
    parts = [layout.label(str(r["label"]) + ":", font_size=17, color=T.MUTED)]
    if rep is not None:
        parts.append(layout.label("reproduced here", font_size=14, color=T.MUTED))
        parts.append(layout.label(rep, font_size=22, color=T.GREEN, weight="BOLD"))
    if pub is not None:
        parts.append(layout.label("paper", font_size=14, color=T.MUTED))
        parts.append(layout.label(pub, font_size=22, color=T.FG, weight="BOLD"))
    if rep is not None and r.get("published") not in (None, 0):
        rel = abs(float(r["reproduced"]) - float(r["published"])) / abs(float(r["published"]))
        if rel <= float(r.get("tolerance", 0.05)):
            parts.append(layout.label("✓", font_size=22, color=T.GREEN, weight="BOLD"))
    row = VGroup(*parts).arrange(RIGHT, buff=0.22)
    box = SurroundingRectangle(row, color=T.MUTED, buff=0.18, corner_radius=0.1)
    box.set_stroke(opacity=0.5).set_fill("#1f2030", opacity=0.85)
    return VGroup(box, row)


def present_claim(scene, entry):
    """Before the scene: the date on the ribbon, the citation, the headline and
    the claim under test."""
    rib = timeline(entry)
    head = header(entry)
    title = _body(entry["headline"], 40, 38, weight="BOLD")
    title.move_to([0, 0.75, 0])
    claim = _body(entry["claim"], 62, 24, color=T.FG)
    claim.next_to(title, DOWN, buff=0.55)
    scene.play(FadeIn(rib), run_time=0.8)
    scene.play(FadeIn(head, shift=RIGHT * 0.15), FadeIn(title, shift=UP * 0.15), run_time=0.8)
    scene.play(FadeIn(claim, shift=UP * 0.1), run_time=0.6)
    timing.hold_to_read(scene, title, claim, settle=1.2)
    scene.play(FadeOut(VGroup(rib, head, title, claim)), run_time=0.6)


def present_contribution(scene, entry):
    """After the scene: what the paper added, on the gauge it moves."""
    rib = timeline(entry)
    head = header(entry)
    adds_h = layout.label("What this paper added", font_size=20, color=T.TEAL, weight="BOLD")
    panel, animate = gauges.build(entry)
    text_w = 36 if panel is not None else 66
    adds = _body(entry["adds"], text_w, 23)
    col = VGroup(adds_h, adds).arrange(DOWN, buff=0.3, aligned_edge=LEFT)
    if panel is not None:
        panel.scale(1.18)
        col.move_to([-3.4, -0.2, 0]).to_edge(LEFT, buff=0.9)
        panel.move_to([3.2, -0.25, 0])
    else:
        col.move_to([0, 0.0, 0])
    res = result_row(entry)
    if res is not None:
        res.move_to([0, -3.2, 0])

    scene.play(FadeIn(rib), FadeIn(head), run_time=0.7)
    scene.play(FadeIn(col, shift=UP * 0.12), run_time=0.7)
    if panel is not None:
        scene.play(FadeIn(panel), run_time=0.6)
        scene.wait(0.6)
        animate(scene)
    timing.hold_to_read(scene, adds, settle=0.8)
    if res is not None:
        scene.play(FadeIn(res, shift=UP * 0.1), run_time=0.6)
        timing.hold_to_read(scene, res, settle=1.4)
    scene.play(*[FadeOut(m) for m in scene.mobjects], run_time=0.6)
