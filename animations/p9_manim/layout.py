"""Reusable on-screen furniture: title cards, citation chips, captions, takeaways,
classy equations, and reading-paced end-of-section beats."""
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    UR,
    DL,
    FadeIn,
    MathTex,
    Rectangle,
    SurroundingRectangle,
    Text,
    VGroup,
    Write,
)

from . import theme as T
from . import timing


def title_card(title, subtitle=None):
    """Centred title (+ optional subtitle) as a VGroup."""
    parts = [Text(title, color=T.FG, font_size=T.TITLE_SIZE, weight="BOLD")]
    if subtitle:
        sub = Text(subtitle, color=T.MUTED, font_size=T.SMALL_SIZE)
        sub.next_to(parts[0], DOWN, buff=0.35)
        parts.append(sub)
    return VGroup(*parts)


def concept_badge(label):
    """Small pinned 'PREFACE 03' style chip, top-left."""
    chip = Text(label, color=T.TEAL, font_size=20, weight="BOLD")
    chip.to_corner(UP + LEFT, buff=0.4)
    return chip


def citation_chip(text):
    """Author-year / arXiv chip, pinned top-right."""
    t = Text(text, color=T.MUTED, font_size=18)
    t.to_corner(UR, buff=0.4)
    return t


def caption(text, font_size=None):
    """Bottom caption line."""
    t = Text(text, color=T.FG, font_size=font_size or T.SMALL_SIZE)
    t.to_edge(DOWN, buff=0.45)
    return t


def takeaway(text):
    """A highlighted one-line takeaway box near the bottom."""
    t = Text(text, color=T.FG, font_size=T.SMALL_SIZE, weight="BOLD")
    box = SurroundingRectangle(t, color=T.TEAL, buff=0.25, corner_radius=0.12)
    box.set_fill(T.TEAL, opacity=0.08)
    g = VGroup(box, t)
    g.to_edge(DOWN, buff=0.5)
    return g


def footer_credit(text="p9-core models"):
    t = Text(text, color=T.MUTED, font_size=14)
    t.to_corner(DL, buff=0.25)
    return t


# ---- classy equations -------------------------------------------------------

def equation(tex, color=None, scale=1.0, t2c=None):
    """A LaTeX equation in the house style: serif math, accent colour, optional
    per-symbol colouring via ``t2c`` (a {substring: color} map)."""
    m = MathTex(tex, color=color or T.TEAL, tex_to_color_map=t2c or {})
    return m.scale(scale)


def equation_card(tex, color=None, scale=1.0, t2c=None):
    """A feature equation seated on a faint rounded card — the classy 'hero'
    treatment for a scene's governing relation."""
    m = equation(tex, color=color or T.FG, scale=scale, t2c=t2c)
    card = SurroundingRectangle(m, color=T.MUTED, buff=0.32, corner_radius=0.16)
    card.set_stroke(T.MUTED, width=1.3, opacity=0.6).set_fill("#222436", opacity=0.7)
    return VGroup(card, m)


def show_equation(scene, mobj, settle=1.0, run_time=1.1):
    """Write an equation in and hold it long enough to read (reading-paced)."""
    scene.play(Write(mobj), run_time=run_time)
    timing.hold_to_read(scene, mobj, settle=settle)
    return mobj


# ---- reading-paced end-of-section beat --------------------------------------

def show_takeaway(scene, text, settle=2.5):
    """Fade in the takeaway and HOLD it long enough to read without pausing the
    video (reading time + a still settle pause). Use at the end of every scene."""
    box = takeaway(text)
    scene.play(FadeIn(box, shift=UP * 0.15))
    timing.hold_to_read(scene, text, settle=settle)
    return box
