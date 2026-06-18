"""Reusable on-screen furniture: title cards, citation chips, captions, takeaways."""
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    UR,
    DL,
    MathTex,
    Rectangle,
    SurroundingRectangle,
    Text,
    VGroup,
)

from . import theme as T


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
