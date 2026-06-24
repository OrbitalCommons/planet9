"""Reusable on-screen furniture: title cards, citation chips, captions, takeaways,
classy equations, and reading-paced end-of-section beats."""
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    UR,
    DL,
    Create,
    FadeIn,
    FadeOut,
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


def label(text, font_size=18, color=None, weight="NORMAL", **kw):
    """A small Text label with even kerning. manim's Text lays glyphs out
    unevenly below ~20px, so we render at a safe size and scale down to the
    requested visual size."""
    internal = max(font_size, 28)
    t = Text(text, font_size=internal, color=color or T.FG, weight=weight, **kw)
    if internal != font_size:
        t.scale(font_size / internal)
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


def explain_equation(scene, parts, explain, color=None, accent=None, scale=1.0,
                     card=True, where=None):
    """Show an equation, then walk through it TERM BY TERM: highlight each term
    and show a plain-language label for it, in turn, before a final reading hold.

    ``parts`` is a list of LaTeX fragments (e.g. ``[r"v^2", "=", "GM",
    r"\\left(\\frac{2}{r}-\\frac{1}{a}\\right)"]``); each becomes an indexable
    submobject (this avoids the command-splitting pitfalls of substring
    isolation). ``explain`` is a list of ``(part_index, "what it means")``.
    Returns the equation group so the caller can fade it out.
    """
    color = color or T.FG
    accent = accent or T.TEAL
    eq = MathTex(*parts).set_color(color).scale(scale)
    eq.move_to(where if where is not None else UP * 1.7)
    group = eq
    if card:
        rect = SurroundingRectangle(eq, color=T.MUTED, buff=0.3, corner_radius=0.16)
        rect.set_stroke(T.MUTED, width=1.3, opacity=0.6).set_fill("#222436", opacity=0.7)
        group = VGroup(rect, eq)
        scene.play(Create(rect), Write(eq), run_time=1.2)
    else:
        scene.play(Write(eq), run_time=1.1)

    for idx, meaning in explain:
        part = eq[idx] if 0 <= idx < len(eq) else None
        lbl = Text(meaning, color=accent, font_size=T.SMALL_SIZE)
        lbl.to_edge(DOWN, buff=0.7)
        if part is not None:
            box = SurroundingRectangle(part, color=accent, buff=0.06, corner_radius=0.05)
            scene.play(part.animate.set_color(accent), Create(box), FadeIn(lbl), run_time=0.45)
            timing.hold_to_read(scene, meaning, settle=0.4)
            scene.play(part.animate.set_color(color), FadeOut(box), FadeOut(lbl), run_time=0.4)
        else:
            scene.play(FadeIn(lbl), run_time=0.4)
            timing.hold_to_read(scene, meaning, settle=0.4)
            scene.play(FadeOut(lbl), run_time=0.35)
    timing.hold_to_read(scene, eq, settle=1.0)
    return group


# ---- reading-paced end-of-section beat --------------------------------------

def show_takeaway(scene, text, settle=2.5):
    """Fade in the takeaway and HOLD it long enough to read without pausing the
    video (reading time + a still settle pause). Use at the end of every scene."""
    box = takeaway(text)
    scene.play(FadeIn(box, shift=UP * 0.15))
    timing.hold_to_read(scene, text, settle=settle)
    return box
