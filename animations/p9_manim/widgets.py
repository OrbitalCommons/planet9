"""Shared scene widgets harvested from repeated scene code.

These factor out the boilerplate that recurred across many scenes: the dark
Axes + labels, scatter "star fields", survey footprint rectangles, and small
data callouts (value rows / ladders). Keeps the per-scene files focused on
their physics rather than re-deriving layout each time.
"""
import numpy as np
from manim import (
    Axes,
    DOWN,
    Dot,
    LEFT,
    Rectangle,
    Text,
    UP,
    VGroup,
)

from . import theme as T


def axes(x_range, y_range, x_length=9.0, y_length=3.9, font_size=16, shift_down=0.5):
    """The film's standard dark Axes (muted, no tips)."""
    ax = Axes(
        x_range=x_range,
        y_range=y_range,
        x_length=x_length,
        y_length=y_length,
        axis_config={"color": T.MUTED, "include_tip": False, "font_size": font_size},
    )
    if shift_down:
        ax.shift(DOWN * shift_down)
    return ax


def labeled_axes(x_range, y_range, x_label=None, y_label=None, y_rotate=False, **kw):
    """Standard Axes plus x/y labels in the house style. Returns (ax, labels)."""
    ax = axes(x_range, y_range, **kw)
    extras = VGroup()
    if x_label:
        extras.add(Text(x_label, color=T.FG, font_size=18).next_to(ax, DOWN, buff=0.25))
    if y_label:
        yl = Text(y_label, color=T.FG, font_size=15)
        if y_rotate:
            yl.rotate(np.pi / 2).next_to(ax, LEFT, buff=0.1)
        else:
            yl.next_to(ax, UP, buff=0.05)
        extras.add(yl)
    return ax, extras


def star_field(n=90, seed=0, x=(-6.0, 6.0), y=(-2.4, 2.2), color=None,
               radius=0.02, opacity=0.5):
    """A faint scatter of background "stars"."""
    rng = np.random.default_rng(seed)
    c = color or T.MUTED
    return VGroup(*[
        Dot([rng.uniform(*x), rng.uniform(*y), 0.0], radius=radius, color=c).set_opacity(opacity)
        for _ in range(n)
    ])


def footprint_rect(width, height, color=None, fill_opacity=0.08, stroke_width=2, center=None):
    """A survey-footprint rectangle (translucent fill + border)."""
    c = color or T.PURPLE
    r = Rectangle(width=width, height=height, color=c, stroke_width=stroke_width)
    r.set_fill(c, opacity=fill_opacity)
    if center is not None:
        r.move_to(center)
    return r


def value_rows(rows, color=None, font_size=18, line_buff=0.18, align=LEFT):
    """A stacked VGroup of one-line value strings (e.g. a scale table or ladder)."""
    color = color or T.FG
    g = VGroup(*[Text(s, color=color, font_size=font_size) for s in rows])
    g.arrange(DOWN, buff=line_buff, aligned_edge=align)
    return g


def callout(title, rows, title_color=None, body_color=None, font_size=18):
    """A titled mini-table: a heading over `value_rows`."""
    head = Text(title, color=title_color or T.TEAL, font_size=font_size + 4, weight="BOLD")
    body = value_rows(rows, color=body_color, font_size=font_size)
    body.next_to(head, DOWN, buff=0.25)
    return VGroup(head, body)
