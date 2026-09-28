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


class _BoxAxes(Axes):
    """Axes that cross at the bottom-left corner of the plotted range rather
    than at zero, so a range that spans or excludes zero does not draw its axes
    through the middle of the data."""

    @staticmethod
    def _origin_shift(axis_range):
        return axis_range[0]


def axes(x_range, y_range, x_length=9.0, y_length=3.9, font_size=16, shift_down=0.5,
         cross_at_zero=False):
    """The film's standard dark Axes (muted, no tips). They cross at the lower-left
    corner of the ranges; ``cross_at_zero=True`` restores manim's crossing at 0."""
    ax = (Axes if cross_at_zero else _BoxAxes)(
        x_range=x_range,
        y_range=y_range,
        x_length=x_length,
        y_length=y_length,
        axis_config={"color": T.MUTED, "include_tip": False, "font_size": font_size},
    )
    if shift_down:
        ax.shift(DOWN * shift_down)
    return ax


def labeled_axes(x_range, y_range, x_label=None, y_label=None, y_rotate=False,
                 numbers=False, **kw):
    """Standard Axes plus x/y labels in the house style. Returns (ax, labels).

    ``numbers=True`` adds tick numbers and places the labels clear of them.
    ``y_rotate=True`` turns the y label to run up the axis (left of the numbers)
    instead of sitting above the axis top."""
    ax = axes(x_range, y_range, **kw)
    if numbers:
        ax.add_coordinates()
    extras = VGroup()
    if x_label:
        extras.add(Text(x_label, color=T.FG, font_size=18)
                   .next_to(ax, DOWN, buff=0.12 if numbers else 0.25))
    if y_label:
        yl = Text(y_label, color=T.FG, font_size=15)
        if y_rotate:
            yl.rotate(np.pi / 2).next_to(ax, LEFT, buff=0.12)
        else:
            yl.next_to(ax.c2p(x_range[0], y_range[1]), UP, buff=0.12)
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


def histogram(ax, edges, counts, color=None, opacity=0.65, base=None):
    """Bars on ``ax`` for ``counts`` between consecutive ``edges``. ``base``
    (same length as counts) stacks the bars on top of another histogram."""
    from manim import Polygon

    c = color or T.TEAL
    g = VGroup()
    for k, n in enumerate(counts):
        y0 = float(base[k]) if base is not None else 0.0
        if n <= 0:
            continue
        a, b = ax.c2p(edges[k], y0), ax.c2p(edges[k + 1], y0 + n)
        bar = Polygon(a, [b[0], a[1], 0], b, [a[0], b[1], 0], stroke_width=0.6, color=c)
        g.add(bar.set_fill(c, opacity=opacity))
    return g


def curve(ax, xs, ys, color=None, stroke_width=3.0):
    """A polyline through data points on ``ax``."""
    from manim import VMobject

    m = VMobject(color=color or T.TEAL, stroke_width=stroke_width)
    m.set_points_as_corners([ax.c2p(x, y) for x, y in zip(xs, ys)])
    return m


def marker_line(ax, x, y_range, text, color=None, font_size=14, side=None):
    """A dashed vertical marker at ``x`` with a label at its top."""
    from manim import RIGHT, DashedLine
    from . import layout

    c = color or T.PURPLE
    line = DashedLine(ax.c2p(x, y_range[0]), ax.c2p(x, y_range[1]), color=c, stroke_width=2)
    lab = layout.label(text, font_size=font_size, color=c)
    lab.next_to(line.get_end(), side if side is not None else RIGHT, buff=0.1)
    return VGroup(line, lab)
