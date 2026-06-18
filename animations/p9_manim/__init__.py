"""Shared library for the Planet Nine manim film.

Importing the package pulls in the theme (which sets the dark background) and
exposes the layout / orbit / data helpers.
"""
from . import theme, layout, orbits, dataio  # noqa: F401
from .theme import (  # noqa: F401
    BG, FG, MUTED, BLUE, GREEN, RED, ORANGE, PURPLE, TEAL, YELLOW, SUN,
)
