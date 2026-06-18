"""Helpers for per-paper scenes: title cards, citation chips pulled from each
crate's Cargo.toml, and a reproduced-result readout.
"""
import os
import re

from manim import DOWN, UP, MathTex, SurroundingRectangle, Text, VGroup

from . import theme as T
from . import dataio


def crate_description(crate):
    """First line of the crate's Cargo.toml `description`, or the crate name."""
    path = dataio.repo_path("crates", crate, "Cargo.toml")
    if os.path.exists(path):
        with open(path) as fh:
            for line in fh:
                m = re.match(r"\s*description\s*=\s*\"(.*)\"", line)
                if m:
                    return m.group(1)
    return crate


def title_block(crate, headline):
    """A scene-opening title: the paper headline + an auto citation chip."""
    title = Text(headline, color=T.FG, font_size=34, weight="BOLD").to_edge(UP, buff=0.6)
    cite = Text(_short_cite(crate_description(crate)), color=T.MUTED, font_size=18)
    cite.next_to(title, DOWN, buff=0.18)
    return VGroup(title, cite)


def _short_cite(desc):
    """Trim a long Cargo description down to an author-year(-ish) chip."""
    # e.g. "Reproduction of Batygin & Brown (2016) 'Evidence for ...'"
    m = re.search(r"Reproduction of (.+?\(\d{4}[a-z]?\))", desc)
    if m:
        return m.group(1)
    return desc[:60]


def result_readout(label, value, color=T.GREEN):
    """A boxed 'reproduced: <value>' readout for the paper's headline number."""
    lab = Text(label, color=T.MUTED, font_size=20)
    val = Text(value, color=color, font_size=30, weight="BOLD")
    val.next_to(lab, DOWN, buff=0.12)
    g = VGroup(lab, val)
    box = SurroundingRectangle(g, color=color, buff=0.25, corner_radius=0.1)
    box.set_fill(color, opacity=0.06)
    return VGroup(box, g)
