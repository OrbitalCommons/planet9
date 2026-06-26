"""Auto-generated discussion & key-learning scenes.

For every paper in ``content.PAPER_TEXT`` this defines a ``Discuss_<crate>`` scene
(hypothesis + arguments, then conclusions). For every prereq learning in
``content.PREFACE_LEARNING`` it defines a ``Learn_<scene>`` scene (the important
takeaway, held for reading plus a pause).

The assembler interleaves each of these right after its parent scene, so no
existing scene file has to be edited.
"""
from manim import Scene

from p9_manim import cards, content, paper


def _make_discuss(crate, data):
    cite = paper._short_cite(paper.crate_description(crate))

    class _D(Scene):
        def construct(self):
            cards.present_discussion(self, cite, data["hypothesis"],
                                     data["arguments"], data["conclusions"])

    name = "Discuss_" + crate.replace("-", "_")
    _D.__name__ = name
    _D.__qualname__ = name
    return name, _D


def _make_learn(key, text):
    class _L(Scene):
        def construct(self):
            cards.present_learning(self, text)

    name = "Learn_" + key
    _L.__name__ = name
    _L.__qualname__ = name
    return name, _L


for _crate, _data in content.PAPER_TEXT.items():
    _n, _c = _make_discuss(_crate, _data)
    globals()[_n] = _c

for _key, _text in content.PREFACE_LEARNING.items():
    _n, _c = _make_learn(_key, _text)
    globals()[_n] = _c
