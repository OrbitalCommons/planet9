"""Auto-generated companion scenes that frame every paper and prereq scene.

For every ledger entry this defines ``Claim_<crate>`` (played before the paper's
scene: the date on the publication ribbon, the citation and the claim under
test) and ``Adds_<crate>`` (played after: what the paper added, on the gauge it
moves, with the reproduced number against the published one). For every prereq
learning in ``content.PREFACE_LEARNING`` it defines ``Learn_<scene>``.

The assembler interleaves these around their parent scene, so scene files stay
about the physics.
"""
from manim import Scene

from p9_manim import cards, content, contribution, ledger


def _named(cls, name):
    cls.__name__ = name
    cls.__qualname__ = name
    return name, cls


def _make_claim(entry):
    class _C(Scene):
        def construct(self):
            contribution.present_claim(self, entry)

    return _named(_C, "Claim_" + entry["crate"].replace("-", "_"))


def _make_adds(entry):
    class _A(Scene):
        def construct(self):
            contribution.present_contribution(self, entry)

    return _named(_A, "Adds_" + entry["crate"].replace("-", "_"))


def _make_learn(key, text):
    class _L(Scene):
        def construct(self):
            cards.present_learning(self, text)

    return _named(_L, "Learn_" + key)


for _entry in ledger.entries():
    for _make in (_make_claim, _make_adds):
        _n, _c = _make(_entry)
        globals()[_n] = _c

for _key, _text in content.PREFACE_LEARNING.items():
    _n, _c = _make_learn(_key, _text)
    globals()[_n] = _c
