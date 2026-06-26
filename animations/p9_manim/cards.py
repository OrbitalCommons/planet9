"""Reusable text "cards" and the discussion / key-learning presenters.

These render the slower, text-heavy beats: a paper's hypothesis + arguments, its
conclusions, and a prereq scene's key takeaway -- each paced by ``timing`` so the
viewer can read it, followed by a static pause.
"""
from manim import DOWN, FadeIn, FadeOut, LEFT, RIGHT, Text, UP, VGroup

from . import theme as T
from . import timing


def _heading(text, color):
    return Text(text, color=color, font_size=30, weight="BOLD")


def bullets(items, color=None, font_size=24, width=58, line_buff=0.28):
    """A left-aligned bulleted list; each item is wrapped to ``width`` chars."""
    color = color or T.FG
    rows = VGroup()
    for it in items:
        wrapped = timing.wrap(it, width=width)
        dot = Text("•", color=color, font_size=font_size)
        body = Text(wrapped, color=color, font_size=font_size, line_spacing=0.8)
        body.next_to(dot, RIGHT, buff=0.18, aligned_edge=UP)
        rows.add(VGroup(dot, body))
    rows.arrange(DOWN, buff=line_buff, aligned_edge=LEFT)
    return rows


def hypothesis_card(hypothesis, arguments, accent=T.BLUE):
    """A 'Hypothesis' heading + claim, then an 'Arguments' bulleted list."""
    head = _heading("Hypothesis", accent).to_edge(UP, buff=0.7)
    claim = Text(timing.wrap(hypothesis, width=56), color=T.FG, font_size=27, line_spacing=0.85)
    claim.next_to(head, DOWN, buff=0.35)
    argh = Text("Arguments presented", color=accent, font_size=22, weight="BOLD")
    argh.next_to(claim, DOWN, buff=0.5).to_edge(LEFT, buff=1.0)
    blist = bullets(arguments, color=T.FG, font_size=23).next_to(argh, DOWN, buff=0.25).to_edge(LEFT, buff=1.0)
    return VGroup(head, claim, argh, blist)


def conclusions_card(points, accent=T.TEAL):
    """A 'Conclusions' heading + bulleted takeaways."""
    head = _heading("Conclusions", accent).to_edge(UP, buff=0.7)
    blist = bullets(points, color=T.FG, font_size=25).next_to(head, DOWN, buff=0.5).to_edge(LEFT, buff=1.0)
    return VGroup(head, blist)


def key_learning_card(text, accent=T.TEAL):
    """A 'Key takeaway' heading + the one important learning, centred."""
    head = _heading("Key takeaway", accent).to_edge(UP, buff=1.1)
    body = Text(timing.wrap(text, width=46), color=T.FG, font_size=30, line_spacing=0.95,
                weight="BOLD")
    body.next_to(head, DOWN, buff=0.6)
    return VGroup(head, body)


def _clear(scene):
    if scene.mobjects:
        scene.play(*[FadeOut(m) for m in scene.mobjects], run_time=0.6)


def present_discussion(scene, citation, hypothesis, arguments, conclusions):
    """Hypothesis + arguments beat, then (after the diagram) the conclusions beat."""
    from . import layout
    _clear(scene)
    chip = layout.citation_chip(citation) if citation else None
    hyp = hypothesis_card(hypothesis, arguments)
    if chip:
        scene.add(chip)
    scene.play(FadeIn(hyp, shift=UP * 0.2))
    timing.hold_to_read(scene, hyp)
    scene.play(FadeOut(hyp))

    con = conclusions_card(conclusions)
    scene.play(FadeIn(con, shift=UP * 0.2))
    timing.hold_to_read(scene, con)
    if chip:
        scene.play(FadeOut(con), FadeOut(chip))
    else:
        scene.play(FadeOut(con))


def present_learning(scene, text):
    """A prereq scene's 'important learning' recap, held still for reading + a pause."""
    card = key_learning_card(text)
    scene.play(FadeIn(card, shift=UP * 0.2))
    timing.hold_to_read(scene, card, settle=3.0)
    scene.play(FadeOut(card))
