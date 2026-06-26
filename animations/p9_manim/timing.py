"""Reading-time estimation and reading-paced holds.

The film should leave text on screen long enough to actually read. ``reading_time``
estimates that from a word count (with a floor and a small penalty for maths/long
words); ``hold_to_read`` waits that long, plus a static "settle" pause so the
viewer gets a beat of no motion.
"""
import textwrap

# Comfortable on-screen reading speed (slower than prose: viewers also parse
# diagrams and equations). Words per minute.
WPM = 170.0
MIN_SECONDS = 1.6
# Extra seconds of stillness after the reading time elapses (the "pause").
SETTLE = 2.0


def _text_of(obj):
    """Best-effort plain text from a string or manim mobject (recursively)."""
    if obj is None:
        return ""
    if isinstance(obj, str):
        return obj
    t = getattr(obj, "text", None)
    if isinstance(t, str):
        return t
    ts = getattr(obj, "tex_string", None)
    if isinstance(ts, str):
        # strip the most common TeX noise so the word count is sane
        for ch in ("\\", "{", "}", "$", "_", "^", "~"):
            ts = ts.replace(ch, " ")
        return ts
    subs = getattr(obj, "submobjects", None)
    if subs:
        return " ".join(_text_of(s) for s in subs)
    return ""


def word_count(*objs):
    return sum(len(_text_of(o).split()) for o in objs)


def reading_time(*objs, wpm=WPM, min_seconds=MIN_SECONDS):
    """Seconds to comfortably read the given text/mobjects."""
    words = word_count(*objs)
    return max(min_seconds, words / wpm * 60.0)


def hold_to_read(scene, *objs, settle=SETTLE, extra=0.0):
    """Wait long enough to read ``objs``, then hold still for ``settle`` seconds."""
    scene.wait(reading_time(*objs) + extra)
    if settle:
        scene.wait(settle)


def wrap(text, width=52):
    """Wrap a long string to a fixed character width for on-screen Text."""
    return "\n".join(textwrap.wrap(text, width=width)) if text else ""
