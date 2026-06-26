"""Load data computed by the Rust crates / plotting scripts.

The film is data-driven: physics numbers come from the workspace, not from
hand-tuned constants in Python. Right now the survey-reach map is read from
``figures/viability.json`` (produced by ``cargo run -p p9-viability``).
"""
import json
import os

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(os.path.dirname(_HERE))


def repo_path(*parts):
    return os.path.join(_REPO, *parts)


def load_json(rel_path):
    """Load a JSON file relative to the repo root; return None if absent."""
    p = repo_path(*rel_path.split("/"))
    if not os.path.exists(p):
        return None
    with open(p) as fh:
        return json.load(fh)


def viability():
    """The mass x distance survey-reach dataset, or None if not generated."""
    return load_json("figures/viability.json")


_ANIM_CACHE = {}


def anim():
    """Real computed values from the reproduction crates (p9-anim-data).

    Regenerate with ``cargo run -p p9-anim-data``. Returns {} if absent so
    scenes can fall back gracefully.
    """
    if "d" not in _ANIM_CACHE:
        _ANIM_CACHE["d"] = load_json("animations/data/anim.json") or {}
    return _ANIM_CACHE["d"]


def section(name):
    """One named section of the anim dataset (or None)."""
    return anim().get(name)


def solar_system():
    """Real solar-system reference scales (bodies with a/e/q/Q/period/speed,
    plus reference magnitudes, temperatures, and survey depths)."""
    return section("solar_system")


def body(name):
    """One solar-system body's real parameters by name (or None)."""
    s = section("solar_system")
    if not s:
        return None
    for b in s["bodies"]:
        if b["name"].startswith(name):
            return b
    return None


def real_etnos(which="kbo"):
    """Real observed-object list as dicts plus the clustering R-bar.

    ``which`` = 'kbo' (6 stable KBOs, Batygin & Brown 2016) or 'etno'
    (10-object Brown 2017 sample). Returns (objects, r_bar) where each object is
    {name, a, e, i_deg, varpi_rad}; falls back to ([], None) if absent.
    """
    import math

    key = "kbo_sample" if which == "kbo" else "etno_sample"
    s = section(key)
    if not s:
        return [], None
    objs = [
        {
            "name": o["name"],
            "a": o["a"],
            "e": o["e"],
            "i_deg": o["i_deg"],
            "varpi_rad": math.radians(o["varpi_deg"]),
        }
        for o in s["objects"]
    ]
    return objs, s.get("r_bar_varpi")
