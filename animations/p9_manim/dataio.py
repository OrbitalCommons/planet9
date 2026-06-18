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
