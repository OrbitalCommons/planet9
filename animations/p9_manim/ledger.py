"""The publication ledger: one dated entry per reproduced paper, and the state
of the field that each entry changes.

Entries live one-per-file in ``animations/ledger/<crate>.yaml`` so a paper's
scene, data export and ledger entry can be edited independently. The ledger
orders the papers in time and replays them to answer, for any paper, "what did
the field know the day before this appeared?" -- the *before* half of every
contribution beat.

Entry schema (YAML)::

    crate: p9-2021-ztf
    cite: Brown & Batygin (2022)
    arxiv: "2110.13117"
    date: "2021-10"            # arXiv posting month
    kind: search               # see KINDS
    headline: First big bite out of the search
    claim: one sentence -- what the paper sets out to show
    adds: one or two sentences -- what the field knew after it that it did
          not know before
    track: excluded            # clustering | orbit | excluded | sample | none
    excluded:   {label: ZTF, cumulative: "@combined_after"}
    clustering: {sigma: 3.8, stance: pro, label: BB16}
    orbit:      {mass: 6.2, a: 380, e: 0.30, i: 16}
    sample:     {n: 14, what: "ETNOs with a > 250 AU"}
    result:     {label: excluded fraction, key: ztf, fmt: pct1, published: 0.564}

Any numeric field may be written ``"@key"``: it is then read from the Rust
export ``anim.json -> papers -> <crate> -> key`` so the film draws what the
reproduction crate computed.
"""
import glob
import os

import yaml

from . import dataio
from . import theme as T

KINDS = {
    "evidence": ("EVIDENCE", T.GREEN),
    "discovery": ("DISCOVERY", T.GREEN),
    "critique": ("CRITIQUE", T.RED),
    "dynamics": ("DYNAMICS", T.BLUE),
    "indirect": ("INDIRECT SIGNATURE", T.ORANGE),
    "search": ("SEARCH", T.PURPLE),
    "forecast": ("FORECAST", T.TEAL),
    "alternative": ("ALTERNATIVE", "#ff9e64"),
    "physical": ("PHYSICAL MODEL", "#e0af68"),
    "review": ("SYNTHESIS", T.FG),
}

TRACKS = ("clustering", "orbit", "excluded", "sample")

_LEDGER_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "ledger")
_CACHE = {}


def _month_index(date):
    """'2021-10' -> fractional year 2021.75 (mid-month)."""
    y, m = str(date).split("-")[:2]
    return int(y) + (int(m) - 0.5) / 12.0


def _resolve(crate, value):
    """Resolve an ``@key`` reference against the Rust export; pass numbers through.
    Returns None when the key (or the export) is missing."""
    if isinstance(value, str) and value.startswith("@"):
        papers = dataio.section("papers") or {}
        return (papers.get(crate) or {}).get(value[1:])
    return value


def _resolve_payload(crate, payload):
    if not isinstance(payload, dict):
        return payload
    return {k: _resolve(crate, v) for k, v in payload.items()}


def entries():
    """Every ledger entry, oldest first (ties broken by crate name)."""
    if "entries" not in _CACHE:
        out = []
        for path in sorted(glob.glob(os.path.join(_LEDGER_DIR, "*.yaml"))):
            with open(path) as fh:
                e = yaml.safe_load(fh)
            e["t"] = _month_index(e["date"])
            for track in TRACKS:
                if track in e:
                    e[track] = _resolve_payload(e["crate"], e[track])
            if "result" in e:
                r = e["result"]
                r["reproduced"] = _resolve(e["crate"], "@" + r["key"]) if r.get("key") else None
            out.append(e)
        out.sort(key=lambda e: (e["t"], e["crate"]))
        _CACHE["entries"] = out
    return _CACHE["entries"]


def entry(crate):
    for e in entries():
        if e["crate"] == crate:
            return e
    raise KeyError(f"no ledger entry for {crate}")


def has(crate):
    return any(e["crate"] == crate for e in entries())


def earlier(crate, track=None):
    """Entries strictly before ``crate`` in ledger order, optionally only those
    carrying a payload on ``track``."""
    out = []
    for e in entries():
        if e["crate"] == crate:
            break
        if track is None or e.get(track):
            out.append(e)
    return out


def kind_label(e):
    return KINDS.get(e.get("kind"), ("PAPER", T.FG))


def span():
    """(first, last) fractional years covered by the ledger, padded to whole years."""
    ts = [e["t"] for e in entries()]
    return int(min(ts)), int(max(ts)) + 1


def format_value(value, fmt):
    """Format a result number: pct0/pct1 (fractions -> percent), sigma, deg, au,
    mearth, mag, or a plain ``{:...}`` format spec."""
    if value is None:
        return None
    if fmt in ("pct0", "pct1", "pct2"):
        return f"{100.0 * value:.{fmt[-1]}f}%"
    if fmt == "sigma":
        return f"{value:.1f}σ"
    if fmt == "deg":
        return f"{value:.1f}°"
    if fmt == "au":
        return f"{value:.0f} AU"
    if fmt == "mearth":
        return f"{value:.1f} M⊕"
    if fmt == "mag":
        return f"V = {value:.1f}"
    if fmt and "{" in fmt:
        return fmt.format(value)
    return f"{value:.3g}"


MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]


def date_label(e):
    y, m = str(e["date"]).split("-")[:2]
    return f"{MONTHS[int(m) - 1]} {y}"
