#!/usr/bin/env python3
"""Render the Planet Nine cross-paper survey.

Consumes the JSON written by `cargo run -p p9-survey --bin survey`. The loader
is a tight mirror of the Rust `schema::SurveyDataset`: every struct is a typed
dataclass and `load()` rejects a file whose schema_version, keys, or field
types don't match — so the Rust and Python sides stay in lock-step (bump
SCHEMA_VERSION on both when the contract changes).

Usage:
    python3 scripts/plot_survey.py [DATA_JSON] [OUT_SVG]
Defaults: crates/p9-survey/p9_survey_data.json -> crates/p9-survey/p9_survey_sky.svg
"""
import json
import sys
from dataclasses import dataclass, fields, is_dataclass
from typing import Any, get_args, get_origin

import numpy as np

SCHEMA_VERSION = 2  # must equal p9_survey::schema::SCHEMA_VERSION

# Tokyo-Night palette (portal dark background #1a1b26)
FG, MUTE = "#c0caf5", "#565f89"
STUDY_COLORS = {  # by substring of study name, oldest->newest
    "2016 Batygin & Brown (nominal)": "#7dcfff",
    "2016 Batygin & Brown (inclined-TNO variant)": "#7aa2f7",
    "2019 Batygin et al.": "#9ece6a",
    "2021 Brown & Batygin": "#f7768e",
    "2024 Siraj": "#e0af68",
}
TELE_COLORS = {"Rubin": "#7dcfff", "Roman": "#bb9af7", "Proposed": "#9ece6a"}


# ---- tight typed schema (mirror of crates/p9-survey/src/schema.rs) ----
@dataclass
class OrbitSolution:
    name: str; citation: str; arxiv: str
    mass_earth: float; albedo: float
    a_au: float; a_sigma_au: float
    e: float; e_sigma: float
    i_deg: float; i_sigma_deg: float
    omega_deg: float; omega_sigma_deg: float
    omega_big_deg: float; omega_big_sigma_deg: float
    note: str


@dataclass
class SkyGrid:
    ra_min_deg: float; ra_max_deg: float; n_ra: int
    dec_min_deg: float; dec_max_deg: float; n_dec: int


@dataclass
class StudyResult:
    solution: OrbitSolution
    prob: list
    peak_ra_deg: float; peak_dec_deg: float
    area68_deg2: float; area95_deg2: float
    dist_p16_au: float; dist_median_au: float; dist_p84_au: float
    v_p16: float; v_median: float; v_p84: float


@dataclass
class Footprint:
    dec_min_deg: float; dec_max_deg: float
    galactic_lat_min_deg: float; coverage_fraction: float


@dataclass
class Telescope:
    name: str; band: str; limiting_mag: float
    footprint: Footprint; space_based: bool; note: str


@dataclass
class StudyDetection:
    study_name: str; detection_prob: float
    prob_in_footprint: float; prob_bright_enough: float


@dataclass
class TelescopeResult:
    telescope: Telescope
    per_study: list[StudyDetection]


@dataclass
class Overlays:
    galactic_plane: list
    ecliptic: list


@dataclass
class Constraints:
    rubin_dec_max_deg: float
    favored_nu_lo_deg: float; favored_nu_hi_deg: float
    favored_arc: list


@dataclass
class SurveyDataset:
    schema_version: int; generated_by: str
    samples_per_study: int; rng_seed: int
    grid: SkyGrid
    studies: list[StudyResult]
    telescopes: list[TelescopeResult]
    overlays: Overlays
    constraints: Constraints


def _coerce(value: Any, typ: Any, path: str) -> Any:
    """Strictly coerce `value` into the annotated type `typ`."""
    if is_dataclass(typ):
        return _from_dict(typ, value, path)
    origin = get_origin(typ)
    if origin is list:
        if not isinstance(value, list):
            raise TypeError(f"{path}: expected list, got {type(value).__name__}")
        (inner,) = get_args(typ) or (Any,)
        if is_dataclass(inner):
            return [_from_dict(inner, v, f"{path}[{i}]") for i, v in enumerate(value)]
        return value  # list of scalars / nested lists (prob grid, polylines)
    if typ is float:
        if not isinstance(value, (int, float)):
            raise TypeError(f"{path}: expected number, got {type(value).__name__}")
        return float(value)
    if typ is int:
        if not isinstance(value, int):
            raise TypeError(f"{path}: expected int, got {type(value).__name__}")
        return value
    if typ is bool:
        if not isinstance(value, bool):
            raise TypeError(f"{path}: expected bool, got {type(value).__name__}")
        return value
    if typ is str:
        if not isinstance(value, str):
            raise TypeError(f"{path}: expected str, got {type(value).__name__}")
        return value
    return value  # list/Any fall-through


def _from_dict(cls: type, d: dict, path: str = "") -> Any:
    if not isinstance(d, dict):
        raise TypeError(f"{path or cls.__name__}: expected object, got {type(d).__name__}")
    want = {f.name for f in fields(cls)}
    got = set(d.keys())
    if want != got:
        raise KeyError(
            f"{path or cls.__name__}: key mismatch; missing {want - got}, extra {got - want}"
        )
    kw = {f.name: _coerce(d[f.name], f.type, f"{path or cls.__name__}.{f.name}") for f in fields(cls)}
    return cls(**kw)


def load(path: str) -> SurveyDataset:
    with open(path) as fh:
        raw = json.load(fh)
    if raw.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(
            f"schema_version {raw.get('schema_version')} != expected {SCHEMA_VERSION}; "
            "regenerate the JSON with the matching p9-survey binary"
        )
    return _from_dict(SurveyDataset, raw)


def study_color(name: str) -> str:
    for key, col in STUDY_COLORS.items():
        if name.startswith(key) or key in name:
            return col
    return FG


def grid_2d(study: StudyResult, grid: SkyGrid) -> np.ndarray:
    return np.asarray(study.prob, float).reshape(grid.n_dec, grid.n_ra)


def levels_68_95(h: np.ndarray) -> list:
    flat = np.sort(h.ravel())[::-1]
    c = np.cumsum(flat)
    return [flat[np.searchsorted(c, q * c[-1])] for q in (0.95, 0.68)]


def main() -> None:
    data_path = sys.argv[1] if len(sys.argv) > 1 else "crates/p9-survey/p9_survey_data.json"
    out_path = sys.argv[2] if len(sys.argv) > 2 else "crates/p9-survey/p9_survey_sky.svg"
    d = load(data_path)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from scipy.ndimage import gaussian_filter
    from matplotlib.colors import LinearSegmentedColormap
    import matplotlib.patheffects as pe

    # Neutral grey "probability cloud" so the coloured overlays read clearly on
    # the dark background (the warm magma map fought the orange/red overlays).
    p9_cloud = LinearSegmentedColormap.from_list(
        "p9cloud",
        [(0.0, (0, 0, 0, 0.0)), (0.35, (0.50, 0.52, 0.60, 0.55)),
         (1.0, (0.90, 0.92, 0.99, 1.0))],
    )
    outline = [pe.Stroke(linewidth=4.0, foreground="#1a1b26"), pe.Normal()]

    g = d.grid
    extent = [g.ra_min_deg, g.ra_max_deg, g.dec_min_deg, g.dec_max_deg]
    primary = next(s for s in d.studies if "2021 Brown" in s.solution.name)

    fig = plt.figure(figsize=(9.6, 7.4))
    fig.patch.set_alpha(0)
    gs = GridSpec(2, 1, height_ratios=[2.4, 1.0], hspace=0.42)

    # ---- panel 1: sky heatmap ----
    ax = fig.add_subplot(gs[0]); ax.set_facecolor("none")
    H = gaussian_filter(grid_2d(primary, g), 1.2)
    im = ax.imshow(H / H.max(), origin="lower", extent=extent, aspect="auto",
                   cmap=p9_cloud, vmin=0.02, zorder=2)
    xc = (np.arange(g.n_ra) + 0.5) * (g.ra_max_deg - g.ra_min_deg) / g.n_ra + g.ra_min_deg
    yc = (np.arange(g.n_dec) + 0.5) * (g.dec_max_deg - g.dec_min_deg) / g.n_dec + g.dec_min_deg
    for s in d.studies:
        col = study_color(s.solution.name)
        if s is not primary:
            Hs = gaussian_filter(grid_2d(s, g), 1.4)
            ax.contour(xc, yc, Hs, levels=levels_68_95(Hs), colors=col,
                       linewidths=[0.9, 1.6], alpha=0.95, zorder=4)
        ax.plot(s.peak_ra_deg, s.peak_dec_deg, "*", ms=14, mfc=col, mec="white",
                mew=0.5, zorder=6)
        ax.plot([], [], color=col, lw=2.2,
                label=f"{s.solution.name.split('(')[0].strip()} — V≈{s.v_median:.1f}, {s.dist_median_au:.0f} AU")

    gp = np.array(d.overlays.galactic_plane); ec = np.array(d.overlays.ecliptic)
    op = np.argsort(gp[:, 0]); ax.plot(gp[op, 0], gp[op, 1], ls="--", color="#9aa5ce", lw=1.8,
                                       alpha=0.9, zorder=3, path_effects=outline,
                                       label="Galactic plane (crowded zone)")
    oe = np.argsort(ec[:, 0]); ax.plot(ec[oe, 0], ec[oe, 1], ls=":", color=FG, lw=0.9,
                                       alpha=0.5, zorder=3, label="Ecliptic")

    # Search-narrowing overlays: Rubin cede line + Cassini favored-ν arc.
    c = d.constraints
    ax.axhspan(g.dec_min_deg, c.rubin_dec_max_deg, color=MUTE, alpha=0.16, zorder=1)
    ax.axhline(c.rubin_dec_max_deg, color="#c0caf5", ls=(0, (6, 3)), lw=1.4, alpha=0.9,
               zorder=5, label=f"Rubin reach (Dec < +{c.rubin_dec_max_deg:.0f}°)")
    # The favored-ν arc is a compact curve that can cross RA = 0/360; draw it as
    # connected segments (split on wrap) in a high-contrast colour with a dark
    # outline so it reads on the grey cloud.
    arc = np.array(c.favored_arc)
    ra_a, dec_a = arc[:, 0], arc[:, 1]
    start = 0
    for k in range(1, len(ra_a) + 1):
        if k == len(ra_a) or abs(ra_a[k] - ra_a[k - 1]) > 180:
            sl = slice(start, k)
            ax.plot(ra_a[sl], dec_a[sl], color="#bb9af7", lw=2.8, zorder=8,
                    path_effects=outline)
            start = k
    ax.plot([], [], color="#bb9af7", lw=2.8,
            label=f"Cassini favoured-ν arc (ν {c.favored_nu_lo_deg:.0f}–{c.favored_nu_hi_deg:.0f}°, B&B orbit)")

    ax.set_xlim(g.ra_min_deg, g.ra_max_deg); ax.set_ylim(g.dec_min_deg, g.dec_max_deg)
    ax.invert_xaxis()
    ax.set_xticks(np.arange(0, 361, 30))
    ax.set_xticklabels([f"{int(t/15)}h" for t in np.arange(0, 361, 30)], color=FG, fontsize=8)
    ax.set_yticks(np.arange(-45, 46, 15))
    ax.set_yticklabels([f"{t:+d}°" for t in np.arange(-45, 46, 15)], color=FG, fontsize=8)
    ax.set_xlabel("Right Ascension", color=FG); ax.set_ylabel("Declination", color=FG)
    ax.set_title("Where is Planet 9? Dwell-time position probability across orbit solutions",
                 color=FG, fontsize=12, pad=8)
    for sp in ax.spines.values(): sp.set_color(MUTE)
    ax.tick_params(colors=MUTE)
    ax.legend(loc="lower left", fontsize=7.2, framealpha=0.0, labelcolor=FG, ncol=2)
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01)
    cb.set_label("relative position probability (2021 MCMC)", color=FG, fontsize=8)
    cb.ax.tick_params(colors=MUTE); plt.setp(cb.ax.get_yticklabels(), color=FG, fontsize=7)

    # ---- panel 2: detection probability per telescope across studies ----
    ax2 = fig.add_subplot(gs[1]); ax2.set_facecolor("none")
    studies = [s.solution.name.split("(")[0].strip() for s in d.studies]
    x = np.arange(len(studies)); w = 0.26
    for j, t in enumerate(d.telescopes):
        col = next((c for k, c in TELE_COLORS.items() if k in t.telescope.name), FG)
        probs = [next(p.detection_prob for p in t.per_study
                      if p.study_name == s.solution.name) for s in d.studies]
        ax2.bar(x + (j - 1) * w, np.array(probs) * 100, w, color=col, alpha=0.9,
                label=f"{t.telescope.name.split('(')[0].strip()} (m<{t.telescope.limiting_mag:.1f})")
    ax2.set_xticks(x)
    ax2.set_xticklabels([s.replace("Batygin & Brown", "B&B").replace("Batygin et al.", "Batygin+")
                         .replace("Siraj, Chyba & Tremaine", "Siraj+") for s in studies],
                        color=FG, fontsize=7.5, rotation=12)
    ax2.set_ylabel("P(detect) [%]", color=FG, fontsize=9)
    ax2.set_ylim(0, 100)
    ax2.set_title("Can a survey catch it? Detection probability = P(in footprint AND brighter than depth)",
                  color=FG, fontsize=10.5, pad=6)
    for sp in ax2.spines.values(): sp.set_color(MUTE)
    ax2.tick_params(colors=MUTE); plt.setp(ax2.get_yticklabels(), color=FG, fontsize=7)
    ax2.legend(loc="upper right", fontsize=7.2, framealpha=0.0, labelcolor=FG)
    ax2.grid(axis="y", color=MUTE, alpha=0.25, lw=0.5)

    fig.savefig(out_path, transparent=True, bbox_inches="tight")
    print(f"wrote {out_path}  (schema v{d.schema_version}, "
          f"{d.samples_per_study} samples/study, {len(d.studies)} studies)")


if __name__ == "__main__":
    main()
