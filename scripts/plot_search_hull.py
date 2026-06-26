#!/usr/bin/env python3
"""Render the Planet Nine search-coverage and reach-hull figures.

Reads figures/search_hull.json (produced by the Rust crate) and writes two
transparent dark-theme SVGs:

  figures/p9_search_sky.svg   -- searched-depth + un-searched-probability sky map
  figures/p9_reach_hull.svg   -- distance x apparent magnitude reach hull

No physics here; this only draws what the JSON contains.
"""
import json
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

# ---- Tokyo-Night dark palette (matches the portal background) ----
BG = "none"  # transparent
FG = "#c0caf5"
MUTED = "#565f89"
BLUE = "#7aa2f7"
GREEN = "#9ece6a"
RED = "#f7768e"
ORANGE = "#e0af68"
PURPLE = "#bb9af7"
TEAL = "#7dcfff"

mpl.rcParams.update({
    "text.color": FG, "axes.labelcolor": FG, "axes.edgecolor": MUTED,
    "xtick.color": FG, "ytick.color": FG, "axes.titlecolor": FG,
    "font.size": 11, "axes.linewidth": 1.0, "figure.facecolor": BG,
    "axes.facecolor": BG, "savefig.facecolor": BG, "legend.framealpha": 0.0,
})

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = json.load(open(os.path.join(HERE, "figures", "search_hull.json")))


def arr(x):
    return np.array([np.nan if v is None else v for v in x], dtype=float)


GRID = DATA["grid"]
N_RA = GRID["n_ra"]
N_DEC = GRID["n_dec"]
RA_C = np.array(GRID["ra_centers_deg"])
DEC_C = np.array(GRID["dec_centers_deg"])


def cell(key):
    """Reshape a flat row-major cell array to (n_dec, n_ra)."""
    return arr(DATA["cells"][key]).reshape(N_DEC, N_RA)


def fig_sky():
    summ = DATA["summary"]
    excluded = 100.0 * summ["excluded_prob"]
    residual = 100.0 * summ["residual_prob"]
    space_ok = 100.0 * summ["residual_space_reachable_prob"]
    space_depth = DATA["hull"]["space_depth"]

    fig, (axt, axb) = plt.subplots(
        2, 1, figsize=(12, 9), sharex=True, constrained_layout=True)

    gal_b = cell("gal_b_deg")

    # --- TOP: searched depth map, coloured discretely by which survey is
    # deepest at each direction (so the narrow 19.5-23.8 depth range reads
    # cleanly instead of washing out). Uncovered cells stay transparent =
    # never searched. ---
    from matplotlib.colors import BoundaryNorm, ListedColormap

    surv = sorted(DATA["surveys"], key=lambda s: s["depth"])
    depths = [s["depth"] for s in surv]
    names = [s["name"] for s in surv]
    palette = [PURPLE, BLUE, TEAL, GREEN, ORANGE][: len(depths)]
    cmap = ListedColormap(palette)
    cmap.set_bad(alpha=0.0)
    edges = (
        [depths[0] - 0.5]
        + [(depths[i] + depths[i + 1]) / 2 for i in range(len(depths) - 1)]
        + [depths[-1] + 0.7]
    )
    norm = BoundaryNorm(edges, cmap.N)

    depth = np.ma.masked_invalid(cell("best_depth"))
    pcm = axt.pcolormesh(RA_C, DEC_C, depth, shading="auto", cmap=cmap,
                         norm=norm, zorder=1)
    cb = fig.colorbar(pcm, ax=axt, pad=0.012, fraction=0.045, ticks=depths)
    cb.ax.set_yticklabels(
        ["{}  ({:.1f})".format(n, d) for n, d in zip(names, depths)])
    cb.set_label("deepest survey here", color=FG)
    cb.ax.yaxis.set_tick_params(color=FG)
    plt.setp(cb.ax.get_yticklabels(), color=FG)
    axt.contour(RA_C, DEC_C, gal_b, levels=[0], colors=[MUTED],
                linestyles="--", linewidths=1.2, zorder=4)
    axt.set_title("What we've searched (coverage hull): deepest optical "
                  "survey per cell  —  dark = never searched", fontsize=12,
                  pad=8)
    axt.set_ylabel("Dec  (deg)")

    # --- BOTTOM: un-searched probability map ---
    residual_grid = cell("residual_prob")
    residual_grid = np.ma.masked_where(
        ~np.isfinite(residual_grid) | (residual_grid <= 0.0), residual_grid)
    pcm2 = axb.pcolormesh(RA_C, DEC_C, residual_grid, shading="auto",
                          cmap="inferno", zorder=1)
    cb2 = fig.colorbar(pcm2, ax=axb, pad=0.012, fraction=0.045)
    cb2.set_label("un-searched posterior probability", color=FG)
    cb2.ax.yaxis.set_tick_params(color=FG)
    plt.setp(cb2.ax.get_yticklabels(), color=FG)

    # posterior contours so the viewer sees where P9 likely is
    prob = cell("prob")
    pmax = np.nanmax(prob)
    if np.isfinite(pmax) and pmax > 0:
        levels = pmax * np.array([0.2, 0.45, 0.7, 0.9])
        axb.contour(RA_C, DEC_C, prob, levels=levels, colors=[TEAL],
                    linewidths=1.0, alpha=0.85, zorder=3)

    # galactic plane
    axb.contour(RA_C, DEC_C, gal_b, levels=[0], colors=[MUTED],
                linestyles="--", linewidths=1.2, zorder=4)
    axb.text(8, 2, "galactic plane", color=MUTED, fontsize=8.5,
             rotation=12, zorder=5)

    # top targets as ring markers
    tgts = DATA["targets"][:15]
    tra = [t["ra_deg"] for t in tgts]
    tdec = [t["dec_deg"] for t in tgts]
    axb.scatter(tra, tdec, s=90, facecolors="none", edgecolors=ORANGE,
                linewidths=1.5, zorder=6, label="top 15 targets")
    t1 = DATA["targets"][0]
    axb.scatter([t1["ra_deg"]], [t1["dec_deg"]], s=140, facecolors="none",
                edgecolors=RED, linewidths=2.0, zorder=7)
    axb.annotate(
        "#1  RA {:.1f}°, Dec {:+.1f}°\nV≈{:.1f}".format(
            t1["ra_deg"], t1["dec_deg"], t1["v_mean"]),
        (t1["ra_deg"], t1["dec_deg"]),
        textcoords="offset points", xytext=(10, 10), color=RED,
        fontsize=9, zorder=8,
        arrowprops=dict(arrowstyle="-", color=RED, lw=0.8))

    axb.set_title("Where to look next (un-searched probability) — teal = P9 "
                  "posterior, rings = best targets", fontsize=12, pad=8)
    axb.set_ylabel("Dec  (deg)")
    axb.set_xlabel("RA  (deg)")
    axb.legend(loc="lower right", fontsize=8.5, labelcolor=FG)

    for ax in (axt, axb):
        ax.set_xlim(0, 360)
        ax.set_ylim(-60, 60)
        ax.set_xticks(np.arange(0, 361, 60))
        ax.grid(True, color=MUTED, alpha=0.12, lw=0.5)

    fig.suptitle("Planet Nine: sky-coverage hull and the un-searched sky",
                 fontsize=15, color=FG, y=0.985)
    fig.text(
        0.5, 0.012,
        "{:.0f}% of the P9 posterior already excluded by past optical "
        "surveys;  {:.0f}% remains,  {:.0f}% of it reachable at "
        "V≤{:.1f} from space.".format(excluded, residual, space_ok,
                                      space_depth),
        color=MUTED, fontsize=9, ha="center")

    out = os.path.join(HERE, "figures", "p9_search_sky.svg")
    fig.savefig(out, transparent=True)
    print("wrote", out)
    return out


def fig_hull():
    h = DATA["hull"]
    dist = np.array(h["distance_au"])
    vcurve = arr(h["v_curve"])
    cloud = h["cloud"]
    cd = np.array([c["dist_au"] for c in cloud])
    cv = np.array([c["v_mag"] for c in cloud])

    allsky = h["allsky_depths"]
    deepest_allsky = max(s["depth"] for s in allsky)
    space_depth = h["space_depth"]
    space_reach = h["space_reach_au"]
    ps1 = next((s for s in allsky if s["name"] == "PS1 3pi"), allsky[-1])

    fig, ax = plt.subplots(figsize=(11, 7.2))

    # posterior cloud
    ax.scatter(cd, cv, s=6, color=BLUE, alpha=0.25, edgecolors="none",
               zorder=2, label="posterior samples")

    # fiducial V(d) curve
    ax.plot(dist, vcurve, color=FG, lw=2.4, zorder=4,
            label="V(d), {:.1f} M⊕".format(h["fiducial_mass_earth"]))

    xmin, xmax = dist.min(), dist.max()

    # shade the already-searched (all-sky) region: brighter than deepest all-sky
    ax.axhspan(0, deepest_allsky, xmin=0, xmax=1, color=ORANGE, alpha=0.10,
               zorder=0)
    ax.text(xmax, deepest_allsky - 0.15, "already searched (all-sky)",
            color=ORANGE, fontsize=8.5, ha="right", va="bottom", alpha=0.95)

    # all-sky survey depth lines + reach guides
    for s in allsky:
        ax.axhline(s["depth"], color=ORANGE, ls="--", lw=1.1, alpha=0.8,
                   zorder=3)
        ax.text(xmin, s["depth"], " {} ({:.1f})".format(s["name"], s["depth"]),
                color=ORANGE, fontsize=8, va="bottom", ha="left")
        ax.axvline(s["reach_au"], color=MUTED, ls=":", lw=0.8, alpha=0.5,
                   zorder=1)

    # space telescope depth line + reach guide
    ax.axhline(space_depth, color=GREEN, ls="-", lw=1.6, alpha=0.9, zorder=3)
    ax.text(xmin, space_depth, " space telescope ({:.1f})".format(space_depth),
            color=GREEN, fontsize=8.5, va="bottom", ha="left")
    ax.axvline(space_reach, color=TEAL, ls=":", lw=1.2, alpha=0.7, zorder=1)
    ax.text(space_reach, ax.get_ylim()[0], " space reach", color=TEAL,
            fontsize=8, rotation=90, va="bottom", ha="right", alpha=0.8)

    ax.set_xlim(xmin, xmax)
    ax.invert_yaxis()  # brighter (smaller mag) at top
    ax.set_xlabel("heliocentric distance  (AU)")
    ax.set_ylabel("apparent magnitude  (brighter ↑)")
    ax.set_title("Planet Nine reach hull: distance × brightness",
                 fontsize=13, pad=10)
    ax.grid(True, color=MUTED, alpha=0.18, lw=0.6)
    ax.legend(loc="upper left", fontsize=9, labelcolor=FG)

    fig.text(
        0.013, 0.012,
        "The space telescope reaches ~{:.0f} AU vs ~{:.0f} AU for the "
        "deepest all-sky survey ({}). Cloud points fainter than the orange "
        "band are the un-searched part of the hull.".format(
            space_reach, ps1["reach_au"], ps1["name"]),
        color=MUTED, fontsize=7.3)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = os.path.join(HERE, "figures", "p9_reach_hull.svg")
    fig.savefig(out, transparent=True)
    print("wrote", out)
    return out


if __name__ == "__main__":
    fig_sky()
    fig_hull()
