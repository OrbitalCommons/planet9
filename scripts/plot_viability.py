#!/usr/bin/env python3
"""Render the Planet Nine "still-viable" parameter-space figures.

Reads figures/viability.json (produced by `cargo run -p p9-viability`) and
writes two transparent dark-theme SVGs:

  figures/p9_viable_region.svg     -- mass x distance detection map
  figures/p9_magnitude_distance.svg -- apparent magnitude vs distance

No physics here; this only draws what the Rust crate computed.
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
DATA = json.load(open(os.path.join(HERE, "figures", "viability.json")))
M = np.array(DATA["mass_earth"])


def arr(x):
    return np.array([np.nan if v is None else v for v in x], dtype=float)


# Colour by survey kind.
KIND_COLOR = {"optical": BLUE, "thermal": ORANGE, "mm": PURPLE}


def fig_viable():
    fig, ax = plt.subplots(figsize=(11, 7.2))
    wide = arr(DATA["wide_envelope_au"])
    deepest = arr(DATA["detection_envelope_au"])
    ytop = 3000.0

    # --- shaded regions ---
    # Excluded regardless of sky position (below the all-sky floor).
    ax.fill_between(M, 60, wide, color=RED, alpha=0.22, zorder=0)
    # Excluded only where the deep, small-footprint surveys looked.
    ax.fill_between(M, wide, deepest, color=ORANGE, alpha=0.16, zorder=0,
                    hatch="////", edgecolor=ORANGE, linewidth=0.0)
    # Still viable everywhere (above the deepest reach of any real survey).
    ax.fill_between(M, deepest, ytop, color=TEAL, alpha=0.13, zorder=0)

    # --- favoured current-distance band ---
    band = DATA["favored_distance_band_au"]
    ax.axhspan(band[0], band[1], color=GREEN, alpha=0.08, zorder=0)
    ax.text(19.6, np.sqrt(band[0] * band[1]), "favoured\ncurrent\ndistance",
            color=GREEN, ha="right", va="center", fontsize=8.5, alpha=0.9)

    # --- individual survey reach curves ---
    for s in DATA["surveys"]:
        y = arr(s["reach"])
        c = KIND_COLOR[s["kind"]]
        if s["status"] == "forecast":
            style = dict(ls=(0, (5, 3)), lw=1.6, alpha=0.75)
        elif s["coverage"] == "wide":
            style = dict(ls="-", lw=2.2, alpha=0.95)
        else:  # partial-footprint data
            style = dict(ls=(0, (3, 1, 1, 1)), lw=1.6, alpha=0.7)
        ax.plot(M, y, color=c, **style)
        # label at the right edge
        if np.isfinite(y[-1]):
            ax.text(M[-1] + 0.15, y[-1], s["name"], color=c, fontsize=7.8,
                    va="center", ha="left", alpha=0.95)

    # --- bold the two envelopes ---
    ax.plot(M, wide, color=RED, lw=2.8, zorder=5)
    ax.plot(M, deepest, color=ORANGE, lw=1.8, ls=":", zorder=5, alpha=0.9)

    # --- nominal orbit markers ---
    for o in DATA["nominal_orbits"]:
        ax.errorbar(o["mass_earth"], o["current_distance_au"],
                    yerr=[[o["current_distance_au"] - o["perihelion_au"]],
                          [o["aphelion_au"] - o["current_distance_au"]]],
                    fmt="o", color=GREEN, ecolor=GREEN, elinewidth=1.2,
                    capsize=3, ms=7, zorder=6, mec=FG, mew=0.8)
        ax.annotate(o["name"].split(" (")[0], (o["mass_earth"], o["current_distance_au"]),
                    textcoords="offset points", xytext=(8, 9), color=GREEN,
                    fontsize=8.0, zorder=7)

    # --- region labels ---
    ax.text(3.0, 150, "EXCLUDED over ~3/4 of sky\n(wide surveys; see search hull)", color=RED, fontsize=10,
            fontweight="bold", ha="left", va="center", alpha=0.95)
    ax.text(15.5, 820, "ruled out only inside\ndeep survey footprints",
            color=ORANGE, fontsize=8.5, ha="center", va="center", alpha=0.95)
    ax.text(6.0, 2100, "STILL VIABLE\n(too faint / too far\nto have been seen)",
            color=TEAL, fontsize=11, fontweight="bold", ha="center", va="center")

    ax.set_yscale("log")
    ax.set_xlim(1, 22.5)
    ax.set_ylim(80, ytop)
    ax.set_xlabel("Planet Nine mass  (Earth masses)")
    ax.set_ylabel("heliocentric distance  (AU)")
    ax.set_title("Planet Nine parameter space not yet eliminated by surveys",
                 fontsize=13, pad=12)
    ax.set_yticks([100, 200, 300, 500, 700, 1000, 2000, 3000])
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.grid(True, which="both", color=MUTED, alpha=0.18, lw=0.6)

    # legend for line styles
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color=RED, lw=2.8, label="all-sky exclusion floor"),
        Line2D([0], [0], color=ORANGE, lw=1.8, ls=":", label="deepest reach (any survey)"),
        Line2D([0], [0], color=FG, lw=2.2, ls="-", label="near-all-sky survey (data)"),
        Line2D([0], [0], color=FG, lw=1.6, ls=(0, (3, 1, 1, 1)), label="deep / partial footprint (data)"),
        Line2D([0], [0], color=FG, lw=1.6, ls=(0, (5, 3)), label="forecast (LSST, SO)"),
        Line2D([0], [0], marker="o", color=GREEN, lw=0, mec=FG, label="nominal orbit (q–Q range)"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=8.2, ncol=1,
              labelcolor=FG, handlelength=2.4)

    fig.text(0.013, 0.012,
             "Reaches computed from p9-core reflected-light & thermal models + the WISE/Simons reproductions. "
             "Optical=blue, thermal-IR=orange, mm=purple.",
             color=MUTED, fontsize=7.3)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = os.path.join(HERE, "figures", "p9_viable_region.svg")
    fig.savefig(out, transparent=True)
    fig.savefig(out.replace(".svg", ".png"), dpi=150, transparent=True)
    print("wrote", out)


def fig_magnitude():
    fig, ax = plt.subplots(figsize=(10, 6.4))
    dist = np.array(DATA["distance_au"])
    cmap = [BLUE, TEAL, GREEN, ORANGE, RED]
    for k, mc in enumerate(DATA["mag_curves"]):
        ax.plot(dist, mc["v_mag"], color=cmap[k % len(cmap)], lw=2.2,
                label=f"{mc['mass_earth']:.0f} M⊕")

    for name, depth in DATA["survey_optical_depths"]:
        ax.axhline(depth, color=MUTED, ls="--", lw=1.0, alpha=0.7)
        ax.text(dist[-1], depth, f" {name} ({depth:.1f})", color=MUTED,
                fontsize=8, va="center", ha="left")

    band = DATA["favored_distance_band_au"]
    ax.axvspan(band[0], band[1], color=GREEN, alpha=0.08)
    ax.text(np.sqrt(band[0] * band[1]), 14.5, "favoured\ncurrent distance",
            color=GREEN, ha="center", fontsize=8.5)

    ax.set_xlim(dist[0], dist[-1] * 1.18)
    ax.invert_yaxis()  # brighter (smaller mag) at top
    ax.set_xlabel("heliocentric distance  (AU)")
    ax.set_ylabel("apparent magnitude  (brighter ↑)")
    ax.set_title("How bright is Planet Nine? Reflected-light magnitude vs distance",
                 fontsize=13, pad=10)
    ax.grid(True, color=MUTED, alpha=0.18, lw=0.6)
    ax.legend(title="mass", loc="lower left", fontsize=9, labelcolor=FG,
              title_fontproperties={"size": 9})
    fig.text(0.013, 0.012,
             "Reflected V magnitude from p9-core planet_apparent_magnitude (Neptune albedo 0.41). "
             "A survey detects P9 where its curve lies above that survey's dashed depth line.",
             color=MUTED, fontsize=7.3)
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    out = os.path.join(HERE, "figures", "p9_magnitude_distance.svg")
    fig.savefig(out, transparent=True)
    fig.savefig(out.replace(".svg", ".png"), dpi=150, transparent=True)
    print("wrote", out)


if __name__ == "__main__":
    fig_viable()
    fig_magnitude()
