# 03 — Layer 2: LSST coverage-to-date → search hull (Phase 2)

Purpose: keep goal G2 — the hull/viability numbers quoted anywhere in this
repo reflect what Rubin has *actually* observed, updated monthly. Today the
hull treats Rubin as a forecast; within a year LSST will have converted a
large share of the southern residual, and our maps must track that
conversion at the cadence it happens.

## Source options, ranked

| Rank | Source | Access | Notes |
|---|---|---|---|
| 1 | Alert-stream-derived coverage | world-public via Fink Data Transfer per-night pulls (Layer 3) | every processed visit emits alerts; (night, band, pointing) is reconstructible from alert metadata. Zero new dependencies once Layer 3 exists; blind to visits whose alert production failed (rare, acceptable) |
| 2 | Published completed-visit metadata (scheduler/`schedview`-style nightly summaries) | public web, format to be confirmed at implementation | cleaner than (1) if a stable machine-readable feed exists; treat as an upgrade, verify during Phase 2 |
| 3 | RSP ConsDB/visit tables | data-rights gated | deliberately not a dependency (non-goal); revisit only if (1) and (2) both prove inadequate |

Decision D-5: build against (1), spike (2) for one afternoon during
implementation, ignore (3).

## Map specification (SR-7)

`rubin_watch/lsst_coverage.json` (schema in design/07):

- HEALPix, **nside 64** (49,152 pixels, ~0.84 deg² each), ring ordering.
- Per band in {u, g, r, i, z, y}, but the hull consumes r (and g as a
  cross-check): `n_visits[pix]`, `last_visit_mjd[pix]`,
  `depth_single_median[pix]` (5σ PSF from alert metadata where available,
  else the band's fiducial: r ≈ 24.3 single visit).
- Derived, computed in Rust at hull-integration time (never stored, so the
  model can improve without re-ingesting):
  - `depth_catalog(pix)` — single-visit depth (a bright mover is caught on
    any one visit; this is the conservative reachability depth),
  - `depth_linked(pix, N)` — the magnitude at which the probability of ≥ N
    detections over the pixel's visit history reaches 95%, via the existing
    `p9-2023-lsst-strategy` P(≥N of M) machinery (`poisson_binomial_tail` +
    per-visit logistic efficiency). N = 3 matches SSP-like linkability and
    our own linker's requirement (SR-5).

## Hull and viability integration

- `p9-search-hull` gains an optional input: if `rubin_watch/
  lsst_coverage.json` exists, add a `RealSurvey`-like entry **"LSST (to
  date)"** whose acceptance is per-pixel membership (n_visits ≥ 3 in r)
  rather than an analytic footprint. Implementation: extend the hull's
  survey layer with a pixel-mask variant alongside `Footprint` — the
  `max_distance_au` cap machinery (from the TESS fix) already generalizes
  the struct; this adds a `PixelMask(HealpixMask)` footprint variant.
  Depth per cell = `depth_linked(pix, 3)` for exclusion (an object must be
  *linkable* to be excluded by non-detection), and the map keeps the
  existing Rubin *forecast* entry separately so "now" vs "10-year" remain
  distinct products.
- `p9-viability` swaps the static "Rubin/LSST (visit)" row's sky fraction
  for the observed fraction (pixels with ≥ 3 r visits / 4π).
- Regeneration cadence: monthly, or immediately after any Layer-1 delta
  that changes the posterior sample (design/06 decision rules).

## Depth model honesty

Two labeled systematics, both documented in the artifact header:

1. Alert-derived depth is the *difference-image* point-source depth; for a
   moving point source on a clean template it is the right number, but
   year-2+ template self-subtraction (design/05 §hazards) degrades the
   effective depth for ultra-slow movers. The coverage artifact therefore
   records `template_epoch_flag` per pixel once templates begin including
   survey-year data — consumers decide the penalty.
2. Near-plane pixels: DIA crowding losses are not modeled; pixels with
   |b| < 10° carry `crowding_flag = true` and the hull continues to treat
   them conservatively (the same deliberate-fidelity choice as the ZTF
   plane treatment, documented at the survey entry).

## Acceptance tests

1. **Synthetic full coverage**: a fixture map with every southern pixel at
   200 visits reproduces the hull's existing Rubin-forecast exclusion
   numbers within 1 point (consistency between "to date" machinery and the
   forecast entry).
2. **Monotonicity**: for fixture months M1 ⊂ M2 (visit sets), excluded% is
   non-decreasing and residual% non-increasing.
3. **Empty map**: absent/empty coverage file → hull output byte-identical
   to today's (the entry is strictly additive).
4. **Depth-linked sanity**: a pixel with exactly 3 visits at r=24.3 yields
   depth_linked(3) < depth_single (needing 3 hits is harder than 1), and
   depth_linked grows with visit count toward the logistic ceiling.
