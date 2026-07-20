# Runbook — coverage map (Layer 2)

Monthly (or after any significant Data Transfer pull): rebuild the observed
LSST coverage map and regenerate the hull/viability products. Design:
`design/03`.

## The run

```bash
# 1. build the map from the alert store (healpy venv)
~/.venvs/fink/bin/python scripts/rubin_watch/build_coverage.py

# 2. validate (schema gate — exit 2 on malformed)
cargo run --release -p p9-rubin-watch -- coverage

# 3. regenerate the consumers
cargo run --release -p p9-search-hull
cargo run --release -p p9-viability
```

The hull prints `including LSST (to date) ... N linkable pixels` when the
map is present; without the file its output is byte-identical to the
pre-coverage hull (tested).

## What the numbers mean

- **covered** pixel: ≥ 1 reconstructed visit touched it (visit = detections
  clustered on band + exposure time, footprint painted to the 1.75° camera
  radius — `--mode detections` gives the strict lower bound instead).
- **linkable** pixel: ≥ 3 r-band visits; only these count for exclusion,
  at the LINKED depth (P(≥3 detections) = 95% via
  `p9_core::analysis::surveys::linked_depth` — shallower than single-visit
  at exactly 3 visits, deeper once visits accumulate).
- Crowding-flagged pixels (|b| < 10°) never enter the hull entry — same
  conservative plane treatment as the analytic surveys.
- Depths are the labeled band fiducials (Light-SSO packets carry no
  per-visit maglim); the builder prints faintest-detection vs fiducial as
  the over-claim check (2026-07-19: r 24.39 obs vs 24.3 fiducial — snug).

## Landing

Commit `rubin_watch/lsst_coverage.json` + regenerated `figures/*.json` on a
`meawoppl/rubin-watch-YYYY-MM-DD-coverage` branch; standard gate; note the
excluded% movement in the PR body (design/06 threshold: mention always,
regenerate-and-report when ≥ 0.5 pt).

## Caveats to keep honest

- Visit reconstruction from selection-A alerts only sees pointings that
  produced SSO detections; sparse nights undercount. Selection-B pulls will
  densify this; until then the map is a lower bound on Rubin's true
  coverage — which biases our residual-probability numbers HIGH (safe
  direction).
- The map indexes with healpy conventions; the Rust side
  (`p9_core::coords::healpix`) is fixture-locked to healpy, and
  astropy-healpix is deliberately not the reference (different boundary
  convention at the exact equator).
