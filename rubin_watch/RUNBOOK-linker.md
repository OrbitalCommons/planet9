# Runbook — slow-mover linker (Layer 4)

Monthly, after the coverage rebuild (or after any selection-B pull). Design:
`design/05`. The linker hunts d ∈ [300, 1000] AU movers — no intra-night
tracklets, 3.4–11.2″/day at opposition — via distance-hypothesis linking.

## The run

```bash
# 1. export tile inputs from the store (unassociated detections + the
#    reconstructed visit list; --include-known only for smoke tests)
~/.venvs/fink/bin/python scripts/rubin_watch/export_linker_input.py

# 2. link each tile (exit 1 = candidates written to rubin_watch/candidates/)
for f in ~/.cache/p9-rubin-watch/linker/tile-*.json; do
  cargo run --release -p p9-rubin-watch -- link --input "$f"
done

# 3. measure completeness on the same real cadence (SR-6)
cargo run --release -p p9-rubin-watch -- link \
  --input ~/.cache/p9-rubin-watch/linker/tile-XXXX.json --calibrate 40
```

Uses the DE421 ephemeris when `~/.cache/starfield/de421.bsp` exists; falls
back to the circular-Earth model with a loud `FALLBACK` tag (fine for
smoke tests, not for candidates).

## How it works (one paragraph)

For each 1/d hypothesis (~580 steps over 300–1000 AU), every detection is
inverted to a heliocentric ecliptic direction; a real body at that distance
collapses to a point drifting at its own mean motion (≤ 0.7″/day) while
stars and wrong-distance movers smear. Drift-aware λ-sweep clustering
proposes member sets with ≥ 3 distinct Chile nights; each set is refit with
the full 5-parameter model over the ENTIRE distance range (the cluster's
hypothesis is a proposal, not a prior — residual-parallax chains can link a
true object at wrong hypotheses), and must pass χ²/dof ≤ 3, beat a
static-source model by Δχ² ≥ 25, and land in a loose implied-H window
before becoming a candidate record.

## Gates and knobs (linker.rs constants)

| Constant | Value | Why |
|---|---|---|
| R range | 300–1000 AU | posterior residual band; below 300 the SSP/community own it |
| INV_D_STEP | 4e-6 AU⁻¹ | ≤ 2″ smear over 90-day windows |
| LINK_TOL / DRIFT_MAX | 3″ / 1.5″ day⁻¹ | transform residual + noise / own motion ceiling with margin |
| MIN_NIGHTS | 3 | over-determines the 5-param fit |
| χ²/dof, static Δχ², H window | ≤ 3, ≥ 25, [−9, 5] | fit quality, star veto, photometric plausibility |

## Validation state (2026-07-19)

- Offline: synthetic 600 AU mover recovered with fitted d within 10%;
  static source and 120 AU decoy rejected; 2-night object rejected;
  **SR-6 calibration harness ≥ 90% for V ≤ 23 with ≥ 3 covered nights** on
  a dense synthetic cadence (all as tests, all green).
- Live: 4,509 real selection-A detections (single night, known asteroids)
  through DE421 → 0 candidates, as physics demands. A real candidate
  search needs selection-B data spanning ≥ 3 nights — the standing portal
  ask.

## Candidate handling

`rubin_watch/candidates/cand-*.json` (status `linked`). Before human
review: run the script-side vetting (SkyBoT/astcheck — to be wired in the
vet script when the first real candidate exists), fetch cutouts, check
archival precovery per design/05. **No MPC submission** until design/09
OQ-1 is resolved.
