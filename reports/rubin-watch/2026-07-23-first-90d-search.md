---
date: 2026-07-23
kind: rubin-watch-linker
window: 2026-04-15..2026-07-14
selection: "NOT b_is_solar_system, Light static packet, diaSource.reliability > 0.9"
topics:
  - ftransfer_lsst_2026-07-23_90173   # batch 21, 3-night probe
  - ftransfer_lsst_2026-07-23_771892  # batch 22, 90-day window
rows_landed: 813196
nights: 45
tiles: 81
tiles_with_3_nights: 76
candidates: 0
---

# First real slow-mover search (90-day window)

## Data Transfer failure diagnosis (batches 8, 20)

Two consecutive selection-B jobs (2026-07-16..18) died server-side without
creating a topic — not even the `_schema` stump every live job writes
first. Reading the portal's Spark job source localized it: each requested
night is gated on the Fink statistics API, and a range with no archived
night makes the driver `sys.exit(1)` silently. The LSST archive stalled at
night **2026-07-14** while livestream alerts kept flowing, so "recent"
nights were unfetchable. `fink_submit_job.py` now pre-checks the range
against the same API and refuses doomed submissions. (Our selection-A
store data was, correspondingly, nights 07-13/14 — not 07-16/18.)

## The 3-night lesson (batch 21)

Selection-B over 07-12..14 (234,734 rows) linked to **nothing — and could
not have**: no nside-8 tile was visited on ≥ 3 distinct nights; Rubin's
pointings across those three nights were essentially disjoint at tile
scale. A 3-night pull is a vacuous null. The design's 90-day rolling
window (D-7) is not an optimization, it is the minimum instrument.

## The 90-day search (batch 22)

813,196 unassociated high-reliability detections over 45 archived nights
(2026-04-15..07-14; max night 159,621 < the 3×10⁵ alarm). Coverage
rebuilt: 1,933 reconstructed visits, 2,743 deg² covered, r-band linkable
865 deg² (was 185); hull "LSST (to date)" entry now 956 linkable pixels.
81 linker tiles exported, 76 spanning ≥ 3 detection nights.

**Result: 0 candidates in 81 tiles** — after one algorithmic amendment
learned live. The ungated run produced 4 "candidates" (fitted d 569–840
AU, χ²/dof 0.2–0.8) that were all kinematically impossible: fitted drift
terms 4–45× the bound-orbit mean motion at their fitted distances, and
implied H of −6.2 to −8.4 (brighter than Neptune, seen 3 times in 90
days). These are chance triplets — the 5-parameter model absorbs whatever
sky motion a random triplet demands into the drift term unless drift is
tied back to the distance it claims. The linker now gates fitted drift at
`DRIFT_MARGIN` (2) × 0.9856/d^1.5 °/day (`unphysical_drift_is_rejected`
locks it); injections draw truth within 1× the ceiling, so completeness
is unaffected (SR-6 test still ≥ 90%).

## Completeness on real cadence (SR-6, 40 injections/tile, seed 20260719)

| tile | visits | nights | live dets | efficiency |
|---|---|---|---|---|
| 0494 | 219 | 10 | 22,571 | 33/40 (0.82) |
| 0381 | 166 | 25 | 20,009 | 34/40 (0.85) |
| 0349 | 166 | 25 | 5,256 | 39/40 (0.97) |

Injected bodies at d ∈ [300, 1000] AU, V ∈ [20, 23.5], own motion within
the physical ceiling, detected through the per-band logistic efficiency on
each tile's REAL visit list. Fitted distances land within ~5% of truth at
high n_det (e.g. injected 970 → fitted 970; 368 → 368).

Miss anatomy on the worst tile (0494, 7/40): 2 injections drew V=20.3 at
900–960 AU (implied H ≈ −9.4, outside the H window — an artifact of the
harness drawing V independent of d, not a linker miss; harness refinement
noted), 2 were faint (V ≥ 23.3, ≤ 3 detections), and 3 landed in fields
whose visits concentrate in < 3 distinct nights (visits cluster by field
per night — a real coverage limitation the efficiency surface must and
does charge against us). Tile 0349, with 25 well-spread nights, reaches
0.97.

## Honest scope

The exclusion statement this run supports: within the r-band linkable 865
deg² (non-crowding), no body at 300–1000 AU brighter than the per-tile
completeness limit was present and moving during 2026-04-15..07-14, at the
measured per-tile efficiency. The per-tile efficiency surface is what the
hull consumes — no claim extends beyond it. Next cycles densify coverage
monthly; the archive stall (07-14) is the current upstream bound — if it
persists past ~1 week it goes to contact@fink-broker.org.
