---
schema_version: 1
run_utc: "2026-07-19T21:30:00Z"
layer: broker
trigger: ["first-lsst-data-transfer"]
degraded: []
---

# Rubin watch — first LSST Data Transfer (Phase-3 volume spike)

Topic `ftransfer_lsst_2026-07-19_391190` (portal job: one night,
`b_is_solar_system`, Light-SSO packets). A first attempt
(`…_554964`, Batch ID 5) failed server-side in Fink's Spark date parsing
before topic creation and was reported to the Fink team; broker metadata
showed no completed LSST transfer by any user between 2026-07-15 and this
job, with two schema-only stumps on 07-17 — their backend was wobbly, not
just our job.

## Acceptance results (design/04)

| Check | Result |
|---|---|
| Volume | 13,214 Kafka messages → **45,185 alert rows** (~3.4 rows/message), 4.7 MB raw — well inside the 10⁵/night budget; no alarm |
| Extraction | 45,185 / 45,185 rows mapped (100%); the Light-SSO packet is **flat** (not diaSource-nested) — `extract_record` now handles both shapes, verified against `lsst.v11_1` |
| Integrity | all diaSourceIds unique; 45 healpix tiles across 2 UTC nights (the observing night straddles the boundary) |
| Contamination | 0% — every row carries `unpacked_primary_provisional_designation` (selection A is fully associated, as designed) |
| Photometry sanity | median 21.1 mag, faintest 24.4 (consistent with single-visit depth); median reliability ≈ 1.0; bands i/r/z/g |

## Ledger cross-match (the point of the whole system)

The night's detections include **2014 JW80** — a vetted object in our
distant-object ledger (a = 138.4 AU, q = 38.0 AU, i = 40.9°, e = 0.725) —
observed at z = 22.04, reliability 0.98, RA 247.18°, Dec −30.81°
(MJD 61234.224). End-to-end: Rubin shutter → Fink transfer → our store →
matched against our census, all machine-side.

## Measured baselines (into design/04)

- Selection-A (SSO-tagged) volume: **~45k rows/night** early-survey.
- Light-SSO packet columns of record: diaSourceId, ssObjectId, phaseAngle,
  diaDistanceRank, snr, psfFlux/Err, scienceFlux/Err, templateFlux/Err,
  band, midpointMjdTai, ra, dec, reliability, packed/unpacked designations,
  schema/broker/science versions (`lsst.v11_1`, fink 5.0rc0 / 8.52.0).

## Next

1. Selection-B job (negate `b_is_solar_system` + quality cuts over the
   priority footprint) — the unassociated-residue volume measurement.
2. Phase 2 coverage: this pull already carries (night, band, tile) — the
   coverage-map builder has its first real input.
3. Weekly transfer cadence per RUNBOOK-fink once selection-B volumes are
   known.
