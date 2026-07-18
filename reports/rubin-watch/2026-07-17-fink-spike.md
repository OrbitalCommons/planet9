---
schema_version: 1
run_utc: "2026-07-18T01:14:00Z"
layer: broker
trigger: ["fink-registration"]
degraded: []
---

# Rubin watch — Fink connectivity spike (Phase 3 opening)

Credentials arrived 2026-07-17 (username `matthew`, group
`matthew_cosmicfrontier`; LSST server `kafka-lsst.fink-broker.org:24499`,
ZTF `kafka-ztf.fink-broker.org:24499`; fink-client ≥ 11 required — design/04
said ≥ 10, amended).

## What was verified live

1. **Registration**: `fink_client_register` for both surveys; configs in
   `~/.finkclient/` (not committed).
2. **Topic inventory via Kafka metadata** — LSST broker: 2,320 topics, of
   which 11 curated livestream topics (SN/extragalactic oriented,
   `fink_uniform_sample_lsst`, `fink_hostless_candidate_lsst`) and the rest
   `ftransfer_lsst_*` Data-Transfer job topics. **No LSST SSO livestream
   topic exists yet** — confirming design/09 OQ-3/R-1: the LSST
   solar-system path today is Data Transfer, exactly the mode Layer 3 was
   designed around. ZTF broker: `fink_sso_ztf_candidates_ztf` and
   `fink_sso_fink_candidates_ztf` are live.
3. **Livestream consumption, both surveys**:
   - ZTF SSO candidates delivered (e.g. ZTF26abgssgf, mag 18.06, emitted
     2026-07-13, consumed 2026-07-18 — inside the ~7-day queue).
   - LSST alerts delivered (`fink_hostless_candidate_lsst`, diaObjectId
     170591519677875016 et al.).

## Implications for the design

- Selection A/B batch pulls proceed via portal-created Data Transfer jobs
  (browser step) + `fink_datatransfer` consumption (scripted). The
  ~7-day topic lifetime means pulls must be consumed promptly once created.
- `fink_hostless_candidate_lsst` overlaps our selection-B concept (no host
  crossmatch): worth a standing subscription as a low-volume real-time
  side channel — a P9-class stationary transient would land there.
- The ZTF SSO topics give us a permanent live testbed for consumer/parquet
  code before LSST SSO volumes exist.

## Next actions

1. Create the first Data Transfer job in the portal (one night,
   `b_is_solar_system`, Light-SSO packets) → run the design/04 acceptance
   tests: volume vs budget, schema round-trip into Rust, contamination
   spot-check. (Browser step — user or supervised session.)
2. Implement the consumption/re-partition loop in
   `scripts/rubin_watch/fink_pull.py` against the ZTF SSO testbed first.
3. Store the credentials email in 1Password (working copies are local
   YAML only).
