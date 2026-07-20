# Rubin Watch

Hitching the planet9 workspace into the live Rubin/LSST data ecosystem: new
distant-object discoveries flow into our statistical machinery, Rubin's
accumulating sky coverage flows into the search-hull/viability maps, and the
alert stream is mined for the ultra-slow movers (d ≳ 200 AU) that sit at the
edge of Rubin's own linking pipeline — the regime where Planet Nine lives.

This directory will hold the operational ledgers and run reports; the system
is specified in `design/`. Nothing here runs yet — the design is the current
deliverable.

## Why now

LSST survey operations began 2026-06-30. The April 2026 early-data release
already put ~11,000 asteroids and ~380 TNOs through the MPC, including
**2025 LS2** (a = 523 AU, e = 0.917, q = 43.2 AU, ϖ = 336°) — a bona-fide
ETNO by our vetted-sample cuts whose longitude of perihelion lands ~76° off
the ϖ ≈ 52° cluster mean. Every object like it moves the clustering, gap,
and posterior statistics this workspace maintains. Our frozen tables go
stale in months at Rubin's discovery rate; this system keeps them current
and turns each discovery into a recomputed verdict.

## Document map

| Doc | Contents |
|---|---|
| [design/00-goals-and-requirements.md](design/00-goals-and-requirements.md) | Goals, non-goals, science requirements, success criteria |
| [design/01-architecture.md](design/01-architecture.md) | System architecture, dataflow, language split, repo layout |
| [design/02-layer1-mpc-watch.md](design/02-layer1-mpc-watch.md) | **Layer 1**: MPC/SBDB distant-object poller (Phase 1) |
| [design/03-layer2-coverage.md](design/03-layer2-coverage.md) | **Layer 2**: LSST coverage-to-date → search hull (Phase 2) |
| [design/04-layer3-broker-intake.md](design/04-layer3-broker-intake.md) | **Layer 3**: Fink/broker alert intake (Phase 3) |
| [design/05-layer4-slow-mover-search.md](design/05-layer4-slow-mover-search.md) | **Layer 4**: distance-hypothesis slow-mover linking (Phase 4) |
| [design/06-statistics-battery.md](design/06-statistics-battery.md) | What recomputes when a new object lands; delta reports |
| [design/07-schemas.md](design/07-schemas.md) | Ledger / candidate / coverage / alert-subset schemas |
| [design/08-operations.md](design/08-operations.md) | Runbooks, cadence, git/PR conventions, credentials, failure modes |
| [design/09-risks-open-questions.md](design/09-risks-open-questions.md) | Risk register, open questions, decision log |

## Phasing at a glance

| Phase | Layer | Needs | Value |
|---|---|---|---|
| 1 | MPC watch | nothing (public, verified live 2026-07-16) | most of the science, ~5% of the effort |
| 2 | Coverage map | Layer-1 plumbing | hull/viability stay honest as LSST eats the southern residual |
| 3 | Broker intake | **Fink registration (user action)** | raw material for Layer 4; coverage cross-check |
| 4 | Slow-mover search | Layers 2+3 | the one search nobody else is guaranteed to run |
| 5 | Streaming/SNAPS | broker APIs maturing | latency reduction only; not on critical path |

## Status

- 2026-07-19 (later) — **Phase 2 (coverage) implemented.** Store →
  `rubin_watch/lsst_coverage.json` (sparse, healpy-convention; Rust
  `ang2pix_ring` fixture-locked to healpy) → hull "LSST (to date)" entry
  with per-pixel linked depth (P(≥3) = 95%) and crowding exclusion →
  viability "to date" row. First real map: 1,187 deg² covered / 155 deg²
  linkable from one night; hull headline correctly unmoved. See
  `RUNBOOK-coverage.md`.
- 2026-07-19 — **First LSST Data Transfer landed** (Phase-3 volume spike
  complete). 45,185 SSO-tagged alert rows from one night, 100% extracted
  into the partitioned store, 0% contamination, all acceptance checks
  green — including a ledger cross-match: Rubin caught our tracked
  distant object **2014 JW80** (a = 138 AU) at z = 22.0. See
  `reports/rubin-watch/2026-07-19-lsst-transfer-spike.md`. (A first
  portal job failed in Fink's server-side date parsing — reported
  upstream; broker metadata showed their LSST transfers stalled for all
  users since 07-15.) Next: selection-B job + coverage-map builder.
- 2026-07-17 — **Fink access live** (Phase 3 opened). Registered both
  surveys with fink-client 11.0; topic inventories pulled via Kafka
  metadata; real alerts consumed from `fink_hostless_candidate_lsst`
  (LSST) and `fink_sso_ztf_candidates_ztf` (ZTF). Key fact: no LSST SSO
  livestream topic exists yet — Data Transfer is the LSST SSO path, as
  designed. See `RUNBOOK-fink.md` and
  `reports/rubin-watch/2026-07-17-fink-spike.md`. Next: first portal
  Data Transfer job (browser step) → volume spike.
- 2026-07-16 (later) — **Phase 1 implemented and seeded.**
  `crates/p9-rubin-watch` (`mpc-sync` subcommand) + `RUNBOOK-mpc.md`;
  ledger seeded with the live census (370 minor planets: 358 SBDB
  ∩/∪ 347 MPC-distant; comets excluded via `sb-kind=a` after their
  discoverer-name aliases poisoned the join). Inaugural delta report at
  `reports/rubin-watch/2026-07-16-mpc.md`: the cut-defined vetted ETNO
  sample is now n = 20 with R̄ = 0.175 (Rayleigh p = 0.55) vs the curated
  10-object baseline's R̄ = 0.649 (2.64σ) — the clustering question now has
  a live, reproducible dashboard. 11 gap-band and 61 Neptune-crossing-class
  curation candidates surfaced; posterior-staleness flag set (16-object
  sample diff). Fink scaffold at `scripts/rubin_watch/fink_pull.py`
  (blocked on registration).
- 2026-07-16 — design directory written. External endpoints verified: MPC
  `distant_extended.json.gz` (updates daily), JPL `sbdb_query.api`
  constraint queries, SBDB single-object lookups.

## Glossary

- **DIA / DIASource** — difference-image analysis; one detection on one
  visit's difference image. The atom of the alert stream.
- **ssObjectId** — Rubin's identifier linking a DIASource to a known solar
  system object; absence = "unassociated".
- **SSP** — Rubin's nightly Solar System Processing (tracklet-less linking,
  MPC submission, MPCORB/SSObject daily products).
- **Broker** — community service receiving the full alert stream (Fink,
  ALeRCE, ANTARES, Babamul, Lasair, Pitt-Google, AMPEL) or a filtered
  downstream (SNAPS, POI).
- **MPEC** — Minor Planet Electronic Circular; discovery/orbit announcements.
- **Tracklet** — intra-night pair/triplet of detections of one mover; P9-class
  objects move too little intra-night to form one (see design/05).
- **Stationary point** — epoch where a distant object's apparent motion
  crosses zero as Earth's reflex reverses; a linking hazard and a false-match
  source (see design/05).
- **ETNO** — extreme trans-Neptunian object; our vetted class is a ≥ 250 AU,
  q ≥ 40 AU (see design/02 for the full class table).
- **PPDB / RSP** — Prompt Products Database / Rubin Science Platform
  (data-rights gated; deliberately not on our critical path).
