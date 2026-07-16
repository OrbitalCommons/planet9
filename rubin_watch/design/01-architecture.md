# 01 — Architecture

## System diagram

```
  EXTERNAL (all world-public unless noted)
  ┌────────────────┐  ┌──────────────────┐  ┌───────────────────┐
  │ MPC            │  │ JPL SBDB          │  │ Fink broker        │
  │ distant_ext-   │  │ sbdb.api (single) │  │ Data Transfer +    │
  │ ended.json.gz  │  │ sbdb_query.api    │  │ Kafka livestream   │
  │ (daily) MPECs  │  │ (constraints)     │  │ (registration req) │
  └───────┬────────┘  └────────┬──────────┘  └─────────┬─────────┘
          │ weekly             │ per-object            │ nightly/monthly
          ▼                    ▼                       ▼
  ┌───────────────────────────────────┐   ┌───────────────────────────┐
  │ LAYER 1: mpc-watch (Rust)         │   │ LAYER 3: broker intake     │
  │ poll → diff → classify → tier     │   │ (Python: fink-client)      │
  │ crates/p9-rubin-watch             │   │ scripts/rubin_watch/       │
  └──────┬───────────────┬────────────┘   └──────┬──────────┬─────────┘
         │               │                        │          │
         ▼               │                        ▼          ▼
  rubin_watch/           │              local alert store   LAYER 2:
  etno_ledger.json       │              (parquet+DuckDB,    coverage map
  (committed)            │               NOT committed)     rubin_watch/
         │               │                        │         lsst_coverage.json
         ▼               ▼                        ▼         (committed)
  ┌──────────────────────────────┐   ┌──────────────────────────┐   │
  │ STATISTICS BATTERY (Rust)    │   │ LAYER 4: slow-mover      │   │
  │ clustering / gap / crossing  │   │ linker (Rust)            │   │
  │ via existing paper crates    │   │ distance-hypothesis fit  │   │
  └──────┬───────────────────────┘   └──────┬───────────────────┘   │
         │                                  │                       │
         ▼                                  ▼                       ▼
  reports/rubin-watch/*.md          rubin_watch/candidates/   p9-search-hull /
  + GitHub issue per delta          *.json + vetting chain    p9-viability regen
```

## Layering principle

Each layer is independently useful and independently shippable, ordered by
value-per-effort. Layer 1 needs no credentials and delivers most of the
science (goal G1). Layer 2 keeps the maps honest (G2). Layers 3–4 are the
custom search (G3) and are the only components with an external-credential
dependency (Fink) or nontrivial storage.

## Language and dependency split

**Decision D-1**: Python owns broker I/O, Rust owns science.

- **Rust** (`crates/p9-rubin-watch`, one binary, subcommands per layer):
  polling MPC/SBDB (reuses `starfield::sbdb::SbdbClient` via
  `p9_core::data::refresh`), ledger management, classification/tiering,
  statistics battery orchestration (calls into the existing paper crates as
  path dependencies), slow-mover linking (reuses `p9_core::coords::observer`
  Earth ephemeris + `candidate_pair` machinery from the IRAS/AKARI work),
  coverage-map ingestion for the hull.
- **Python** (`scripts/rubin_watch/`): `fink-client` consumption (it is
  Python-only), parquet writing, DuckDB spot queries, plotting. Mirrors the
  existing `scripts/plot_*.py` convention: Python at the I/O and viz edges,
  Rust for every number that ends up in a claim.

Interchange format between the two: parquet (alert subsets) and JSON
(everything else), schemas pinned in design/07. No Python computes a
statistic that lands in a report.

## Repo layout additions

```
crates/p9-rubin-watch/          # Rust binary: mpc-sync | coverage | link | battery
scripts/rubin_watch/            # fink_pull.py, store.py, plot_candidates.py
rubin_watch/                    # committed state (this dir)
  README.md, design/            #   docs (this design)
  etno_ledger.json              #   Layer-1 ledger        (committed)
  lsst_coverage.json            #   Layer-2 map           (committed)
  candidates/                   #   Layer-4 vetted output (committed)
reports/rubin-watch/            # delta reports, one per run (committed)
~/.cache/p9-rubin-watch/        # alert store, cutouts     (NOT committed)
```

Committed-vs-cache rule: anything a statistic depends on is committed
(ledger, coverage, candidates); bulk intermediates (alert parquet, cutouts)
are cache, reconstructible from the brokers, and capped (design/08).

## Dataflow triggers

| Event | Detected by | Triggers |
|---|---|---|
| New/changed distant object at MPC | Layer 1 weekly diff | classify → tier → statistics battery → delta report + issue |
| Provisional object gains 2nd opposition | Layer 1 (arc/opposition fields) | tier promotion → battery re-run flagged "promotion" |
| Month boundary | cron/manual | Layer 2 coverage rebuild → hull + viability regen → delta if excluded% moved ≥ 0.5 pt |
| New Fink night pulled | Layer 3 | store append → Layer 4 incremental link over affected tiles |
| Linker candidate survives auto-vetting | Layer 4 | candidate record + cutout fetch + human-review issue |
| Frozen-table drift (existing sbdb-refresh guard) | Layer 1 piggyback | drift report appended to delta |

## Design tenets (inherited from the workspace)

1. Offline determinism: every test runs from fixtures; network only in
   explicitly invoked binaries (SR-8).
2. One source of truth: survey geometry/depths come from `p9_core::analysis::
   surveys`; photometry from `p9_core::analysis::photometry`; no local copies.
3. Honest tiers and labels: provisional vs vetted, conditional priors named
   at every use (SR-2, SR-3).
4. Append-only ledgers with full provenance; regeneration is always possible
   from source APIs + git history.
5. Everything lands via branch + PR + the standard full-suite gate.
