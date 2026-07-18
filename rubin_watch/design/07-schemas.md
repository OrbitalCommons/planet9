# 07 — Schemas

All committed artifacts are JSON with a top-level `schema_version` (integer)
and `generated_utc`. Consumers reject versions they don't know (fail loud,
never guess). Bulk alert data is parquet with the column table below.
Migration rule: bumping `schema_version` requires a migration note in this
file and a fixture for the old version in the parser tests.

## `rubin_watch/etno_ledger.json` (Layer 1)

```jsonc
{
  "schema_version": 1,
  "generated_utc": "2026-07-23T04:00:00Z",
  "source_stamps": {                       // provenance of the last sync
    "mpc_distant_last_modified": "…",
    "sbdb_query_fetched": "…"
  },
  "objects": [
    {
      "id": "rw-0001",                     // stable internal key (never reused)
      "primary_desig": "2025 LS2",
      "aliases": ["2025 LS2"],             // grows on numbering/renames
      "tier": "provisional",               // provisional | vetted | retired
      "class_tags": ["etno_vetted_class"], // design/02 class table
      "discovery": { "date": "2025-06-…", "source": "rubin_ssp|other|unknown" },
      "elements": {                        // best current (SBDB full precision)
        "epoch_jd": 0.0, "a_au": 523.0, "e": 0.917, "q_au": 43.2,
        "i_deg": 12.6, "omega_deg": 192.0, "omega_big_deg": 144.0,
        "varpi_deg": 336.0,
        "source": "sbdb", "solution_id": "2", "fetched_utc": "…"
      },
      "elements_alt": null,                // populated when MPC↔SBDB discrepant
      "quality": { "arc_days": 350, "n_oppositions": 1,
                    "condition_code": null, "n_obs": null },
      "flags": [],                         // e.g. "discrepant"
      "history": [                         // append-only event log
        { "utc": "…", "event": "ingested", "detail": "mpc_distant 2026-07-16" }
        // "refreshed" | "promoted" | "demoted" | "reclassified" | "renamed" | "retired"
      ]
    }
  ]
}
```

Ledger invariants (tested): ids unique and stable; every alias maps to
exactly one id; `history` append-only; serialization is canonical (sorted
keys, fixed float formatting) so no-op runs are byte-identical.

## `rubin_watch/lsst_coverage.json` (Layer 2)

```jsonc
{
  "schema_version": 1,
  "generated_utc": "…",
  "healpix": { "nside": 64, "ordering": "ring" },
  "window": { "first_night": "2026-06-30", "last_night": "2026-07-31" },
  "source": "fink_alert_metadata",         // + optional visit-feed cross-check
  "bands": {
    "r": {
      "n_visits":            [ /* 49152 × u16 */ ],
      "last_visit_mjd":      [ /* 49152 × f32, 0 = never */ ],
      "depth_single_median": [ /* 49152 × f32, 0 = never */ ]
    }
    // g, i as available
  },
  "flags": {
    "template_epoch": [ /* 49152 × bool */ ],   // design/03 systematic 1
    "crowding":       [ /* 49152 × bool */ ]    // |b| < 10°
  }
}
```

Derived quantities (`depth_linked(pix, N)` etc.) are computed by consumers,
never stored (design/03).

## Alert-subset parquet (Layer 3, NOT committed)

Partitioning `night=YYYYMMDD/tile=<nside8 pixel>/part-*.parquet`:

| column | type | notes |
|---|---|---|
| `alert_id` | u64 | dedup key |
| `dia_source_id` | u64 | latest-wins on broker reprocess |
| `dia_object_id` | u64 | |
| `ss_object_id` | u64 nullable | null ⇒ unassociated (selection B) |
| `mjd` | f64 | midpoint TAI as served |
| `ra`, `dec` | f64 | degrees, ICRS |
| `ra_err`, `dec_err` | f32 | arcsec |
| `band` | utf8 (1 char) | |
| `psf_mag`, `psf_mag_err` | f32 | |
| `reliability` | f32 | broker/DIA real-bogus |
| `visit`, `detector` | u64 nullable, u16 nullable | coverage derivation (LSST only) |
| `ss_name` | utf8 nullable | associated object designation when the packet carries one (ZTF `ssnamenr`) |
| `object_id` | utf8 nullable | survey object key (ZTF objectId / LSST diaObjectId as string) |
| `fink_class` | utf8 nullable | crossmatch label at intake time |
| `topic`, `survey` | utf8 | provenance |

(Selection membership is implied by the pull that produced the partition —
topic/provenance columns replace the earlier `selection`/`ingested_at`
fields; file mtimes carry ingest time. Implemented and validated against
live ZTF SSO alerts 2026-07-17.)

## `rubin_watch/candidates/<id>.json` (Layer 4)

```jsonc
{
  "schema_version": 1,
  "candidate_id": "cand-2026-08-0001",
  "status": "auto_vetted",   // linked → auto_vetted → human_review → rejected|confirmed
  "fit": {                    // the 5-parameter parallactic solution
    "lambda0_deg": 0.0, "beta0_deg": 0.0, "inv_d_au": 0.00167,
    "dlam_dt_deg_day": 0.0, "dbeta_dt_deg_day": 0.0,
    "d_au": 600.0, "d_au_interval": [560, 640],
    "chi2": 0.0, "dof": 0, "arc_days": 0.0, "n_nights": 0
  },
  "members": [ { "alert_id": 0, "mjd": 0.0, "ra": 0.0, "dec": 0.0,
                 "band": "r", "psf_mag": 0.0 } ],
  "vetting": {                // every stage recorded, pass or fail
    "static_veto_delta_chi2": 0.0,
    "astcheck": { "run_utc": "…", "matches": [] },
    "skybot":   { "run_utc": "…", "matches": [] },
    "photometric_consistency": { "pass": true, "note": "…" },
    "cutouts": { "fetched": 0, "artifact_rejects": 0 },
    "precovery": [ { "survey": "DES", "mjd": 0.0, "predicted_mag": 0.0,
                     "depth": 23.8, "found": null } ]
  },
  "score": 0.0,
  "window": { "tile": 0, "mjd_start": 0.0, "mjd_end": 0.0,
              "stationary_flag": false },
  "history": [ { "utc": "…", "event": "linked" } ]
}
```

## Delta-report frontmatter (`reports/rubin-watch/*.md`)

```yaml
---
schema_version: 1
run_utc: …
layer: mpc | coverage | broker | linker
trigger: [ "new:2025 LS2", "promoted:…", "monthly" ]
degraded: []           # e.g. ["sbdb"] per design/02 failure table
battery: { B1: {...}, … }   # machine-readable copy of every table in the body
artifacts_regenerated: [ "figures/search_hull.json" ]
issues_opened: [ 402 ]
---
```

The frontmatter duplicates the body's numbers so downstream tooling never
parses prose.

## Priority-footprint polygon export

`rubin_watch/priority_footprint.json` — polygon set (RA/Dec vertex lists)
exported by `p9-search-hull` alongside `search_hull.json`, marking the
80%-of-ν-weighted-residual region + DES-hole cells. `schema_version`,
`derived_from` (hull git hash), vertex arrays. Layer 3's selection B and
Layer 4's tile list both consume this file, so the footprint updates
automatically whenever the hull moves.
