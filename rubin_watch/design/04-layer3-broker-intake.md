# 04 — Layer 3: broker alert intake (Phase 3)

Purpose: land a filtered, schema-stable subset of Rubin alerts on local disk
so Layer 4 can link across nights and Layer 2 can derive coverage. Alerts
are world-public with no proprietary period; the practical tap is a broker.

## Broker selection

| Broker | Role here | Basis |
|---|---|---|
| **Fink** (primary) | Data Transfer (batch, by-night, server-side filtered) + later Kafka livestream | most SSO-developed full-stream broker; LSST API/portal gained SSO support 2026; registration = username via email, PLAINTEXT Kafka, `fink-client` pip package |
| **ALeRCE** (secondary) | per-object API cross-checks (ssObjectId, stamp-classifier movers) | working SSO notebooks; independent classifier for candidate sanity checks |
| **SNAPS** (watch) | future: dedicated SSO enrichment/outlier alerts | the one solar-system downstream broker; public-access paper Apr 2026, API "coming soon" — integrate when it ships (design/09 OQ-4) |
| Lasair | none | explicitly does not handle solar-system alerts |
| Pitt-Google/AMPEL/ANTARES/Babamul | none for now | no SSO advantage over Fink for our use |

**Registration (user action, blocks this phase):** request Fink Data
Transfer + livestream credentials with an exclosure address; store via
1Password (`op inject` pattern, design/08). Everything below assumes it.

## Intake modes

### Batch (primary, Phase 3)

Nightly-or-weekly pulls via Fink Data Transfer:

- **Selection A — SSO-tagged**: the `b_is_solar_system` block (Fink's
  SSO-associated alerts) over the full sky. Cheap, small, powers coverage
  cross-checks and keeps a local record of every *known*-object detection
  in our footprint (useful for injection-calibration and for the
  known-object veto in Layer 4).
- **Selection B — unassociated residue in the priority footprint**: alerts
  with no ssObjectId, no Fink stellar/variable crossmatch, reliability high,
  19.5 ≤ mag ≤ 24.0, inside the **priority footprint** = the hull's
  80%-of-ν-weighted-residual region (~4,900 deg²; RA 0–120° ∪ 330–360°,
  Dec −40..+25 core, exported from `search_hull.json` as a polygon set —
  regenerated whenever the hull regenerates) **plus** the DES-hole cells.
- Exact server-side filter capability is Fink's moving target; anything not
  expressible server-side is applied client-side in `fink_pull.py` before
  writing (the record budget below is post-filter).

### Streaming (Phase 5, explicitly deferred)

`fink-client` Kafka consumption of the same selections once Fink's LSST
real-time SSO filters stabilize. Buys latency we do not need: a 300–1000 AU
mover is findable a week later. Not on the critical path; revisit when
SNAPS/Fink announce stable SSO topics.

## Volume budget

Order-of-magnitude, to be re-measured in the first week of operation:
~10M alerts/night over ~9,600 deg² visited ⇒ ~10³/deg²/night raw. Fink's
non-SSO, non-variable, high-reliability residue historically prunes 1–2
orders of magnitude; our footprint intersects a fraction of any night's
visits. Budget: **≤ 10⁵ selection-B records/night worst case, ~10⁴
typical** (≈ 5–50 MB/night as parquet). Alarm (design/08) if a night
exceeds 3× budget — that signals a filter regression, not science.

## Record schema (subset — full spec in design/07)

Keep exactly what linking and vetting need; drop cutouts (fetch on demand
for candidates only):

`alert_id, dia_source_id, dia_object_id, ss_object_id?, mjd, ra, dec,
ra_err, dec_err, band, psf_mag, psf_mag_err, reliability, visit, detector,
fink_class?, ingested_at`

## Storage

- **Format**: parquet, partitioned `night=YYYYMMDD/tile=<healpix-nside8>/`
  (nside 8 ⇒ ~53 deg² tiles — matches linker work units).
- **Store root**: `~/.cache/p9-rubin-watch/alerts/` — *never committed*;
  reconstructible from Fink by night.
- **Query layer**: DuckDB directly over the parquet tree (no server, no
  schema migration machinery). Python writes, Rust reads (arrow/parquet
  crate) — interchange pinned by the schema doc.
- **Retention**: 120-day hot window (linker needs ≤ 90-day windows + slack);
  older partitions archived to a tarball or dropped (re-pullable).

## Known-object veto data

Layer 4 needs "is this position a known object?" beyond ssObjectId:
maintain a weekly ephemeris cache of all MPC distant + local-region objects
(SkyBoT cone queries via astroquery for candidate positions, plus
`astcheck` against the MPCORB file for bulk pre-filtering). Both are
existing, community-standard tools; wrapping them is `scripts/rubin_watch/`
territory, and every veto is recorded on the candidate (design/07).

## Acceptance tests

1. One pulled night lands schema-valid (validated against design/07) with
   volume in budget; DuckDB count == record count.
2. Selection-B contamination: SkyBoT spot-check of 200 random records →
   < 5% known objects (measures the ssObjectId+crossmatch veto quality).
3. Rust reader round-trip: parquet written by Python is read by the linker
   with identical values (bit-exact for f64 columns).
4. Replay: deleting a night's partition and re-pulling reproduces identical
   record counts (broker-side determinism check; if it fails, document —
   brokers may reprocess — and key the store on alert_id dedup instead).
