# 02 — Layer 1: MPC watch (Phase 1)

The highest-value, zero-credential layer: a weekly poller that discovers
newly designated distant objects, classifies them into the sample classes
our analyses consume, tiers them by orbit quality, and hands deltas to the
statistics battery (design/06).

## Sources (all verified live 2026-07-16)

| Priority | Source | Endpoint | Role |
|---|---|---|---|
| Primary | MPC distant-object extended file | `https://www.minorplanetcenter.net/Extended_Files/distant_extended.json.gz` (Last-Modified updates daily; few-MB gzip, ~5–6k Centaurs/TNOs) | discovery diff: the census of designated distant objects |
| Secondary | JPL SBDB constraint query | `https://ssd-api.jpl.nasa.gov/sbdb_query.api` with `sb-cdata={"AND":["a\|GT\|150","q\|GT\|30"]}` etc. | cross-check census (JPL syncs from MPC within days); 82 objects match a>150∧q>30 today |
| Per-object | JPL SBDB single lookup | `https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=<desig>` (already integrated via `starfield::SbdbClient`) | full-precision elements, uncertainties, arc, opposition count, orbit-solution id |
| Tertiary | MPC Orbits API | `https://data.minorplanetcenter.net/api/get-orb` (single designation only) | element cross-check in `mpc_orb` format when SBDB and MPC file disagree |
| Announce | Distant-object MPECs | MPC MPEC index | earlier-than-file discovery notice; parse opportunistically, never required |

Caveats to encode: SBDB serves ~3 significant figures for weakly constrained
TNOs (documented in `p9_core::data::refresh`); the MPC file and SBDB can
disagree during the days after designation — the ledger stores both when they
differ beyond the refresh-guard tolerances and flags the row.

## Class table (SR-1)

An object can carry multiple tags. Cuts follow the crates that consume them.

| Tag | Cut | Consuming analyses |
|---|---|---|
| `etno_vetted_class` | a ≥ 250 AU ∧ q ≥ 40 AU | apsidal clustering (p9-2026), 2021-orbit R̄, hull sample vetting |
| `etno_napier_class` | a ≥ 230 AU ∧ 30 < q < 40 AU | napier-critique debiasing sample |
| `gap_band` | 50 ≤ q ≤ 75 AU | perihelion-gap statistic (p9-2021-perihelion-gap) |
| `gap_context` | q ≥ 40 AU ∧ 100 ≤ a < 250 AU | gap-population context table |
| `high_inclination` | i ≥ 60° | inclined-TNO analyses (p9-2016-inclined-tnos) |
| `retrograde` | i > 90° | Drac-class comparisons |
| `neptune_crossing_class` | a > 100 AU ∧ i < 40° ∧ q < 30 AU | p9-2024-neptune-crossing sample rules |
| `sedna_like` | q ≥ 60 AU | IWG/Sedna-population notes |
| `watch_only` | a ≥ 150 AU, none of the above | census tracking |

## Tier rules (SR-2)

| Tier | Criteria (all required) | Statistics treatment |
|---|---|---|
| `provisional` | designated, but: 1 opposition, or arc < 400 d, or SBDB condition code ≥ 6, or e ≥ 1 solution | battery runs *with and without*; never enters a frozen table; delta report shows both |
| `vetted` | ≥ 2 oppositions ∧ arc ≥ 400 d ∧ condition code ≤ 5 | eligible for frozen-table adoption via a human-reviewed PR (never auto-committed into `p9_core::data::*`) |
| `retired` | designation superseded/merged, or object dropped below all class cuts after orbit revision | kept in ledger with pointer to successor; excluded from battery |

Rationale: Rubin's first-year orbits ride the a–e fit degeneracy hard
(the same drift our sbdb-refresh guard documents on multi-decade objects —
worse at 1 opposition). Example to encode as a fixture: 2025 LS2 is
`provisional` today (arc 350 d, soln 2) despite its spectacular a = 523 AU.

## Poller algorithm

```
mpc-sync [--dry-run] [--ledger rubin_watch/etno_ledger.json]

1. fetch distant_extended.json.gz  → parse → filter by SR-1 superset cut
2. fetch sbdb_query census (same cuts) → union of designations
3. diff against ledger designations (aliases resolved first — see below)
4. for each NEW designation:
     a. SBDB single lookup → full elements, epoch, arc, n_opp, condition code
     b. classify (class table) ; tier (tier rules)
     c. append ledger record (schema design/07), provenance = {mpc_file_date,
        sbdb_soln, fetched_at}
5. for each KNOWN designation:
     a. refresh elements if MPC file Last-Modified newer than ledger record
     b. detect tier transitions (promotion/demotion) and class changes
        (orbit revisions move objects across cuts — record as history event)
6. piggyback: run the existing sbdb-refresh drift check across all frozen
   tables; append any drift to the run report
7. emit delta = {new[], promoted[], reclassified[], drifted[]}
8. if delta nonempty and not --dry-run:
     run statistics battery (design/06) → write reports/rubin-watch/
     YYYY-MM-DD-mpc.md → open issue "rubin-watch: <n> new / <m> promoted"
9. exit code 0 = clean/no-op, 1 = delta emitted, 2 = source failure
```

**Alias handling.** Provisional designations get replaced by numbered ones
(e.g. 2004 VN112 → (474640) Alicanto). The ledger keys on a stable internal
id; `aliases[]` accumulates every MPC/SBDB spelling; the differ resolves
through aliases before declaring anything "new". Mis-resolution here is the
main dedup risk — the fixture suite must include a rename.

**Idempotence.** Running twice against unchanged sources must produce no
ledger change and no report (byte-identical ledger, exit 0). This is the
first acceptance test.

## Failure handling

| Failure | Behavior |
|---|---|
| MPC file unreachable / malformed | fall back to SBDB census alone; mark run `degraded:mpc` in report header |
| SBDB unreachable | ledger diff against MPC file only, elements refresh skipped; `degraded:sbdb` |
| Both unreachable | exit 2; no ledger mutation; issue only if 3 consecutive failures (ops, design/08) |
| Element disagreement MPC↔SBDB beyond refresh tolerances | ledger stores both under `elements` and `elements_alt`; row flagged `discrepant`; battery uses SBDB (full precision) and the report notes the flag |
| MPC format change | parser is fixture-tested; unknown fields ignored, missing required fields → exit 2 with the raw record logged |

## Acceptance tests (offline, fixtures)

1. **Seed reproduction**: fixture census (frozen copy of the 2026-07-16
   SBDB query) → ledger with expected counts per class tag.
2. **Idempotence**: second run, unchanged fixtures → no-op.
3. **New-object path**: inject synthetic a=400/q=45 object → classified
   `etno_vetted_class`, tiered `provisional` (1 opp), delta report contains
   with/without clustering numbers.
4. **2025 LS2 record**: fixture of its SBDB response → a=523, e=0.917,
   q=43.2, ϖ=336.0°, tier `provisional`, tags {etno_vetted_class}.
5. **Rename**: fixture rename event → no duplicate, alias appended, no delta.
6. **Demotion**: orbit-revision fixture moving q 41→38 → class change
   etno_vetted → napier, history event recorded, battery re-run flagged.

## First live run (Phase-1 exit)

Seed the ledger from the real census, ingest the ~380-TNO Rubin batch as it
appears in the distant file, and produce the inaugural delta report: the
vetted-clustering battery with the post-Rubin sample, explicitly including
the 2025 LS2 with/without comparison (ϖ = 336° vs cluster mean ≈ 52°,
circular separation ≈ 76° ≈ 2.2 circular-σ — expected to *reduce* R̄ and κ;
quantify in the report, not here).
