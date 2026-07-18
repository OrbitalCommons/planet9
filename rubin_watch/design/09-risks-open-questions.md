# 09 — Risks, open questions, decisions

## Risk register

| ID | Risk | L×I | Mitigation |
|---|---|---|---|
| R-1 | Fink LSST SSO filtering less capable than advertised (server-side filters immature in 2026) | M×M | selection filters degrade gracefully to client-side in `fink_pull.py`; volume budget assumes worst case; ALeRCE as secondary tap |
| R-2 | MPC file format drift breaks the parser | M×L | fixture-tested parser; unknown fields ignored; exit-2 + ops issue on required-field loss; SBDB census is an independent path |
| R-3 | 1-opposition Rubin orbits swing (a–e degeneracy) and churn class tags | H×M | tier system is the design answer: provisional never enters frozen tables; class-change history events make churn visible instead of silent |
| R-4 | Template self-subtraction blinds the linker from year 2 | H×M | dipole-aware static veto + negative-flux intake (design/05); `template_epoch_flag` in coverage; year-1 window is clean and is where the calibration baseline gets built |
| R-5 | Stationary-point false links / missed links | M×M | ephemeris-aware windows, purity flags, Δχ² static veto, decoy injections (80–150 AU) in calibration |
| R-6 | Broker outage or alert-schema version bump | M×L | store is re-pullable by night; schema alarms fail loud; 120-day retention gives slack to re-pull |
| R-7 | Volume surprise (footprint intersects deep-drilling fields; DDF cadence ~10× normal) | M×L | 3× volume alarm; DDF field list excluded from selection B if they dominate (they are ~5 pointings, negligible P9 prior weight) |
| R-8 | Our exclusion claims outrun measured completeness | L×H | SR-6 hard-gates any exclusion statement on the injection surface; hull integration uses `depth_linked`, not single-visit depth |
| R-9 | Scope creep toward re-implementing SSP | M×M | non-goals list; the linker's domain statement (design/05 §honest scope) is the fence |
| R-10 | Statistical churn: weekly deltas re-litigate headline numbers noisily | M×M | regeneration thresholds (design/06) — 0.1σ / 0.5 pt gates; reports show movement without regenerating below threshold |

## Open questions

| ID | Question | Blocking | Next step / owner |
|---|---|---|---|
| OQ-1 | Policy on deriving astrometry from public alert packets and submitting discoveries to the MPC (attribution, duplicate-submission etiquette vs Rubin SSP) | Layer 4 external reporting only | ask on the Rubin community forum + MPC helpdesk *before* any candidate leaves the repo; until answered, candidates stop at internal reports |
| OQ-2 | Does the user hold Rubin data rights (US astronomer)? | nothing (design avoids RSP) | if yes, ConsDB visit tables upgrade Layer 2 (option 3 → option 1); check at leisure |
| OQ-3 | ~~Exact Fink LSST server-side filter surface~~ **Answered 2026-07-17**: no LSST SSO livestream topic exists (11 curated topics, SN/extragalactic + uniform-sample + hostless); LSST SSO data flows via Data Transfer jobs only, with `b_is_solar_system`/Light-SSO packets applied job-side. Remaining sub-question: which selection-B quality cuts the DT job builder accepts — measured in the first pull | — | first portal Data Transfer job |
| OQ-4 | SNAPS API availability (public-access paper Apr 2026, API "coming soon") | Phase 5 only | quarterly check; adopt for SSO enrichment when live |
| OQ-5 | Published machine-readable completed-visit feed (schedview/almanac) as Layer-2 source 2 | quality upgrade only | one-afternoon spike during Phase 2 |
| OQ-6 | Negative-flux DIA sources: does Fink transfer carry them with usable reliability scores? | R-4 mitigation from year 2 | test in Phase 3 with a known slow mover (e.g. Sedna: d ≈ 80 AU, in template, dipole signature expected) |
| OQ-7 | Alert astrometric error model (per-band, per-airmass) sufficiency for the 5-param fit χ² | Phase-4 thresholds | measure empirically from selection-A residuals of known TNOs — a free calibration set |

## Decision log

| ID | Decision | Rationale | Revisit if |
|---|---|---|---|
| D-1 | Python for broker I/O, Rust for all science | fink-client is Python-only; repo culture computes claims in Rust | a maintained Rust Kafka/Avro path for Fink appears |
| D-2 | MPC-first phasing (Layer 1 before any broker work) | most science value, zero credentials, feeds existing machinery | never — it's already the cheapest path |
| D-3 | No MPC submissions pending OQ-1 | policy unknown; reputational risk asymmetric | OQ-1 answered |
| D-4 | Parquet + DuckDB store, no database server | single-operator repo, reconstructible cache, SQL when needed | multi-user or >100 GB scale |
| D-5 | Coverage from alert metadata (+ visit-feed spike) | no data-rights dependency | OQ-2 yes AND OQ-5 inadequate |
| D-6 | HEALPix nside 64 coverage / nside 8 linker tiles | 0.84 deg² resolves survey edges below the hull's 3° grid; 53 deg² tiles bound linker memory | hull grid changes |
| D-7 | 90-day primary linking window, monthly step | B⊥ ≤ 2 AU keeps hypothesis grid ~480 steps; ≥3 nights realistic per WFD cadence | measured cadence in priority footprint differs materially |
| D-8 | Committed ledgers, cached bulk | reproducibility of claims vs repo bloat | ledger approaches ~10 MB (then split by year) |
| D-9 | Priority footprint auto-derived from the hull (polygon export) | single source of truth; footprint follows the science | never |

## Explicit dependencies on user action

1. ~~Fink registration~~ **Done 2026-07-17** (username matthew, group
   matthew_cosmicfrontier; connectivity verified both surveys). Remaining
   user-adjacent step: create the first portal Data Transfer job (browser
   login) — after that everything is scripted.
2. OQ-2 data-rights check (optional upgrade).
3. Frozen-table adoptions and candidate confirmations are human-reviewed
   PRs by design — the automation proposes, the human disposes.
