# 00 — Goals and requirements

## Problem statement

The workspace's scientific claims rest on frozen observational inputs (the
vetted ETNO table, the perihelion-gap sample, the survey-coverage hull) and
on statistics computed from them (κ/R̄/Rayleigh clustering, gap p-values,
posterior sky maps, exclusion fractions). Rubin/LSST began survey operations
on 2026-06-30 and will roughly decuple the known TNO population, while its
own coverage progressively converts our "residual un-searched probability"
into either a detection or an exclusion. Without automation, every one of
our headline numbers decays toward stale.

Separately, Rubin's nightly Solar System Processing is optimized for the
bulk population. A Planet-Nine-class body (d ≈ 300–1000 AU) moves
3.4–11.2″/day at opposition and 0.08–0.26″ within a standard ~33-minute
visit pair — it forms **no intra-night tracklet** and sits at the slow tail
of SSP's linking sensitivity. It would, however, be a perfectly ordinary
V ≈ 20.5–23 unassociated DIASource on every visit. Nobody's pipeline is
guaranteed to connect those dots; ours can be.

## Goals

- **G1 — Statistical currency.** Any newly designated distant object
  (class table in design/02) is ingested, classified, tiered, and reflected
  in a recomputed statistics battery within **one week** of its MPC orbit
  posting, with a written delta report and a tracking issue.
- **G2 — Coverage currency.** The search-hull and viability products carry
  an "LSST (to date)" survey entry rebuilt **monthly** from Rubin's actual
  observed footprint and depth, so the residual-probability numbers quoted
  anywhere in this repo are never more than a month behind the sky.
- **G3 — Independent slow-mover search.** A distance-hypothesis linking
  search over unassociated alerts in our priority footprint (the ~5,000 deg²
  holding 80% of the ν-weighted residual; RA 0–120°, Dec −40..+25 core),
  sensitive to d ∈ [300, 1000] AU and V ≤ 23.5, with a measured injection/
  recovery efficiency — not an assumed one.
- **G4 — Provenance and honesty.** Every ingested element set carries its
  source, epoch, arc, and uncertainty; every statistic distinguishes
  *provisional* from *vetted* inputs; every conditional prior stays labeled.
  Same culture as REPRODUCTION_NOTES — the automation must not launder
  1-opposition orbits into headline claims.

## Non-goals

- **Not** a re-implementation of Rubin SSP, heliolinc, or difference
  imaging. We consume DIASources; we never touch pixels except candidate
  cutouts.
- **Not** a full-stream consumer. ~10M alerts/night is Fink's job; we take
  filtered subsets (design/04 budgets ≤ ~10⁵ records/night worst case).
- **Not** a transient/variable science system. Anything with a stellar/host
  crossmatch or an ssObjectId is dropped at intake (Layer 4 keeps
  unassociated-only).
- **Not** an MPC submitter — not until the policy question (OQ-1,
  design/09) about submitting astrometry derived from public alert packets
  is resolved. Candidates stop at a vetted internal report.
- **Not** dependent on Rubin data rights. The design runs entirely on
  world-public products (alerts, MPC, SBDB). The RSP/PPDB path is an
  optimization documented in design/09, never a dependency.

## Science requirements

| ID | Requirement | Driver |
|---|---|---|
| SR-1 | Track objects with a ≥ 150 AU or q ≥ 40 AU or (q > 30 AU and a ≥ 100 AU) or i ≥ 60° | superset of every sample any crate consumes (design/02 class table) |
| SR-2 | Distinguish provisional (1-opposition) from vetted (multi-opposition) orbits; statistics run on both, reported separately | the a–e fit degeneracy documented in `p9_core::data::refresh`; Rubin's first-year orbits will swing |
| SR-3 | Clustering deltas reported under all three labeled selection stand-ins (`p9_core::analysis::selection::VarpiSelection`) | the verdict is selection-model dependent (core verdict-sensitivity test) |
| SR-4 | Slow-mover search covers apparent rates 2–14″/day near opposition and is ephemeris-aware through stationary points | d ∈ [300, 1000] AU ⇒ 3.4–11.2″/day at opposition, → 0 at stationary epochs |
| SR-5 | Slow-mover search depth V ≤ 23.5 (single-visit), ≥ 3 nights in a 90-day window | posterior residual is V ≈ 20.5–23 at d 400–900 AU; 3 nights over-determines the 5-parameter fit |
| SR-6 | Injection/recovery efficiency ≥ 90% for synthetic P9 draws with V ≤ 23 and ≥ 3 covered nights | claimed exclusions from non-detection must be backed by measured completeness |
| SR-7 | Coverage map resolution ≤ 1 deg² (HEALPix nside 64), per band, with visit counts and depth | matches the hull's 3° grid without aliasing survey edges |
| SR-8 | All tests offline (fixtures); every network touch behind an explicit binary/runbook invocation | repo-wide rule; the gate must never hit the network |

## Success criteria (per phase)

- **Phase 1** — first ingest reproduces the current SBDB census (82 objects
  with a > 150 AU, q > 30 AU as of 2026-07-16), classifies 2025 LS2 into the
  vetted-candidate class, and emits a delta report with the clustering
  battery run with-and-without it. Second run on unchanged inputs is a
  no-op (idempotence).
- **Phase 2** — `search_hull.json` regenerates with an "LSST (to date)"
  entry; the excluded fraction is monotonically non-decreasing across
  monthly runs; a synthetic full-southern-coverage input reproduces the
  hull's existing Rubin-forecast numbers within tolerance.
- **Phase 3** — one night of Fink Data-Transfer output lands in the local
  store with schema-validated records; volume within budget; known-object
  contamination after filtering < 5% (spot-checked against SkyBoT).
- **Phase 4** — SR-6 met on injections; false-candidate rate ≤ ~10 per
  month reaching human review; zero known objects surviving the vetting
  chain in a blind test month.
