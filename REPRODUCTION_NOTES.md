# Reproduction Notes — Residuals vs Published Values

Companion to [MODELING_REVIEW.md](MODELING_REVIEW.md). After the modeling fixes, every crate
computes its headline quantities from real (reduced-scale, seeded) calculations and pins them
against the published numbers in regression tests. This document records the places where the
honest computation does **not** exactly match the paper, what the residual is, the best current
explanation, and where the assertion lives. None of these are papered over in code — each test
asserts the computed value within a stated tolerance and documents the gap.

A residual here means one of three things: (1) our reduced-scale or simplified model genuinely
can't capture an effect the paper's full pipeline does, (2) our input data differs slightly from
the paper's (sample membership, element epochs), or (3) the paper's own stated parameters don't
reproduce its stated result (a genuine reproduction discrepancy worth flagging upstream).

---

## Disagreements

### 1. p9-2024-oort-selfgrav — J2 energy share: computed ~7% vs paper's "<2%" at q = 75 AU

The most interesting residual in the repo: using **the paper's own published fit parameters**
(Batygin & Nesvorný 2024) and the corrected planetary quadrupole, the J2 share of the secular
Hamiltonian evaluates to ~7% at q = 75 AU, not the paper's stated <2%. The published *trend*
(rising toward ~25% at q = 30 AU) reproduces in shape. Either the paper used a different
normalization for the disk potential, or its <2% quote refers to a different decomposition.
- Pinned: `crates/p9-2024-oort-selfgrav/src/hamiltonian.rs:152` (NOTE comment), test at `:209-222`
  asserts the computed ~7% and the monotone 75→50→30 AU trend.

### 2. p9-2024-panstarrs — combined exclusion: computed 0.807 vs published 0.785

The per-orbit OR over one shared seeded reference population, run through all three survey
models, gives ZTF 0.608 / DES-unique 0.024 / PS1-unique 0.175 / combined 0.807 vs the published
0.564 / 0.050 / 0.171 / 0.785. The ZTF over-exclusion (see #3) propagates: a hotter ZTF model
absorbs orbits DES would otherwise uniquely exclude, which is why DES-unique under-shoots while
the combined total over-shoots. The published-table arithmetic (0.564 + 0.050 + 0.171 = 0.785)
is preserved exactly via the shared `p9_core::analysis::surveys` table as a cross-check.
- Pinned: `crates/p9-2024-panstarrs/src/combined_exclusion.rs:227-252` (table identity to 1e-10;
  population-based values within stated tolerances).

### 3. p9-2021-ztf — exclusion fraction: computed 0.609 vs published 0.564

Our per-orbit model (sky position from elements, δ > −30° footprint, logistic
magnitude efficiency, linking factor) over-excludes by ~4.5 points. The paper injected synthetic
P9s into the actual ZTF epoch/pointing/depth stream; our static footprint + global efficiency is
optimistic about cadence losses. Untuned — no free parameter was adjusted to close the gap.
- Pinned: `crates/p9-2021-ztf/src/exclusion.rs:74-98` (asserted within 0.08 of 0.564, fed by a
  seeded p9-2021-orbit reference population).

### 4. p9-2021-orbit — clustering confidence: computed ≈0.989 vs published 0.996

Two known input differences: the vetted `p9_core::data::etno` table carries 10 objects vs the
paper's 11-object 2021 sample, and the survey-bias null uses our simplified longitude-suppression
model rather than the paper's per-survey pointing histories.
- Pinned: `crates/p9-2021-orbit/src/statistical_measures.rs:184-185` (within 0.01 of 0.996).

### 5. p9-2016-constraints — observed 6-KBO R̄: 0.84 from catalog elements vs ~0.96 implied by the paper's quote

The paper's "ϖ = 71° ± 16°" phrasing implies R̄ ≈ 0.96, but computing R̄ directly from the six
objects' catalog elements gives 0.84 (circular σ ≈ 34°). The ±16° quote appears to describe the
tightly clustered core rather than the full six-object circular dispersion. Acceptance tests use
the data-derived 0.84.
- Pinned: `crates/p9-2016-constraints/src/clustering_metric.rs:83,213-216`.

### 6. p9-2021-oort-cloud — 67%/88% ϖ-confinement fractions: deliberately NOT asserted

The reduced-scale secular model (P9 as a static dwell-time-weighted Gauss ring) produces the
qualitative result — the P9 run confines the distant belt more than the P9-free control, robust
across 5 seeds — but the injected-IOC subset trends *aligned* rather than the paper's 67%
*anti-aligned*. Reproducing the sign of the IOC confinement requires the apsidal precession of
P9 and the giants (full N-body), which the static-ring model cannot represent. The published
numbers are documented as out of reach of this model rather than asserted.
- Documented: `crates/p9-2021-oort-cloud/src/injection_simulation.rs:437-443`.

### 7. p9-2025-iras-akari — expected chance pairs: λ ≈ 19 vs the paper's 13 observed

The post-cut Poisson estimate (conditioned on post-cut source counts with a galactic-latitude-
dependent IRAS density model) lands at the right order of magnitude but ~50% high. Sensitive to
the assumed post-cut catalog sizes and the disk/isotropic density split, both documented
assumptions. The false-alarm probability for ≥1 surviving pair, P = 1 − e^(−λ), is the key new
statistic and is reported alongside.
- Documented: `crates/p9-2025-iras-akari/src/candidate_search.rs:202-241`.

### 8. p9-2025-perturbation — Hansen asymptotic prefactors corrected against quadrature

Not a paper residual but a deviation from formulas the original code claimed: validating the
asymptotic forms against the convergence-controlled numerical Hansen coefficients forced a
recalibration of X^(−4,1) from 1/8 to 0.71 (the original prefactor was ~5.7× off) and X^(−4,3)
from 1/12 to 1/9.5. If those literals trace to a published table, the source should be re-checked.
- See: `crates/p9-2025-perturbation/src/hansen.rs` (cross-validation tests against
  `p9_core::analysis::hansen::hansen_coefficient`).

---

## Agreements within tolerance (for completeness)

| Quantity | Computed | Published | Pinned at |
|---|---|---|---|
| DES recovery fraction | 0.877 | 0.87 | `p9-2022-des/src/recovery_analysis.rs` (±0.04, seed 2022, n=6000) |
| DES per-color-model recoveries | within 0.052 | Table 1 | `p9-2022-des/src/color_models.rs` (tolerance 0.12) |
| Neptune-crossing ζ null moments | −7.38 / σ 1.79 | −7.2 / 1.8 | `p9-2024-neptune-crossing/src/hypothesis_test.rs` |
| B&B 2019 combined clustering p | 0.0026 (MC σ ≈ 4e-4) | 0.002 | `p9-2019-clustering/src/clustering_analysis.rs:273-283` |
| Brown 2017 combined bias p | ≈2.3e-4 | 2.5e-4 | `p9-2017-bias` `#[ignore]`d 10⁶-iteration test (within 3×) |
| B&B 2016 joint significance | order 7e-5 band | 0.007% | `p9-2016-evidence` `#[ignore]`d 10⁷-trial MC (5e-6–1e-3 band) |
| F5 resonance implied-a₉ peak | ~661 AU | 660 AU | `p9-2018-resonance/src/probability_analysis.rs` |
| Critical period ratio (2017 dynamics) | in 0.1–0.15 | 0.1–0.15 | `p9-2017-dynamics/src/resonance.rs` |
| BMN21 q_crit(500 AU) | 41.4 AU | 41.4 AU | `p9-2021-stability` (exact formula identity, K(a, q_crit) = 1) |
| Review-paper critical a (5 M⊕, 500 AU, e=0.25) | 250 AU | 200–300 AU band | `p9-2019-review/src/parameter_survey.rs` |
| IRAS/AKARI candidate 47.46′ → distance | ~694 AU (DE421 ephemeris), ~675 AU (circular fallback) | 500–700 AU | `p9-2025-iras-akari/src/orbital_constraints.rs` |
| vZLK evolutionary timescale (nominal IOC) | > 4.5 Gyr | ≫ 4.5 Gyr | `p9-2024-oort-selfgrav/src/vzlk.rs` |

---

## Tuned/assumed parameters (documented, not derived)

These are free parameters chosen to match an observable, or assumptions standing in for data we
do not ingest. Each is commented at the definition site; listed here so nobody mistakes them for
derived quantities:

- `p9-2022-des` `night_usability = 0.29` (`survey_model.rs:162`) — the single calibration knob
  scaling per-night completeness so end-to-end recovery matches the paper; stands in for
  weather/chip-gap/masking losses we don't model per-exposure.
- `p9-2017-bias` galactic-plane suppression (center crossing σ=40°/floor 0.02, anticenter
  σ=15°/floor 0.5) and north/south asymmetry — tuned so the combined bias-corrected p lands near
  the published 2.5e-4; the real quantity requires actual survey pointing histories.
- `p9-2019-clustering` / `p9-2021-orbit` survey-bias null: longitude suppression near
  λ ≈ 95°/275° + δ > −30° coverage — same caveat.
- `p9-2025-clustering` selection weighting w(ϖ) = 1 + 0.3cos(ϖ − 60°) in the bias-resampling MC.
- `p9-2021-oort-cloud` `IOC_TRAPPING_FRACTION_ASSUMED` (Brasser et al. 2006 provenance) —
  replaces the previously invented √ρ scaling.
- `p9-2021-orbit` posterior emulator: (a, e) sampled jointly with an assumed ρ = 0.85
  correlation; the crate is documented as a posterior-summary emulator, not an MCMC reproduction.
- `p9-2024-neptune-crossing` discovery distances: **now per-object data, not a heuristic.**
  `FIRST_OBS_DISTANCES` (`hypothesis_test.rs`) stores each object's heliocentric distance at its
  SBDB `first_obs` epoch (first observation of the fitted arc; equals the discovery date except
  where precovery was later linked), computed 2026-06-12 by propagating the object's own SBDB
  osculating elements (`p9_core::data::refresh::first_obs_distance`). The fetched values span
  1.0q–1.7q (median ≈ 1.2q), so the old r ≈ 1.15q approximation was slightly low; it remains as
  the documented fallback (`approximate_discovery_distances`) for objects absent from the table.
  The table is offline/const; `--features sbdb-refresh` + the `#[ignore]`d live test recompute
  and pin it. The same feature gates `p9_core::data::refresh`, which diffs the frozen element
  tables (`data::etno::BROWN_2017_SAMPLE`, the 17-object neptune-crossing sample) against live
  SBDB lookups with documented per-element tolerances — the guard against the fabricated-table
  failure mode (MODELING_REVIEW.md §5 N1). Live check 2026-06-12: the neptune-crossing table
  matches JPL exactly; the Brown-2017 table intentionally pins paper-epoch solutions and shows
  documented a–e fit-degeneracy drift vs current JPL (see `etno.rs` module docs).
- Circular-Earth observer model — **now a fallback, not the primary geometry**. The survey
  models (`p9-2025-iras-akari` proper motion/parallax, `p9-2022-des`/`p9-2024-panstarrs`
  per-night apparent positions) take an explicit `p9_core::coords::observer::EarthProvider`:
  the primary `EphemerisEarth` uses real DE421 Earth states (via starfield) with light-time
  retardation (~3–4 days at 500–700 AU) and stellar aberration; `CircularEarth` (mean-longitude
  1 AU circle) remains as the documented analytic fallback for kernel-less machines and
  speed-critical MC loops. Kernel-gated tests pin the two within the documented error
  (≤ ~1′ for the IRAS/AKARI separations, < 0.3° of apparent position for DES/PS1 footprint
  work). Survey epochs are starfield `Time`s: IRAS mid-epoch 1983-06-24 and AKARI FIS
  mid-epoch 2006-12-31 (real window midpoints, refining the nominal 1983.5/2006.5; the
  mid-to-mid baseline becomes 23.52 yr instead of 23.0), DES night grid anchored at
  2013-08-31, PS1 epoch grid anchored at 2009-06-02.

  Regression shifts from the ephemeris path (documented magnitudes):
  - IRAS/AKARI implied distances move up ~3.4%: d(69.6′) 511 → 528 AU, d(47.46′) 675 → 694 AU,
    d(42′) 733 → 758 AU (longer real baseline + the 16-month AKARI window); separation windows
    widen accordingly (max at 500 AU 71.7′ → 75.3′). All pins stay inside the paper's
    500–700(+) AU interpretation; both model paths are asserted.
  - DES recovery and the combined ZTF/DES/PS1 exclusion are insensitive: at the seeded
    regression scale the ephemeris path reproduces the analytic values to ≤ 1e-4
    (recovery 0.8769 on both paths; combined exclusion 0.8090 vs 0.8091), since the
    real-vs-circular Earth difference is < 0.3° of apparent position against multi-degree
    footprint structure.

## Known numerical issues

- ~~`p9-core` `kepler_drift` (universal-variable iteration): convergence-plateau debug
  assertion can fire at extreme states (reproduced at a = 800 AU, e = 0.95, M ≈ 5.97,
  dt = 4e5 d)~~ — resolved. Root cause: with multi-orbit dt the universal Kepler residual is
  assembled from terms ~ √gm·|dt| whose rounding noise, divided by f′ = r (small near
  perihelion at high e), exceeds the 1e-15·s step tolerance, so safeguarded Newton dithers in
  a few-ulp limit cycle until the iteration budget runs out. Fixed in
  `p9-core/src/integrator/kepler_step.rs` by (1) reducing dt modulo the orbital period for
  bound orbits, bracketing the universal anomaly within one revolution; (2) capping hyperbolic
  brackets where the Stumpff cosh would overflow and forcing bisection while the residual is
  ≫ its target scale (Newton only creeps additively in the exponential tail); (3) accepting
  convergence when the residual reaches its rounding-noise floor, the bracket collapses to
  relative machine precision, or the iterate enters a bit-exact period-2 cycle.
  Cross-validated against starfield's SPICE prop2b port in
  `p9-core/tests/starfield_oracle.rs` (which also pins the element conversions to
  `starfield::elementslib` as a permanent oracle). p9-2021-oort-cloud's element-space drift
  was kept — it is exact two-body propagation and the natural representation for that secular
  model, not a workaround.
- ~~p9-core still lacks an ecliptic↔equatorial transform~~ — resolved: `p9_core::coords::sky`
  wraps starfield's framelib (ECLIPJ2000/GALACTIC SPICE matrices); the local conversions in
  `p9-2022-des/src/sky.rs`, `p9-2021-ztf/src/sky.rs`, `p9-2025-iras-akari` and the hand-rolled
  galactic pole in `p9-core/src/forces/galactic_tide.rs` were deleted in favor of it.
- No shared simulation snapshot driver in p9-core; the 2016-crate run loops were harmonized in
  place but remain three copies.
