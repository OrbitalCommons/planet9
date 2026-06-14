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

### 9b. p9-2016-cassini-ranging — favored true anomaly 135.5° (proxy) vs published 117.8°

The crate computes the real differential (tidal) acceleration a 10 M⊕ P9 on the Brown & Batygin
orbit (Fienga Table 1: a=700, e=0.6, i=30°, ω=150°, Ω=113°) exerts on Saturn vs the Sun, then
forms the Earth–Saturn range signature ρ(t) over Saturn's REAL Cassini-epoch arc (two-body
propagation of the DE440 J2000 state, 2004–2017) and reduces it to a pre-fit RMS (paper's Fig 2
blue) and a post-fit RMS after removing the quadratic the INPOP orbit fit absorbs (Fig 2 red).
The post-fit curve reproduces the paper's qualitative structure: two excluded lobes on the
perihelion-facing arc (matching the forbidden zones [−65°,85°] ∪ [−130°,−100°]) and a low far-side
gap. The mass scaling is exactly linear and the distance falloff is steep (~1/r³ tidal), both pinned.

The honest residual: a template-free perturbation calculation places the absolute post-fit minimum
near aphelion (the paper's "uncertainty zone", where P9 is too faint to constrain). Applying the
paper's Fig 6 detectability logic — restrict to positions whose pre-fit signal exceeds the INPOP
precision floor σ₀ ≈ 74.9 m, on the outbound post-perihelion arc — yields a favored ν ≈ 135.5°
(r ≈ 783 AU), within the test's ±30° tolerance of the published v = 117.8°₋₁₀⁺¹¹ and adjacent to the
quoted favored interval [108°, 129°]. The exact 117.8° additionally requires anti-correlating the
induced signature with the measured Cassini O−C residual template (the data-driven mean-anomaly fit),
which a template-free proxy cannot reproduce; the post-/pre-perihelion tie is broken by scanning the
outbound arc, documented at the function site.
- Pinned: `crates/p9-2016-cassini-ranging/src/perturbation.rs` (`favored_zone_near_published_true_anomaly`,
  `perihelion_arc_is_excluded_relative_to_favored_zone`, `postfit_curve_has_two_excluded_lobes`,
  mass/distance scaling tests).

### 9. p9-2021-orbit — full inference pipeline: Reduced scale validates machinery, not the posterior

The paper's models 1–4 are now implemented from scratch (issue #12): N-body simulation grid
over (m₉, a₉, e₉, i₉) (`sim_grid.rs`), conditional KDE likelihood of the 10 vetted ETNOs with
the paper's 5%→30% bandwidth law (`kde.rs`), Matérn-3/2 GP emulator with grid-searched
hyperparameters and closed-form LOO cross-validation (`gp.rs`), and a Goodman–Weare
affine-invariant ensemble sampler with R̂/ESS diagnostics (`mcmc.rs`), wired in `pipeline.rs`.

**What the Reduced scale achieves** (12-point grid × 200 particles × 10 Myr; seeds 2021/777;
measured 2026-06-12, release build, 64 cores, ~68 s wall / 12.9 CPU-min):

- Pipeline health: MCMC acceptance 0.369, R̂ = [1.033, 1.038, 1.041, 1.053],
  ESS = [688, 647, 680, 594] (32 walkers × 2000 steps), all posterior mass inside the
  uniform prior box.
- GP LOO relative error 0.239 (pinned < 0.35 in the `#[ignore]`d end-to-end test). The
  paper's < 5% applies to its 121-point grid; 12 points in a 4-D box cannot reach that.
- The grid-point log-likelihoods span only ~3.7 nats (−30.5 to −26.7): after 10 Myr the
  secular anti-aligned clustering the inference keys on is only partially developed, so the
  posterior is close to the prior. The published medians (6.2 M⊕, 380 AU, e = 0.21, 16°)
  lie inside the Reduced posterior's 95% intervals — asserted as a *loose containment
  sanity check*, explicitly not a reproduction (the intervals are nearly as wide as the
  prior box).
- The posterior medians consumed downstream remain the published ones via the `posterior`
  summary module (re-documented as the published-posterior representation).

**What Paper scale needs** (`PipelineConfig::paper()`, `#[ignore]`d): the 121-point grid
surrogate (the paper's manual grid is unpublished; we use a regular grid over the same
extent, 8 masses × 15 (a₉, e₉, i₉) combos + 1 central point) × 16,800 particles × 4 Gyr,
plus the paper's MCMC settings (100 walkers, 20,890 steps, burn 260, thin 42 → ≈49,100
samples). At the measured throughput (~265 ns per particle-step) that is ≈3,000 CPU-days
(~1.5 months wall on a 64-core machine, less as removal thins the disk) — not run here.
Documented physics reductions even at Paper scale: J+S+U as the orbit-averaged J2 field
with Neptune+P9 direct (dt = 3000 d) instead of all four giants direct at dt = 300 d, and
the P9 orientation angles (ϖ₉, Ω₉) profiled out of the KDE likelihood over a 15° scan
rather than sampled as MCMC dimensions. The cluster-scattering Fréchet prior is not
implemented (uniform prior only; needs the Batygin & Brown 2021b scattering suite).

### 10. p9-2019-selfgrav-disk — forced apse aligned (Δϖ = 0), not anti-aligned

Sefilian & Touma (2019) show a massive eccentric apsidally-aligned debris disk shepherds ETNOs
into apsidal confinement. The Laplace–Lagrange operator built here (test particle vs a set of
confocal eccentric rings, softened `b_{3/2}` coefficients) reproduces the libration island and
the linear-in-M_disk precession rate, but with a **uniform-eccentricity** disk the forced
equilibrium is apsidally *aligned* (Δϖ = 0), not the paper's anti-alignment. The anti-aligned
forcing requires the disk's eccentricity/density gradients (sign of the cross-term B), which the
minimal model omits. The shepherding-into-a-libration-island mechanism itself reproduces; the
sign is documented rather than asserted.

The ~10⁸–10⁹ yr secular period (paper Sec. 3) is recovered for a fiducial disk softening of ~0.2·a
(a dynamically warm disk's finite radial/vertical thickness); a razor-thin disk precesses ~10×
faster. The softening is a labelled disk-thickness parameter.
- Pinned: `crates/p9-2019-selfgrav-disk/src/secular.rs` (linear-mass scaling, period order of
  magnitude, libration vs shepherding); cross-checked against `p9_core::analysis::secular`
  numerical ring averaging in `cross_check.rs` (analytic A within 5%).
### 11. p9-2019-ossos-scattering — scattering/detached boundary and the P9 lift are analytic/secular, not N-body

Kaib et al. (2019, OSSOS XV) constrain P9 with the observed actively-scattering TNO population.
This crate reproduces the two pieces of physics from p9-core machinery rather than an N-body
suite, with documented reductions:

- **Scattering/detached boundary.** The divide is the Chirikov 2:j resonance-overlap criterion
  from `p9_core::analysis::resonance` (K = 1 root, `critical_perihelion`). The honest output is
  that the Neptune 2:j chain only overlaps (q_crit > 0) above a ≈ 226 AU; below that there is no
  chaos at any q. q_crit(a) rises monotonically, passing the paper's q ≈ 37 AU scattering-disk
  scale between a = 400 and 500 AU and continuing to climb (≈ 41 AU at 500, ≈ 53 AU at 800 AU).
  The published "q ≈ 37 AU" is a single scattering-disk definition, not the a-dependent boundary;
  it is kept as `published::SCATTERING_Q_BOUNDARY_AU` and used only as a scale cross-check.
  Pinned: `crates/p9-2019-ossos-scattering/src/boundary.rs` (q_crit ≡ p9-core to 1e-12,
  K(a, q_crit) = 1 to 1e-10, the 37 AU crossing between a = 400/500).

- **P9 perihelion lifting.** Modelled as a coplanar test-particle secular cycle: the conserved
  Hamiltonian is P9's numerical Gauss-ring term (`p9_core::analysis::secular`) plus the giants'
  axisymmetric J2 precession term (`p9_core::forces::j2_secular::combined_j2_jsu`). The
  single-perturber Gauss-ring Hamiltonian scales as GM₉, so its level-curve *amplitude* is
  P9-mass-independent (a real property of linear secular theory — only the precession *rate*
  scales with mass). To recover the OSSOS-relevant mass dependence the lift is **rate-limited**:
  the secular eccentricity-forcing rate |de/dt| ∝ GM₉ (centred finite difference of ∂H/∂ϖ) drives
  q up over a sculpting time, capped at the cycle amplitude. The baseline time
  `SECULAR_SCULPTING_DAYS = 50 Myr` is a documented model choice (a coherent scattering episode
  before a Neptune encounter resets the orbit), not a derived quantity; it places the nominal-P9
  lift in the rate-limited (un-saturated) regime so the detached fraction responds monotonically
  to P9's mass. At this scale a nominal 6.2 M⊕ P9 over-detaches the seeded large-a scattering
  population (detached fraction → ~1.0 vs 0 with no P9), which is the qualitative OSSOS XV result:
  such P9 configurations are disfavored because they would deplete the observed scattering
  population. Not asserted as a quantitative reproduction — no published scalar exists for the
  detached fraction at fixed time.
  Pinned: `crates/p9-2019-ossos-scattering/src/planet_nine.rs` (rate ∝ GM₉ to 1e-6; lift > 0 for
  real P9, ≈ 0 for none; heavier ≥ lighter), `src/population.rs` (P9 raises detached fraction by
  a non-zero amount > 0.5; heavier P9 detaches ≥ as many).
### 12. p9-2018-wise-search — at W1 (3.4 µm) a cold P9 is detected by *reflected* light, not thermal emission

Meisner et al. (2018) frame the WISE/NEOWISE shift-and-stack search as a thermal-infrared one.
Our forward model exposes a physics reality that the "thermal-IR" label obscures at W1: at the
band's 3.4 µm and a plausible Planet Nine effective temperature of ≈40 K, hν/kT ≈ 107, so the
blackbody thermal flux is ~10⁻⁵⁵ W/m²/Hz/sr — utterly negligible. The thermal channel only
overtakes reflected sunlight above ≈175 K (at 15 M⊕, 400 AU), far hotter than any Planet Nine.
So the W1 detectability this crate computes is, for a cold body, set by **reflected sunlight**;
the thermal term is retained and would dominate for an implausibly warm planet. (The real search's
W2 4.6 µm band is more thermally favorable; the shared survey table pins the deeper, bluer W1.)
- Consequence for the exclusion: the reflected-light W1 reach is ≈186 AU (5 M⊕) to ≈222 AU
  (18 M⊕), only clipping the near edge of the Brown & Batygin posterior (q ≈ 300 AU), so the
  computed exclusion fraction over that posterior is a few percent (≈2% at 6.2 M⊕, ≈7% at 18 M⊕)
  — far below the optical ZTF 56%. This is an honest finding, not a tuned number; it increases
  monotonically with assumed mass as required.
- Pinned: `crates/p9-2018-wise-search/src/thermal_model.rs` (Planck/reflected model, W1 zero
  point 309.54 Jy), `detectability.rs` (finite max-distance bisection), `exclusion.rs`.
### 13. p9-2025-planet-y — warp amplitude: sub-degree for most of the allowed box, "a few degrees" only at the Earth-mass/close corner

The forced-inclination warp in the 50–80 AU band is computed from the linear Laplace-plane
balance `tan(2 i_forced) = A_out sin(2 i_Y) / (A_in + A_out cos(2 i_Y))`, with `A_in` the
giant-planet quadrupole (reusing `p9_core::forces::j2_secular::effective_j2` for JSUN) and
`A_out` Planet Y's exterior ring torque (`¼ n (m_Y/M) α² b_{3/2}^{(1)}(α)`). Across the
published parameter box the warp amplitude (tilt rise from 50→80 AU) is:

| m_Y (M⊕) | a_Y (AU) | warp amplitude |
|----------|----------|----------------|
| 0.05     | 100      | 0.28°          |
| 0.3      | 150      | 0.13°          |
| 1.0      | 200      | 0.13°          |
| 1.0      | 100      | 3.6°           |

So the paper's "few-degree" warp is reached only at the **Earth-mass, 100 AU corner**; a
Mercury-mass object, or an Earth-mass one out at 200 AU, warps the band by only a few tenths
of a degree. This is the honest reason the data signal is *tentative*. It also means the
torque crossover (`A_in = A_out`) sits at ≈ 95–180 AU — *exterior* to the 50–80 AU band — so
within the band the warp is the small rising shoulder of the transition, not full alignment
with Planet Y's plane. The forced tilt is strictly bounded by i_Y. No free parameter was
tuned; the band, mass range, and a_Y range are the paper's labelled constants in
`laplace_plane::published`.
- Pinned: `crates/p9-2025-planet-y/src/laplace_plane.rs` —
  `test_warp_feature_in_50_80_band` (0.01–1.0° for 0.3 M⊕ @ 150 AU),
  `test_warp_amplitude_few_degrees_for_earth_mass` (2–6° for 1 M⊕ @ 100 AU),
  `test_warp_amplitude_scales_with_mass` (∝ m_Y), `test_outer_dominates_near_planet_y`
  (crossover near a_Y, not in-band), `test_no_perturber_plane_stays_flat`.

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
