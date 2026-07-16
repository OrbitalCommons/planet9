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

### 2. p9-2024-panstarrs — combined exclusion: computed 0.713 vs published 0.785

The per-orbit OR over one shared seeded reference population, run through all three survey
models, gives ZTF 0.466 / DES-unique 0.036 / PS1-unique 0.211 / combined 0.713 (ephemeris path
0.466 / 0.036 / 0.206 / 0.708) vs the published 0.564 / 0.050 / 0.171 / 0.785. The population
now follows BB21's stated catalog assumptions (r₉ = (m₉/3)R⊕, per-object albedo U(0.2, 0.75));
the remaining ZTF undershoot (see #3) propagates into the combined total. The published-table
arithmetic (0.564 + 0.050 + 0.171 = 0.785) is preserved exactly via the shared
`p9_core::analysis::surveys` table as a cross-check.
- Pinned: `crates/p9-2024-panstarrs/src/combined_exclusion.rs:227-252` (table identity to 1e-10;
  population-based values within stated tolerances).

### 3. p9-2021-ztf — exclusion fraction: computed 0.477 vs published 0.564

Our per-orbit model (sky position from elements, δ > −30° footprint, logistic
magnitude efficiency, linking factor), fed by the BB21-faithful reference population
(r₉ = (m₉/3)R⊕, per-object albedo U(0.2, 0.75)), under-excludes by ~9 points. The dominant
remaining model difference is the depth semantics: we treat V = 20.5 as a 50% logistic
midpoint, while BB21 define it as the 95%-completeness depth for ≥7 detections — an
effectively deeper survey (issue #255). (Before the population fix the model computed 0.609:
a reference population ~0.6 mag brighter than BB21's stated assumptions more than compensated
for the shallow efficiency reading.) Untuned — no free parameter was adjusted to close the gap.
- Pinned: `crates/p9-2021-ztf/src/exclusion.rs` (asserted within 0.12 of 0.564, fed by the
  seeded BB21-faithful reference population).

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

### 7. p9-2025-iras-akari — expected chance pairs: λ ≈ 22 vs the paper's 22 raw window pairs

The Poisson estimate is conditioned on the paper's Fig. 4 post-cut counts (2,479 IRAS × 108
AKARI-MUSL) over the RETAINED sky (~79% after the |b| ≥ 10° + bulge cuts) with no galactic pair
enhancement — the plane-concentrated density term does not apply to catalogues the cuts have
already de-planed. This lands at λ ≈ 22, matching the paper's 22 raw separation-window pairs
(13 survive their colour-consistency cut). The previous λ ≈ 19-vs-13 discrepancy documented
here was mostly a modeling error (full-sky area × ⟨w²⟩ ≈ 1.59 enhancement), not the count
assumptions. The false-alarm probability for ≥1 surviving pair, P = 1 − e^(−λ), remains the key
new statistic reported alongside.
- Pinned: `crates/p9-2025-iras-akari/src/candidate_search.rs`
  (`lambda_with_paper_counts_matches_raw_window_pairs`).

### 8. p9-2025-perturbation — Hansen asymptotic prefactors corrected against quadrature

Not a paper residual but a deviation from formulas the original code claimed: validating the
asymptotic forms against the convergence-controlled numerical Hansen coefficients forced a
recalibration of X^(−4,1) from 1/8 to 0.71 (the original prefactor was ~5.7× off) and X^(−4,3)
from 1/12 to 1/9.5. If those literals trace to a published table, the source should be re-checked.
- See: `crates/p9-2025-perturbation/src/hansen.rs` (cross-validation tests against
  `p9_core::analysis::hansen::hansen_coefficient`).

### 9. p9-2016-cassini-ranging — favored true anomaly 135.5° (proxy) vs published 117.8°

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

### 10. p9-2021-orbit — full inference pipeline: Reduced scale validates machinery, not the posterior

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

### 11. p9-2019-selfgrav-disk — forced apse aligned (Δϖ = 0), not anti-aligned

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

### 12. p9-2019-ossos-scattering — scattering/detached boundary and the P9 lift are analytic/secular, not N-body

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

### 13. p9-2018-wise-search — at W1 (3.4 µm) a cold P9 is detected by *reflected* light, not thermal emission

Meisner et al. (2018) frame the WISE/NEOWISE shift-and-stack search as a thermal-infrared one.
Our forward model exposes a physics reality that the "thermal-IR" label obscures at W1: at the
band's 3.4 µm and a plausible Planet Nine effective temperature of ≈40 K, hν/kT ≈ 107, so the
blackbody thermal flux is ~10⁻⁵⁵ W/m²/Hz/sr — utterly negligible. The thermal channel only
overtakes reflected sunlight above ≈175 K (at 15 M⊕, 400 AU), far hotter than any Planet Nine.
So the W1 detectability this crate computes is, for a cold body, set by **reflected sunlight**;
the thermal term is retained and would dominate for an implausibly warm planet. (The real search's
W2 4.6 µm band is more thermally favorable; the shared survey table pins the deeper, bluer W1.)
- Consequence for the exclusion: the reflected-light W1 reach is ≈263 AU (5 M⊕) to ≈313 AU
  (18 M⊕), only clipping the near edge of the Brown & Batygin posterior (q ≈ 300 AU), so the
  computed exclusion fraction over that posterior is ≈7% at 6.2 M⊕ and ≈16% at 18 M⊕
  — far below the optical ZTF 56%. This is an honest finding, not a tuned number; it increases
  monotonically with assumed mass as required. (Reach/exclusion re-pinned after fixing the
  geometric-albedo /4 error in `p9_core::analysis::thermal::reflected_flux_jy`, which had made
  all reflected-light magnitudes 1.5 mag too faint; the two core reflected-light paths are now
  cross-pinned against each other on Neptune.)
- Pinned: `crates/p9-2018-wise-search/src/thermal_model.rs` (Planck/reflected model, W1 zero
  point 309.54 Jy), `detectability.rs` (finite max-distance bisection), `exclusion.rs`.

### 14. p9-2025-planet-y — warp amplitude: sub-degree for most of the allowed box, "a few degrees" only at the Earth-mass/close corner

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

### 15. p9-2016-sheppard-etnos — new-object ω clustering is a 2-object corroboration, not a standalone test

Applying the paper's own clustering thresholds (a > 150 AU and q > 35 AU, with the outer Oort
cloud body 2014 FE72 treated separately) to the six newly discovered objects leaves only two in
the clustering sample: 2014 SR349 (ω = 341.3°, aligned ~1° from the published 340° center) and
2013 FT28 (ω = 40.2°, the anti-aligned object). Two angles is too thin for a standalone Rayleigh
test (computed p ≈ 0.25 despite R̄ = 0.87), so the headline ω/ϖ statistics are pinned on the
combined sample (new discoveries + `p9_core::data::etno::BROWN_2017_SAMPLE`), where the new
objects act as corroborating data points exactly as the paper frames them. The combined ω mean is
344.4° (R̄ 0.45, Rayleigh p 0.036) and the combined ϖ mean is 41.8° (R̄ 0.43, Rayleigh p 0.048).
Caveat documented in code: 2014 SR349 and 2013 FT28 also appear in the Brown (2017) table, so the
16-object concatenation double-counts them; `combined_sample_dedup()` (14 objects) confirms the
ϖ clustering and the ~340° ω center survive deduplication.
- Pinned: `crates/p9-2016-sheppard-etnos/src/clustering.rs` (`computed_headline_values_are_pinned`,
  `dedup_combined_still_clusters`, `combined_omega_clusters_near_published_340`).

### 16. p9-2017-resonance-hopping — the analytic n:1 chain overlaps only near P9; hopping at 300–500 AU is an N-body/penetration effect

Becker, Adams et al. (2017) integrate distant TNOs with the nominal Planet Nine (a₉ = 700 AU,
m₉ = 10 M⊕) and find the objects hop between adjacent P9 mean-motion resonances rather than
staying locked in one, because neighbouring resonances overlap. This crate reproduces the
analytic backbone with `p9_core::analysis::resonance`: the n:1 / n:2 resonance locations
(`a = a₉(q/p)^{2/3}`, exact to 1e-9), each resonance's pendulum libration half-width, and the
Chirikov overlap parameter `K = Σδa / Δa` between neighbours.

The honest scale finding: for a **nominal 10 M⊕ P9** the widely-spaced n:1 (and n:2) resonances
out in the 200–500 AU TNO band do **not** overlap (K plateaus at ≈ 0.08 at the 3:1; asserted
isolated in `n_over_1_resonances_isolated_for_nominal_p9`). The neighbour overlap parameter
crosses K = 1 only on the **first-order p:(p−1) chain that crowds toward P9** — the α/(1−α)
Laplace amplitude diverges as a → a₉ — at a computed inner edge `a_hop ≈ 590 AU` (between the
4:3 and 3:2 resonances), the Wisdom (1980) resonance-overlap zone. So the hopping Becker+ see
at smaller a is not adjacent n:1 resonances literally touching; it is the highly eccentric
(e ≈ 0.9) TNOs whose orbits **penetrate** the near-P9 overlap zone at aphelion during their
secular cycles, sampling the chaotic web — which their full N-body integration captures and a
static pendulum-width chain does not. The crate classifies an ETNO as a hopper when its aphelion
Q = a(1+e) reaches a_hop; 2007 TG422 and 2013 RF98 qualify, matching the paper's migrating
class. Sedna also crosses the Q threshold but Becker et al. report it (with 2012 VP113) as
NON-migrating — the aphelion proxy is necessary-not-sufficient and Sedna is its documented
false positive (resonant phase protection, which only the paper's N-body captures, keeps it
locked). No constant is tuned to a target a_hop — the strength `S = α·α/(1−α)·e` is the
bare leading Laplace/eccentric scaling, and K ∝ √μ exactly.
- Pinned: `crates/p9-2017-resonance-hopping/src/chain.rs`
  (`resonance_locations_match_kepler_relation` to 1e-9,
  `overlap_rises_monotonically_and_crosses_hopping_threshold`,
  `below_threshold_isolated_above_overlapping`, `n_over_1_resonances_isolated_for_nominal_p9`,
  `overlap_a_exponent_brackets_validated_chirikov_slope` cross-check vs the validated core
  Chirikov), `src/tno.rs` (`eccentric_distant_tnos_penetrate_overlap_zone`).

### 17. p9-2018-secular-dynamics — critical a ≈ 218 AU vs the paper's a ≳ 250 AU libration threshold

Li, Hadden, Payne & Holman (2018) give a semianalytic secular theory in which a test ETNO's
longitude of perihelion librates (anti-aligned apsidal confinement) above a critical semi-major
axis (their abstract: pericenter libration is maintained for a ≳ 250 AU at q ≈ 30–80 AU). This
crate builds the coplanar conserved secular Hamiltonian as the giants' orbit-averaged J2
quadrupole (`p9_core::forces::j2_secular::combined_j2_jsu`) plus P9's exact numerical Gauss-ring
average (`p9_core::analysis::secular::numerical_secular_hamiltonian` — never re-implemented), then
linearises it in the eccentricity vector to read off the forcing `B` (∝ GM₉·e₉) and the giants'
free precession `A`. The headline results reproduce cleanly:

- Forced eccentricity `e_forced = |B/A|` rises monotonically with α = a/a₉ and **exactly linearly
  with P9's mass** (B ∝ GM₉, A mass-independent) — pinned at a = 100 AU over 2→20 M⊕.
- The apse circulates below and librates (anti-aligned, B < 0) above a computed `a_crit ≈ 218 AU`.

The honest residual: `a_crit ≈ 218 AU` falls ~30 AU short of the paper's ≳ 250 AU. Two documented
reductions account for it: (1) the model is the **coplanar (i = 0)** sector of Li et al.'s theory
(it omits the inclination/Kozai coupling that the full paper also treats), and (2) the free-
precession denominator `A` is the **giants' J2 term alone** — P9's own axisymmetric quadrupole
curvature is a small, quadrature-noisy contribution at the linearisation eccentricity and is
excluded from `A` (it is retained in the full `hamiltonian`). Neither tolerance was tuned; the
nominal P9 (10 M⊕, a₉ = 500 AU, e₉ = 0.6) is the paper's labelled favoured box. The analytic
secular forcing is cross-checked against the exact ring average to ~1% (the hexadecapole O(α²)
truncation residual), and the giants' precession against the closed form (3/2) n J2 (R/a)² (1−e²)⁻²
to ~1e-8.
- Pinned: `crates/p9-2018-secular-dynamics/src/secular_model.rs`
  (`test_critical_semimajor_axis_separates_regimes` — a_crit within 80 AU of 250;
  `test_forced_eccentricity_increases_with_mass`/`_with_alpha`;
  `test_analytic_quadrupole_matches_numerical_ring`; `test_free_precession_matches_closed_form`).

### 18. p9-2021-perihelion-gap — gap statistic reproduced, debiased completeness is not

Oldroyd & Trujillo (2021) report a deficit in the perihelion distribution of distant TNOs at
q ≈ 50–65 AU separating the Neptune-coupled scattering ETNOs (low q) from the detached
inner-Oort-cloud Sednoids (high q). The crate bins the perihelia of the vetted
`p9_core::data::etno::BROWN_2017_SAMPLE` plus a documented `EXTENDED_DISTANT_TNOS` table
(re-transcribed from JPL SBDB 2026-07 after several rows were found to match no real orbit
solution — the old table had miscoded 2015 KH163 *into* the gap at q ≈ 58 vs its real 39.9,
and the genuine gap object 2021 RR205 *out* of it at q ≈ 89 vs its real 55.6) and computes a
dip statistic (gap-window count density ÷ mean flank density). Two honest epochs:

- **Paper epoch** (`paper_epoch_perihelia()`, excludes the post-submission discovery
  2021 RR205): the published two-sided deficit reproduces — dip ratio ≈ 0.12, gap density
  below both flanks, minimum-density bin at q = 52.5 AU inside the published 50–65 AU window.
- **Current knowledge** (`observed_perihelia()`): 2021 RR205 (q = 55.6) sits squarely in the
  window and 2015 TG387's current solution (q = 64.8) falls at its upper edge, so the window
  holds 2 objects — still a deficit relative to the dense low flank (dip ≈ 0.25) but no longer
  emptier than the (thin: Sedna, 2012 VP113) high flank. The gap has partially filled in since
  publication, which matters directly for the sednoid/IOC boundary in the search-space work.

The significance is now computed against a fitted single-continuous null (shifted exponential,
MLE; `synthetic::continuous_null_dip_p_value`): the paper-epoch deficit has exceedance
p ≈ 0.040 (~2σ, matching the paper's moderate-significance framing), weakening to p ≈ 0.11 for
the current sample — the previously missing load-bearing number (the old tests only contrasted
a constructed bimodal draw with a uniform control). A seeded two-population synthetic draw
reproduces the same gap (dip ratio < 0.6) while a single continuous uniform control does not
(dip ratio ≈ 1) — concretizing the paper's "two populations, not one continuous distribution"
argument. The high-q scattering edge is tied to `p9_core::analysis::resonance::critical_perihelion`
(q_crit(a) stays below the gap top for ETNO semi-major axes and reaches 50 AU only at a > 700 AU),
and the detached IOC onset to q ≳ 65 AU.

This is partly observational: the published gap is established after inverting survey detection
biases against faint, distant, high-q objects, which this reduced model does not do. We reproduce
the *gap statistic* (a real deficit at the published q and a two-vs-one-population contrast), not
the full debiased completeness argument, and we do not assert the paper's ∼20% relative-abundance
prediction as a computed number (kept as `published::GAP_RELATIVE_ABUNDANCE` for reference).
- Pinned: `crates/p9-2021-perihelion-gap/src/distribution.rs`
  (`paper_epoch_sample_has_a_perihelion_deficit_at_50_to_65`,
  `current_sample_gap_partially_filled_by_post_paper_discoveries`,
  `observed_gap_center_matches_published_window`), `src/synthetic.rs`
  (`two_population_shows_a_gap`, `single_continuous_population_has_no_gap`), `src/boundaries.rs`
  (`scattering_boundary_sits_below_the_gap`).

### 19. p9-2022-uranus-tilt — single-resonance capture asymptotes to ~90°, not the observed 98°

The Colombo-top reduction (one nodal frequency, Henrard & Murigande 1987 / Ward & Hamilton 2004,
the analytic backbone of Lu & Laughlin 2022) drives the spin axis along the resonant Cassini-2
branch, whose obliquity increases monotonically toward but never past **90°** as α/|g| grows.
Our adiabatic sweep (|g| swept down by ~10² as Planet Nine migrates) cleanly carries Uranus' spin
from θ ≈ 20° to ≈85–89°, independent of the sweep rate (adiabatic invariance). The paper's headline
**98°** (and its 105.6° maximum, 81% within 5%) is reached only in the *full N-body* runs, where the
planet's complete multi-frequency nodal spectrum and libration overshoot push the axis past the
single-resonance 90° asymptote. Reproducing the >90° overshoot requires the full forced-element
spectrum the reduced model omits, so the headline test pins the computed ≈85° within ~18° of 98° and
documents the gap rather than tuning it away. The two first-principles pieces that *do* match the
paper exactly are pinned tightly: the present-day spin-axis precession constant α = 0.0466 arcsec/yr
(paper 0.045, satellite q/l terms from Eq. 3) and the Cassini critical ratio (α/g)_crit = 1.74 at
I = 20° (Eq. 8). The monotone "needs a primordially fast α" requirement is also reproduced: across a
fixed P9-set |g| band the present-day α reaches only ~67°/weak engagement while a 100× enhanced α
drives the spin near 90°.
- Pinned: `crates/p9-2022-uranus-tilt/src/spin_axis.rs` (`alpha_matches_paper_present_value`),
  `crates/p9-2022-uranus-tilt/src/cassini.rs` (`critical_ratio_matches_paper_form`),
  `crates/p9-2022-uranus-tilt/src/resonance_capture.rs` (`adiabatic_sweep_drives_high_obliquity`,
  `present_day_alpha_cannot_tilt_uranus_fully`).

### 20. p9-2016-commensurabilities — angle grouping is significant; the commensurability statistic does NOT beat a random control

de la Fuente Marcos, de la Fuente Marcos & Aarseth (2016) argue the ETNOs are dynamically grouped
in their orbital angles AND near mean-motion commensurabilities. Computing both from the vetted
`p9_core::data::etno` sample against a seeded uniform control splits the claim cleanly:

- **Angle grouping reproduces and is strong.** Using `p9_core::analysis::circular`, the longitude
  of perihelion ϖ has R̄ ≈ 0.45 and a small-n-corrected Rayleigh p ≈ 0.02; the node Ω is the
  weaker grouping (p < 0.2). A seeded uniform-angle Monte Carlo gives an exceedance of a few
  percent, agreeing with the analytic Rayleigh p to within a factor of ~3. This is the robust,
  distribution-free signal.

- **The commensurability statistic does NOT survive a random control.** The observed ETNO period
  ratios sit close to small-integer ratios in absolute terms (mean pairwise distance to nearest
  p/q ≈ 0.04 at max integer 9, planet-scan best ≈ 0.03 at a₉ ≈ 600 AU), but a uniform random
  population drawn over the same semi-major-axis range is *at least as commensurate* ~97% of the
  time (pairwise), and the a₉ scan finds a comparably commensurate planet for random sets too
  (exceedance ≈ 0.53 at max integer 9, worse at low order). This is the Bailey, Brown & Batygin
  (2018) degeneracy made quantitative: the small-integer ratio grid is dense enough that "near a
  commensurability" is not a rare event. The tests assert the honest direction — observed small in
  absolute terms, random control NOT beaten — rather than manufacture significance.

  Resonance-location arithmetic cross-checks to 1e-9 against
  `p9_core::analysis::resonance::resonance_semi_major_axis` (an object placed exactly in an
  interior p:q resonance has period ratio exactly p/q and commensurability distance 0).
- Pinned: `crates/p9-2016-commensurabilities/src/clustering.rs`
  (`test_observed_angles_cluster_significantly`, `test_observed_beats_uniform_control`,
  `test_mc_matches_rayleigh_order_of_magnitude`) and `src/commensurability.rs`
  (`test_resonance_location_arithmetic_cross_checks`,
  `test_pairwise_commensurability_does_not_beat_random_control`,
  `test_planet_commensurability_localizes_but_not_significant`).

### 21. p9-2020-secular-octupole — coplanar octupole libration center is aligned, not anti-aligned; K_OCT calibrated to the ring average

Köhne & Batygin (2020) extend the orbit-averaged ETNO–P9 secular theory to octupole order. This
crate reproduces the MECHANISM on top of `p9_core::analysis::secular`: the coplanar quadrupole
(reused `coplanar_quadrupole`) is a function of `e` alone (apsidal circulation, no Δϖ structure),
and adding the octupole term `+C_oct e e₉ cos Δϖ` opens a Δϖ libration island. The headline test
shows the SAME initial condition circulating under quadrupole-only and librating once the octupole
is on. The dimensionless octupole strength ε_oct = 4·K_OCT·α·e₉/(1−e₉²) is ≈ 1.4 for the prototype
ETNO (a = 250 AU, e₉ = 0.6), confirming the octupole is not a small correction for the relevant
orbits.

Two honest residuals:
- **Libration center is APSIDALLY ALIGNED (Δϖ = 0), not anti-aligned.** The coplanar octupole
  elliptic fixed point is at Δϖ = 0, e* = ε_oct/3 (where −3C_quad e balances the octupole slope).
  The observed ETNOs are ANTI-aligned with Planet Nine; that sign requires the inclined /
  perihelion-detached / resonant dynamics, not the coplanar secular fixed point. Documented as out
  of reach of the coplanar reduction (the same caveat as p9-2019-selfgrav-disk §11).
- **K_OCT = 1.05 is calibrated, not textbook.** The analytic octupole coefficient is anchored to
  the cos(Δϖ) Fourier amplitude of `p9_core`'s EXACT numerical ring average at small α (the
  convergent regime): the fit gives K_OCT ≈ 1.05 in the e→0, e₉→0 limit, with the e₉ dependence
  captured exactly by (1−e₉²)^(−5/2). The textbook pure-octupole multipole coefficient is 15/16 =
  0.9375; the ~12% excess is the hexadecapole-and-higher contribution to the same harmonic that the
  exact ring average folds in. Analytic vs ring-average agree to <6% at α ≈ 0.03 (the residual is
  the O(e²) correction to the leading-in-e amplitude).
- Pinned: `crates/p9-2020-secular-octupole/src/octupole.rs`
  (`analytic_octupole_crosschecks_ring_average_at_small_alpha`, strength/scaling tests),
  `src/libration.rs` (`octupole_produces_libration_island`, `quadrupole_only_circulates`,
  `libration_center_is_aligned_apse`), `src/precession.rs` (quad-vs-octupole precession).

### 22. p9-2016-secular-resonance — circulation is reached by detuning, not by raising e (the resonance is broad)

Beust (2016) explains the ETNO apsidal clustering as capture in a secular resonance with Planet
Nine: where the test particle's free apsidal precession matches P9's, the relative longitude of
perihelion Δϖ = ϖ − ϖ₉ librates rather than circulates, in an **apsidally anti-aligned** (Δϖ = π),
**high-eccentricity** island. The crate reproduces this from the `p9_core::analysis::secular`
Gauss-ring average: it extracts the apsidal Fourier harmonics of H(e, Δϖ) (a coherent DFT over
Δϖ — far better conditioned than finite-differencing the angle, since the apsidal forcing c₁ is a
~10⁻⁴ fraction of the axisymmetric bulk), builds the second-fundamental-model-of-resonance
pendulum K(Δϖ, Γ) = a₀(Γ) − g₉·Γ + c₁*·cos Δϖ, integrates it (RK4, K conserved to ~10⁻⁶), and
locates the anti-aligned libration island. The computed center is Δϖ = π and the resonant e sits
in the high-e regime, matching the abstract; the resonance/island strengthens linearly with P9's
GM (c₁ ∝ GM, pinned to a 4× ratio at fixed e). Every object in the vetted Brown (2017) sample sits
at a semi-major axis that admits an anti-aligned island under the nominal P9.

Honest finding: at the nominal P9 (10 M⊕, 700 AU, e = 0.6) the resonance is so **broad** that its
separatrix spans essentially the entire physical eccentricity range at a fixed ETNO semi-major
axis (the action half-width δΓ_sep corresponds to Δe ≈ 0.3 at a = 400 AU). There is therefore no
circulating regime reachable by changing the particle's eccentricity alone — the circulating
(non-clustered) case is the *non-commensurate* one, reached by detuning g₉ off the free-precession
band. The crate demonstrates both routes to circulation that do exist: detuning g₉ (the physical
non-resonant particle, `resonant_librates_detuned_circulates`) and a large-amplitude launch across
the angle separatrix at the narrower a = 250 AU resonance (`small_amplitude_librates_large_circulates`).
Beust's published libration *amplitude/period* numbers are configuration-specific phase-portrait
read-offs not quoted as single scalars in the abstract, so the crate pins the qualitative
structure (anti-aligned center, high-e island, mass scaling, libration vs circulation) and the
self-consistency of the pendulum rather than a published amplitude figure.

Honest residual (issue #226): the high-e placement of the tuned portrait island is an INPUT,
not a derived result. Solving the physical commensurability — g₉ from P9's own apsidal
precession under the giants' quadrupole (≈2×10⁻¹³ rad/day, a ~10¹¹-yr apse period) with the
giants' J2 included in the particle's free precession — puts the resonance at e* ≈ 0.24 for the
nominal P9 at a = 250 AU (`portrait::physical_island`, pinned): the free precession is
retrograde at low e (P9 ring dominating), crosses the tiny g₉ once on its ascending branch, and
never returns to it at high e. Reproducing Beust's high-eccentricity anti-aligned islands as a
*derived* feature requires his full frame treatment (mutual precession), beyond this reduced
pendulum; the tuned portrait is retained as the paper's-regime demonstration with that caveat
documented at the definition.
- Pinned: `crates/p9-2016-secular-resonance/src/portrait.rs`
  (`island_exists_for_nominal_p9_and_is_anti_aligned`, `circular_perturber_has_no_island`,
  `resonance_strengthens_with_p9_mass`, `every_etno_admits_an_anti_aligned_island`),
  `src/libration.rs` (`resonant_librates_detuned_circulates`,
  `small_amplitude_librates_large_circulates`, `libration_center_is_anti_aligned`,
  `amplitude_grows_with_displacement`), `tests/cross_check.rs`
  (`hamiltonian_matches_p9_core_ring_average` to 1e-12, `island_robust_to_quadrature_refinement`).

### 23. p9-2018-bp519 — secular inclination pumping reaches BP519's i from a moderately-excited start, not from genuinely low i

Becker et al. (2018) argue that 2015 BP519 (i ≈ 54°, e ≈ 0.92, a ≈ 449 AU) is naturally produced
by an inclined Planet Nine that secularly pumps the inclination of scattered objects. This crate
reproduces the pumping with a single-perturber doubly-averaged secular model: P9 as a numerical
Gauss ring (`p9_core::analysis::secular::numerical_secular_hamiltonian`), the test particle's
(e, i, ω, Ω) advanced by RK4 under Hamilton's equations in Delaunay variables (∂H by finite
difference). An inclined P9 drives a secular (Kozai-Lidov/octupole) cycle that pumps a high-e detached
particle's inclination from a mutual i₀ = 40° **up to ≈ 56°** — squarely in BP519's band — while
the orbit stays prograde and detached (e < 1); without P9 (giant-planet J2 ring only, which has no
inclination torque) the inclination is exactly flat. At this reduced scale (N_QUAD = 12, ~70 Myr)
the inclination is lifted largely through the node/octupole channel with the eccentricity nearly
unchanged over the recorded window, rather than the textbook quadrupole e–i anti-correlation.

The honest residual: the **single-perturber secular average cannot pump from a *genuinely* low
inclination** (i ≲ 20°). Kozai libration only raises i for mutual inclinations near/above the
critical angle (≈ 39°, where Θ = √(1−e²)cos i is small enough), and within a cycle the inclination
is bounded — starting below critical, i only oscillates downward. Reproduced exhaustively in
calibration: low-i (5–20°) starts stay low at every (a, e, ω, P9-mass) tried. So the demonstration
starts the particle at a *moderately excited* mutual inclination (i₀ = 40°, representing an
ecliptic-plane scattered object seen in the frame of a P9 tilted by i₉ ≈ 40°, plus prior
excitation), from which P9 carries it to BP519-like values. Bridging genuinely-low to 54° is the
Gyr-scale resonant + Neptune-scattering N-body process of Batygin & Brown (2016), beyond a clean
single-perturber secular average. A second modelling consequence: because the Gauss-ring
Hamiltonian uses P9's plane as the reference, the particle's `i` argument *is* the mutual
inclination, so varying P9's own `i` field has no effect — the sin(i₉) dependence is encoded as
the initial mutual inclination i₀ ≈ i₉, and the test varies that.
- Pinned: `crates/p9-2018-bp519/src/pumping.rs` (`p9_pumps_inclination_to_bp519_like` ≥ 50°,
  `no_p9_control_leaves_inclination_flat` range < 2°,
  `moderate_mutual_inclinations_all_pump_into_bp519_band`,
  `pumping_rate_grows_with_p9_mass`, `pumping_is_kozai_e_i_exchange`,
  `hamiltonian_conserved_along_flow` |ΔH/H| < 1e-7); BP519 elements pinned in
  `src/bp519.rs`; ϖ-cluster context in `src/clustering.rs`. Note: after fixing a sign error in
  the dω/dt Delaunay chain rule (the flow now conserves H to rounding), max inclination is NOT
  monotone in the starting mutual inclination — both 25° and 40° starts pump into the ≥55° band
  (to ~71°/~64°); the previously pinned monotone growth was an artifact of the broken flow.

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
| S&T 2016 combined ω cluster center | 344.4° | ~340° | `p9-2016-sheppard-etnos/src/clustering.rs` (within published spread) |
| S&T 2016 combined ϖ Rayleigh p | 0.048 | significant | `p9-2016-sheppard-etnos/src/clustering.rs` (< 0.1) |

---

## Tuned/assumed parameters (documented, not derived)

These are free parameters chosen to match an observable, or assumptions standing in for data we
do not ingest. Each is commented at the definition site; listed here so nobody mistakes them for
derived quantities:

- `p9-2022-des` `night_usability = 0.255` (`survey_model.rs`) — the single calibration knob
  scaling per-night completeness so end-to-end recovery matches the paper; stands in for
  weather/chip-gap/masking losses we don't model per-exposure.
- `p9-2017-bias` galactic-plane suppression (center crossing σ=40°/floor 0.02, anticenter
  σ=15°/floor 0.5) and north/south asymmetry — tuned so the combined bias-corrected p lands near
  the published 2.5e-4; the real quantity requires actual survey pointing histories.
- `p9-2019-clustering` / `p9-2021-orbit` survey-bias null: longitude suppression near
  λ ≈ 95°/275° + δ > −30° coverage — same caveat.
- `p9-2025-clustering` selection weighting w(ϖ) = 1 + 0.3cos(ϖ − 60°) in the bias-resampling MC.
- `p9-2017-ossos-bias` discovery-longitude windows: now the PUBLISHED OSSOS block centres
  (Bannister et al. 2018 Table 1; two clumps, λ ≈ 7–48° and 203–240°). With the real bimodal
  layout the bias-induced R̄ from isotropy is ≈ 0.29 (pinned) — weaker than the ~0.5–0.65 a
  previously invented single-side window layout manufactured. The window half-widths remain
  simplified block-clump extents (no per-field depth/cadence history).
- `p9-2021-napier-critique` composite selection function (A₁ = 0.90, ϕ₁ = 52°, A₂ = 0.09,
  ϕ₂ = 52°) — an assumed stand-in for the OSSOS/DES/S&T pointing histories, with the lobe
  placed on the observed cluster direction and an order-ten contrast chosen so the debiasing
  flip lands in the published band. The crate's sensitivity tests document the contingency:
  the R̄-only consistency statistic is phase-blind (amplitude alone decides the flip), and the
  direction-aware statistic shows the flip failing outright with the lobe rotated 90° off the
  cluster. Deriving the selection from the real footprints is the survey-geometry
  consolidation work.
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

### 23b. p9-2024-siraj-orbit — MAP 7.0 M⊕ from the vetted 10-object sample vs the published 4.4 M⊕ (stable-21)

The ridge normalization is a one-point calibration frozen to the paper's own headline
(κ_REF ≈ 1.09 from 2.7σ at n = 21 → the published 4.4 M⊕/290 AU³); the input sample's von Mises
κ̂ then sets S_obs relative to that anchor. The vetted `p9_core::data::etno` 10-object ϖ sample
is more concentrated than Siraj, Chyba & Tremaine's stable-21 (κ̂ ≈ 1.7), so the inferred MAP is
κ̂/κ_REF ≈ 1.6× heavier at the same a_p — a sample difference, not a modeling residual. A
previous version calibrated from the live sample's κ̂, which cancelled the data identically and
returned the published orbit for any input; falsifiability tests (uniform input ⇒ collapsed
mass, tighter input ⇒ heavier MAP) now pin the corrected behavior.
- Pinned: `crates/p9-2024-siraj-orbit/src/lib.rs` (`test_map_tracks_sample_concentration`),
  `src/posterior.rs` (`test_inference_is_falsifiable`).

### 23c. p9-2016-inclination-instability — C_GROWTH is an ansatz, not a derived eigenvalue

The growth-rate prefactor (γ = (M_d/M_⊙)·n / C_GROWTH with C_GROWTH = 2π, i.e. one e-fold per
disc secular time) is anchored to Madigan & McCourt's N-body measurement that the inclination
e-folds in about one t_sec — it is NOT derived here. Deriving it requires the full N-ring
linearised secular operator (the collective eccentric-disc mode), beyond this reduced crate.
Tests pin the paper's order-unity τ/t_sec band and the M_d/a scalings against independent
p9-core machinery rather than the constant against itself; the suppression-criterion crossover
(giants' precession vs disc frequency) is the crate's genuinely computed content.
- Pinned: `crates/p9-2016-inclination-instability/src/growth.rs`
  (`efolding_time_in_the_published_order_unity_band`), `src/suppression.rs`.

### 23d. p9-2026-apsidal-clustering — tuned synthetics demonstrate the mechanism; the real sample lands at σ ≈ 2.6

The stable-21/25 synthetic samples are TUNED (their concentrations were solved from the
published 2.7σ → 1.9σ outputs), so those headline values are reproduced by construction — they
pin only the dilution mechanism. The genuine data test runs the same von Mises
log-likelihood-ratio estimator on the vetted 10-object ϖ sample (`p9_core::data::etno`):
κ ≈ 1.73, σ ≈ 2.64 — independently at the paper's ~2.7σ scale (pinned). The asymptotic Wilks
p = exp(−Λ) is cross-checked against the small-n-corrected Rayleigh series from p9-core on the
same sample (agreement within ~0.3σ; the Wilks form is the more optimistic). The paper's exact
stable-21/25 element tables are not transcribed.
- Pinned: `crates/p9-2026-apsidal-clustering/src/samples.rs`
  (`real_etno_sample_is_significantly_clustered`).

### 24. p9-2025-stacking — analytic matched-filter / metric reproduction, not a pixel pipeline

Geringer-Sameth et al. (2025) build a full image-stacking search with an orbit-
parameter Fisher metric. This crate reproduces the *analytic backbone*: the √N
matched-filter SNR and 1.25·log10(N) depth gain, the leading-order rate-error
metric g_vv = T²/(24θ²) (cross-checked against the exact Gaussian-overlap stack
to ~1%), the resulting trial-orbit count (~10^8–10^9 for a ZTF six-year
baseline), and the logarithmic look-elsewhere penalty (~6.5σ / ~0.3 mag for 10^9
trials). Reductions, documented and pinned: the metric is the linear on-sky
(rate-only) sub-block rather than the full 5–6-element bound-orbit metric, and
the ~27th-mag reach is the √N law extrapolated to ~1.6e5 frames rather than a
pixel-level coadd. The headline scalings and the stacking-gain-≫-trials-penalty
conclusion reproduce; absolute trial counts depend on the adopted rate range /
PSF / SNR tolerance, kept as labelled inputs.
- Pinned: crates/p9-2025-stacking/src/{matched_filter,orbit_metric,significance}.rs

### 25. p9-2022-iras-candidate — chance-association count and the back-derived candidate flux

Rowan-Robinson (2022) publishes the *fitted* candidate distance (225 ± 15 AU)
and mass (3–5 M⊕) but no explicit candidate flux in Jy, and does not tabulate
the subset sizes used in the cross-match. Two honest residuals follow:

1. **Back-derived flux.** The reference 60/100 µm fluxes (≈0.35 / 0.86 Jy) are
   the values *self-consistent* with the published 225 AU / 4 M⊕ under this
   crate's blackbody + inverse-square model at an assumed T_eff = 40 K (a cold
   internally-heated ice giant; solar equilibrium at 225 AU is only ~17 K). They
   are labelled as back-derived, not quoted from the paper. The test pins the
   round-trip (flux → 225 AU → flux) and that the implied distance lands within
   the published ±15 AU, not a flux the paper never stated.
2. **Chance-association count.** The all-FSC-against-all-FSC Poisson estimate
   runs to ~10⁵–10⁶ (a loose upper bound, pinned). Restricting to the
   ~3000 unidentified 60 µm primaries cross-matched against the sparse
   single-HCON Reject File (~5×10⁴ sources) over the 2–35 arcmin window gives
   ~4×10³ raw geometric coincidences — the right order for the paper's "several
   hundred candidate associations", which are those *after* Scanpi vetting. The
   subset sizes are order-of-magnitude reference values, so the test pins the
   order of magnitude (hundreds to a few thousand), not a precise count.
- Pinned: crates/p9-2022-iras-candidate/src/{distance.rs,chance.rs}; reference
  constants in src/lib.rs (`REF_CANDIDATE`).

### 26. p9-2023-lsst-strategy — for the *nominal-orbit* P9 population the binding LSST constraint is the footprint, not depth or linking

Schwamb et al. (2023) argue that LSST's cadence/footprint tuning drives the
discoverable fraction of slow distant movers. This crate runs a seeded
reference Planet Nine population (the `p9-2021-orbit` posterior emulator +
`p9_core::analysis::photometry`) through a real LSST detection model built on
the reused `p9-survey` Rubin/LSST footprint preset (δ ∈ [−75°, +12°], |b| > 10°,
coverage 0.85) and the shared single-visit (r ≈ 24.5) / ten-year-stack
(r ≈ 27.0) depths kept as labelled `published::` constants. The per-orbit
discovery probability is `footprint(δ, b) · coverage · P(link | V)`, with the
linking probability the binomial survival of ≥ N detections out of the field's
visits at the single-visit recovery efficiency.

All four published directions reproduce as *monotone* effects: the discoverable
fraction is in (0, 1); requiring more visits-for-linking never raises it;
shrinking the declination extent never raises it; the ten-year co-add never
lowers it. The honest finding behind the magnitudes: the nominal-orbit P9
population is **bright** — V p50 ≈ 20.4, p99 ≈ 24.1 (dwell-weighted current
positions, BB21-faithful photometry), essentially all brighter than the 24.5
single-visit depth — so per-visit
ε ≈ 1 and a single LSST visit already reaches essentially every nominal
solution. The **footprint is therefore the binding constraint** (the dec/|b|
gate moves the fraction by tens of percent), while the depth and linking levers
move it only at the sub-percent level on the full population because almost
nothing sits near the per-visit limit. The depth/linking sensitivity Schwamb
et al. emphasize is genuine but localised to the **faint distant tail**: at a
shallow per-visit cadence (depth ≈ 21, e.g. dark time split across many
filters) both levers bite hard (> 0.05 swings), pinned as
`linking_lever_bites_hard_near_the_depth_limit` and
`stack_depth_increases_discovery_near_the_limit`. No constant is tuned to a
target fraction — there is no single published discoverable-fraction scalar to
hit; the model reproduces the *sensitivities*, not a number.
- Pinned: `crates/p9-2023-lsst-strategy/src/strategy.rs`
  (`baseline_matches_published_reference_constants`,
  `single_visit_depth_matches_p9_survey_rubin_preset`,
  `binomial_survival_matches_known_values`, linking/depth monotonicity),
  `src/lib.rs` headline tests (probability-in-(0,1), the four monotone levers,
  the two faint-tail sensitivity tests, seed stability).
