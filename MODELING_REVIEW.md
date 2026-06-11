# Modeling Review — Planet Nine Workspace

Full-repo review of the modeling (physics fidelity, numerical methods, statistics) across all
23 crates, conducted June 2026. Findings are grouped by theme, then by crate. Each finding has
file:line references and a High/Medium/Low impact tag.

---

## Executive summary

The workspace has a solid, reproducible scaffold (seeded RNGs in most crates, consistent
workspace structure, a real WHM/Bulirsch-Stoer integrator stack in p9-core). But the review
found five repo-wide problems that matter more than any single bug:

1. **Verified physics errors** in core machinery — a sign error in the hyperbolic orbit
   conversion, a resonant-angle definition with p and q swapped, a √2 inconsistency in the
   Chirikov chaotic-boundary formulas, a wrong distance law for reflected sunlight, and a
   dimensionally invalid timescale formula. These produce actively wrong output, not just
   simplified output.
2. **Claimed-but-dead physics.** The giant-planet J2 quadrupole (`p9-core/src/forces/j2_secular.rs`),
   the galactic tide, the hybrid close-encounter integrator, and the democratic-heliocentric
   transforms all exist in p9-core and are **never called by any crate**. Several paper crates
   document these effects as included; the simulations actually run bare Sun + P9.
3. **Hard-wired paper results masquerading as computations.** Many crates carry a
   `paper_values()` struct holding the published numbers, while the compute path is a cartoon
   that is never validated against them (and in several cases *cannot* reproduce them).
   Tests then assert the hard-coded constants equal themselves.
4. **Observational data errors.** The 2024-neptune-crossing TNO table contains objects that are
   not TNOs (2016 QV89 is a ~30 m near-Earth asteroid); the 2018-resonance KBO table has
   scrambled semi-major axes that contradict the vetted table in p9-2017-bias.
5. **Statistics gaps.** The survey-bias-aware null hypothesis — the entire point of Brown 2017
   and Brown & Batygin 2019 — is unimplemented (nulls are uniform), and two crates use
   unseeded RNG in Monte Carlo paths.

### Top 10 fixes by impact

| # | Fix | Where |
|---|-----|-------|
| 1 | Sign error in hyperbolic `elements_to_cartesian` (sin ν negated; verified by round trip) | `p9-core/src/types.rs:291` |
| 2 | Resonant angle has p↔q swapped — libration detection can never fire | `p9-2016-evidence/src/resonance.rs:77` |
| 3 | Wire the J2-effective giant-planet quadrupole into the WHM kick (4+ crates run without it) | `p9-core/src/integrator/whm.rs:37-44`, `forces/j2_secular.rs` |
| 4 | Replace fabricated TNO table (2016 QV89, 2014 LU28 are not TNOs) with SBDB-checked sample | `p9-2024-neptune-crossing/src/observed_tnos.rs:43-237` |
| 5 | Reflected-light magnitude law: use m ≈ H + 5 log₁₀(rΔ), not the stellar 5 log₁₀(d/10) | `p9-2022-des/src/recovery_analysis.rs:67-68` |
| 6 | Chirikov exponent inconsistency: `exp(−q²/4a_N²)` vs `q_crit` derived from `exp(−q²/2a_N²)` — two public APIs disagree on the chaos boundary | `p9-2021-stability/src/chirikov.rs:21,45-52` |
| 7 | Hybrid integrator double-counts planet forces and freezes planets during BS steps; then nothing uses it anyway | `p9-core/src/integrator/hybrid.rs:60-89`, `bulirsch_stoer.rs:22-39` |
| 8 | Thermal model radius (1.86 R⊕ at 10 M⊕ vs documented ~3.5 R⊕) underestimates IR flux ~3.5×, flipping the IRAS detectability conclusion | `p9-2025-iras-akari/src/thermal_model.rs:41-43` |
| 9 | Proper-motion model uses total speed as transverse rate and omits geocentric parallax — cannot reproduce the paper's 500–700 AU window | `p9-2025-iras-akari/src/proper_motion.rs:42-46` |
| 10 | Bias-aware clustering nulls: make `perihelion_bias` depend on ϖ and feed it into the 2019 clustering null | `p9-2017-bias/src/bias_function.rs:30-61`, `p9-2019-clustering/src/clustering_analysis.rs:67-100` |

---

## 1. p9-core (foundation)

### High

- **H1. Hyperbolic branch sign error in `elements_to_cartesian`** — `types.rs:291`.
  `sin_nu = -(e²-1).sqrt()·sinh(H)/(e·cosh(H)-1)`; the leading minus is wrong
  (correct: `sin ν = +√(e²−1)·sinh H/(e·cosh H−1)`). Hyperbolic states are mirrored about
  perihelion (M → −M; inbound/outbound swapped). Verified numerically: round trip converts
  M = 0.7 to M = −0.7, and r·v has the wrong sign. Corrupts any escaping-orbit work
  (Oort, flybys). The elliptic branch (`types.rs:286`) is correct.
- **H2. Hybrid integrator double-counts planetary forces** — `integrator/hybrid.rs:60-85` +
  `bulirsch_stoer.rs:22-39`. BS-flagged particles still receive the dt/2 planet kicks, then
  `bs_step` integrates Sun + direct planet force over the full dt — the perturbation is applied
  ~twice, and the BS force omits the indirect (heliocentric-frame) term the kicks include.
  Chambers (1999) splits with a smooth changeover K(r): kick gets (1−K)·F, BS gets Sun + K·F.
- **H3. Planets frozen during BS encounter steps** — `bulirsch_stoer.rs:22-39`. Static planet
  positions over dt = 300 d (Jupiter moves ~25°). Kepler-drift each planet to the substep epoch.
- **H4. "WHM" scheme is not symplectic; democratic-heliocentric machinery is dead code** —
  `integrator/whm.rs:47-74`, `coords/democratic_helio.rs` (zero callers). Heliocentric positions
  + heliocentric velocities are non-canonical, so the splitting conserves no shadow Hamiltonian;
  secular energy drift at O((m/M)²)/step over 4 Gyr runs. Complete the DLL98 scheme: barycentric
  velocities, direct-only kick, plus the missing Σp²/2M_sun solar-drift step (transforms already
  written, never invoked).
- **H5. Quadrupole secular Hamiltonian has unphysical Δϖ dependence; octupole claimed but
  absent** — `analysis/secular.rs:29-58`. The doubly-averaged quadrupole is axisymmetric in the
  apsidal angle; the `−5e²cos²Δϖ` term is wrong, and the anti-aligned libration islands only
  appear at octupole order (∝ e_p α³), which is not implemented. Worse, the multipole expansion
  diverges for the relevant α ≈ 0.6–0.9; B&B 2016 used numerical Gauss-ring averaging for
  exactly this reason. Fix: numerical double averaging of the exact 1/Δ for α ≳ 0.3.
- **H6. Force modules never wired into any integrator** — `whm.rs:37-44`, `forces/j2_secular.rs`,
  `forces/galactic_tide.rs`. Zero workspace callers of `combined_j2_jsu`,
  `secular_j2_acceleration`, `galactic_tide_acceleration`. Every paper crate using
  `WhmIntegrator` integrates without the solar-quadrupole JSU approximation. Add an
  extra-acceleration hook to the kick step.

### Medium

- **M1. Galactic tide applied along ecliptic z** — `forces/galactic_tide.rs:22-37`. Galactic
  plane is inclined ~60° to the ecliptic; rotate to the galactic frame and back. Radial tide
  terms (G₁, G₂) also missing.
- **M2. J2-effective JSU omits the monopole GM boost** — `forces/j2_secular.rs:27-45`. Absorbing
  J+S+U requires adding ΣGM_JSU ≈ 1.28e-3·GM_sun to the central mass; mean motions/secular
  frequencies otherwise off by ~6e-4. Also J2+J4 ring expansion converges slowly at q ≈ 30 AU
  vs a_U = 19.2 AU; document the validity limit.
- **M3. dt = 300 d under-resolves Jupiter** when the four giants are direct
  (`types.rs:258-269` + `initial_conditions/planets.rs:114-121`): 14.4 steps/orbit vs ≳20 needed.
  Enforce dt ≲ 150–200 d with Jupiter in `bodies`, or default to J2-averaged JSU + direct Neptune.
- **M4. BS step never reduces dt on convergence failure** — `bulirsch_stoer.rs:128-130` silently
  returns "best estimate" at full step exactly where it matters (deep encounters). Halve and
  recurse; at minimum return a convergence flag.
- **M5. Unguarded Newton Kepler solvers** — `types.rs:451-477`, `integrator/kepler_step.rs:54-74`.
  Fixed 50 iterations, no convergence check; hyperbolic starter ignores e (use Danby's
  `ln(2M/e + 1.8)`); no bracketing fallback. Five Kepler-solver copies exist across the
  workspace (see §8) — consolidate one robust solver here.
- **M6. Energy diagnostic isn't the conserved quantity** — `whm.rs:97-147`. Heliocentric
  velocities omit the ½|Σmᵢvᵢ|²/M_sun barycentric correction; test-particle energies pollute the
  planetary budget; total angular momentum never computed anywhere. Compute barycentric E and L.
- **M7. Hybrid encounter detection samples only the step start** — `hybrid.rs:57-58`. Particles
  can cross the 3 R_Hill sphere entirely between checks; test predicted minimum separation over
  [t, t+dt].
- **M8. Scattered-disk rejection sampling silently biases the (a, q) wedge and breaks fixed N** —
  `initial_conditions/scattered_disk.rs:99-104,150-160`. Loop-resample until valid; document the
  resulting distribution.

### Low

- Planet Nine radius formula off by ~4 orders of magnitude (`types.rs:221-222` → ~0.6 km at 10 M⊕).
- Planetary J2/J4 fields stored but never applied (`types.rs:103-105`; `kick.rs:124` uncalled);
  GR absent (fine at q ≳ 30 AU, but state it).
- Circular/equatorial edge cases: prograde assumption `v.x > 0` breaks retrograde equatorial
  orbits (`types.rs:398-405`); use `h.z.signum()`.
- Stumpff functions overflow for large hyperbolic steps (`kepler_step.rs:102-117`); use 4-fold
  argument reduction.
- Constants verified sound against DE440 (GM ratios, PC_AU, ρ_MW, ring J2/J4 coefficients).
  Note `GM_EARTH_MOON` (system) vs `EARTH_MASS_SOLAR` (Earth-only) deserve a comment.
- Invariant test coverage is thin: no element round-trip test (a hyperbolic one would have caught
  H1), energy test too short/loose to catch H4, no L conservation test, no dt-scaling
  (order-verification) test, `tests/` empty.

---

## 2. 2016 crates (evidence, constraints, obliquity, inclined-tnos)

### High

- **H1. Resonant angle p↔q swapped** — `p9-2016-evidence/src/resonance.rs:77`. With
  `a_res = a_p(q/p)^{2/3}`, the stationary angle is φ = q·λ − p·λ₉ + (p−q)ϖ; the code computes
  p·λ − q·λ₉ − (p−q)ϖ, which circulates at n₉(p²−q²)/q even for a perfectly resonant particle.
  `is_librating` will classify every real resonance as circulating — the paper's proposed
  confinement mechanism is undetectable as coded. An integrated 2:1 test orbit would catch it.
- **H2. Multipole-expanded secular Hamiltonian invalid for orbit-crossing geometry** —
  `p9-2016-evidence/src/octupole.rs:39-70,128-152` (+ p9-core H5). Test particles cross P9's
  orbit (α up to 0.79); the expansion doesn't converge. Octupole coefficients
  (A = −35/8 e² + 25/8 etc., cos 3Δϖ harmonic) match no standard form — verify or remove the
  "Mardling (2013)" attribution. The 3D treatment is a bolted-on cos(i). Phase-portrait topology
  will not match the paper. Fix: numerical double averaging (2D Gauss-Legendre over both M's).
- **H3. The headline statistic (joint ϖ + pole clustering, P ≈ 0.007%) is never computed** —
  `p9-2016-evidence/src/kbo_elements.rs:176-195`. Only arithmetic means of angles exist (and
  arithmetic, not circular, means at that). Add circular means + a seeded MC that draws 6 uniform
  (ϖ, pole) pairs and validates against 0.007%.
- **H4. Scattered-disk sims omit the giant-planet quadrupole entirely** —
  `p9-2016-evidence/src/scattered_disk_sim.rs:59` (planar nominal is bare Sun + P9). B&B 2016
  absorbed J–U (planar: +N) into an enhanced solar J2 that sets KBO apsidal precession.
  Precession rates, libration centers, survival fractions all quantitatively off. p9-core's
  `combined_j2_jsu` is dead code. Affects p9-2016-constraints and `phase_portrait.rs:174-176` too.
- **H5. Solar spin-down diverges at early times** — `p9-2016-obliquity/src/solar_model.rs:43-49`.
  `ω ∝ √(t_age/t)` uncapped: at the first step ω ≈ 300× present (beyond breakup);
  `P_SUN_INITIAL_DAYS = 10` only applies at exactly t ≤ 0. Since spin-orbit coupling is largest
  early, the obliquity excitation and the whole required-i₉ survey are inflated.
  Fix: `ω = min(ω(P=10 d), ω_present√(t_age/t))`.
- **H6. Giant planets collapsed to one L-conserving ring** — `solar_model.rs:97-115`,
  `secular_hamiltonian.rs:258,277`. C_sun_gp ∝ Σmᵢ/aᵢ³ is −42% and C_gp_9 ∝ Σmᵢaᵢ² is −45% vs
  the per-planet sums Bailey+ 2016 used — a near-2× systematic in the coupling that shifts the
  required-i₉ contours. Sum couplings per planet.
- **H7. Inclined-TNO sim never uses the hybrid integrator** —
  `p9-2016-inclined-tnos/src/simulation.rs:139-156`. Docstring claims WH/BS hybrid; code
  instantiates plain `WhmIntegrator` (changeover/epsilon config dead). Fixed-step WHM through
  Neptune encounters — the exact mechanism producing Drac/Niku analogs — gives O(1) errors per
  encounter. Also dt = 300 d with Jupiter direct (see core M3).

### Medium

- **M1. p9-2016-constraints survey is scaffolding** — `parameter_grid.rs:107-119`: `GridResult`
  is never filled by any runner; the crate cannot produce the paper's allowed-region figure.
- **M2. Acceptance criteria invented** — `parameter_grid.rs:121-132`, `clustering_metric.rs:146`,
  `inclination_survey.rs:73-86`. Thresholds (n_survivors ≥ 7, R̄ > 0.5, pole cuts) are not
  traceable to the paper; they determine the crate's entire conclusion. Document or recast.
- **M3. Rayleigh p-value approximation poor at n = 6** — `clustering_metric.rs:88` uses
  p ≈ exp(−nR̄²); use the standard small-n series. `confinement_probability` also never compares
  against the observed 6-KBO R̄ (~0.96).
- **M4. Fabricated KBO σ's; clone-stability screening unimplemented** — `kbo_elements.rs:21-144`.
  "6 dynamically stable KBOs" is asserted in data, not reproduced.
- **M5. "Coplanar" phase portraits use an inclined P9 (i₉ = 30°, Ω₉ = 100°) with i = 0
  particles** — `phase_portrait.rs:38-60,168-180`; Δϖ is a broken angle at 30° mutual
  inclination. Also all particles start at perihelion (phase-correlated) and no J2 term.
- **M6. Duplication**: snapshot/run loop copy-pasted 3× (with *inconsistent* snapshot triggers);
  Δϖ wrapping 5×; R̄ computation 3×; two different "paper nominal" P9 configs in the same repo
  (`InclinedTnoConfig::nominal` vs `P9Params::inclined_tnos_2016`). → p9-core circular-statistics
  module + one simulation driver.
- **M7. Kozai-Lidov detection is a static cartoon** — `kozai_lidov.rs:38-67`. Single-epoch
  (a, i) cut; the conservation test (std < 0.2 on √(1−e²)cos i) rejects exactly the
  non-conserving trajectories that produce orbit flips. Classify by e–i anti-correlation in time.

### Low

- Mass-radius R = 3.0 R⊕ M^0.27 overestimates P9 size ~40% (bigger than Neptune at 10 M⊕) —
  `p9-2016-constraints/src/detection_limits.rs:43-46`; skews detectability optimistic.
- Octupole symmetry test asserts |ΔH| > 1e-20 (vacuous) — `octupole.rs:191`.
- High-q threshold doc/code mismatch (100 vs 60 AU) — `parameter_grid.rs:113` vs
  `clustering_metric.rs:60-63`.
- No test pins precession direction/rate against the analytic two-vector frequency
  (obliquity crate); RK4 stages renormalized mid-step (formally breaks 4th order, no
  convergence test).
- `clustering_statistics` mixes linear and circular statistics on bimodal Δϖ.
- Seeds are threaded well but never recorded in serialized outputs.

---

## 3. 2017–2018 crates (bias, dynamics, kuiper-belt, resonance)

### High

- **H1. p9-2018-resonance: J2 quadrupole claimed, never applied** — `simulation.rs:3-5,123-125`.
  "Sun with J2 (representing giant planets)" comment, but `bodies = vec![p9]` and WHM never calls
  any J2 force. Without known-planet precession, the resonance-occupation census is from the
  wrong dynamical problem.
- **H2. p9-2017-dynamics: P9–KBO interaction term has the wrong form** —
  `hamiltonian.rs:104-115`. A Δϖ-dependent term at quadrupole scaling, ∝ e², independent of e₉ —
  predicts apsidal confinement even for circular P9, which is impossible; apsidal coupling first
  appears at octupole order ∝ e·e₉/(1−e₉²)^{5/2}·cos Δϖ. The quadrupole term itself
  (`(1+1.5e²)(1+1.5(cos²i−1))`) doesn't match the standard Legendre form. Add a test that the
  Δϖ term vanishes as e₉ → 0.
- **H3. p9-2017-dynamics: the paper's central experiment is absent** — `simulation.rs:157-162`.
  `quick_test_simulation` returns *initial* conditions; no integrator is ever invoked in the
  crate. None of the published statistics (survival fractions, 38% high-i excursions) can be
  produced. `ObservabilityCriteria` configured but unused.
- **H4. p9-2017-bias: bias has no ϖ dependence — the bias-corrected null is effectively
  uniform** — `bias_function.rs:30-61`. `perihelion_bias` ignores `_varpi`; the only angular
  dependence is a binary latitude penalty. Brown 2017's bias is dominated by *longitudinal*
  coverage (galactic-plane crossings at λ ≈ 95°/275°, opposition seasons), which is what couples
  bias to ϖ. Also `h_mag` stored but unused (detectability should be H + m_lim via r⁴, not a hard
  r_max = 90 AU).
- **H5. p9-2018-resonance: `OBSERVED_KBO_AXES` scrambled** — `probability_analysis.rs:17-28`.
  E.g. Sedna 472 (vs 506 used in p9-2017-bias), 2013 RF98 780 (actual ~325), "2003 CR105" at
  164 AU. Every implied-a₉ peak is displaced tens of AU. The vetted table in
  `p9-2017-bias/src/kbo_sample.rs` should be the single source (→ p9-core).
- **H6. p9-2018-kuiper-belt: production config numerically inadequate and unrunnable** —
  `simulation.rs:64,141,155`. dt = 300 d with Jupiter direct (14 steps/orbit); 4 Gyr × 400
  particles ≈ 2×10¹² particle-steps; q ∈ (30, 36) AU population crosses Neptune with no
  encounter handling (hybrid integrator unused, `hybrid_changeover_hill` dead). Fix: J2-averaged
  JSU + direct Neptune + P9, dt ≈ P_N/20 ≈ 3000 d, hybrid scheme.

### Medium

- **M1. Resonance membership by semi-major-axis proximity at one snapshot** —
  `simulation.rs:184-186`, `resonance_catalog.rs:209-230`. With the extended catalog the 2%
  windows overlap and everything "matches" something; `is_librating` exists but the pipeline
  never records φ time series. Invalidates `p_simple` / the "<5%" headline.
- **M2. "Farey F5" includes ratios far outside F₅** — `resonance_catalog.rs:98-100`
  (`farey_resonances(5, 30)` admits 29:5, 30:1). Inflating the baseline flattens the contrast
  the crate exists to demonstrate. Also the main test assertion is vacuous
  (`ext_peak < f5_peak * 1.5`).
- **M3. Hidden hard cutoff a₉ ∈ [500, 800] disguised as mass marginalization** —
  `probability_analysis.rs:94-106` (`200 + 30·m̄/m̄·…` = 500, dead `_mass_min/max`); equal kernel
  weight σ = 30 AU per resonance guarantees the plateau "result" partly by construction.
- **M4. p9-2017-bias rejection sampler distorts the null** — `clustering_test.rs:117-129`.
  `weight.min(1.0).max(0.01)`: the floor inflates near-zero-bias directions by orders of
  magnitude; the cap saturates high-bias ones (the weight is a ratio-to-mean and routinely > 1).
  Sample against the raw bias with a true maximum bound.
- **M5. MC resolution can't reach the paper's 0.025%** — `clustering_test.rs:132-134` (1000
  iterations; need ≥ 4×10⁵). Add an `#[ignore]`d long test pinning p_combined.
- **M6. Phase portraits hold i fixed while varying e** — `phase_space.rs:67-90`, violating the
  conserved H_z ∝ √(1−e²)cos i. Parameterize by H_z; also generate the (θ = 2Ω−ϖ−ϖ₉) portrait
  the crate defines but never uses. `find_equilibria` uses the e-step as the Δϖ finite-difference
  step (wrong variable's scale).
- **M7. dt = P_N/8 ≈ 7524 d too coarse for q = 30 AU perihelion passages** —
  `p9-2018-resonance/src/simulation.rs:48-54`; pericenter step artifacts read as spurious chaos
  in a resonance-census crate. dt ≈ 2000–3000 d.
- **M8. Alignment classified from one snapshot** — `p9-2018-kuiper-belt/src/population_analysis.rs:46-57`.
  Circulators pollute both aligned and anti-aligned bins; classify by Δϖ(t) across the snapshots
  already collected.

### Low

- 2000 CR105 (a = 228) in the "a > 230" sample with a test loosened to a > 220 to paper over it.
- Latent unit bug: `precession_rate_j2` returns rad/day where rad/yr documented
  (`p9-2017-dynamics/src/hamiltonian.rs:25-26,76-80`), masked by a 0.0 default.
- Pendulum half-width omits the O(1) eccentricity dependence at e ≈ 0.7–0.9; the
  `CRITICAL_PERIOD_RATIO` is hardcoded rather than derived from overlap.
- Inclination marginalization integrates a sin i prior over [0, π] for a population with
  i ≲ 30°.
- Snapshot trigger `t >= next + interval` skips the first snapshot.
- `jacobi_constant` exists in p9-core and is unused by any crate — a free invariant test.

---

## 4. 2019–2021 crates (clustering, review, oort-cloud, orbit, stability, ztf)

### High

- **H1. p9-2021-stability: Chirikov exponent contradiction** — `chirikov.rs:21`,
  `resonance_chain.rs:45` use exp(−q²/4a_N²); the headline `critical_perihelion`
  (`chirikov.rs:45-52`) is only the K = 1 root if the overlap parameter carries exp(−q²/2a_N²).
  At a = 500 AU: q_crit = 41.4 AU, yet `chirikov_overlap_parameter(500, 45) ≈ 1.5` → `is_chaotic`
  and `stability::classify` disagree about the same orbit. Resonance-width √ of the pendulum
  potential confirms /2 is correct. Factor exp(q²/4a_N²) ≈ 1.4 at q = 35 — not cosmetic.
- **H2. `hansen_x_neg3_0` dimensionally wrong** — `hansen.rs:25-32`: `j·exp(−(1−e)²)` where
  (1−e) = q/a, not the documented q/a_N; off ~4× at e = 0.9 with the opposite e-trend. Public
  API of the disturbing function (currently uncalled).
- **H3. p9-2021-ztf: exclusion collapses to a brightness cut** — `exclusion.rs:36-46`.
  `p_detect = 0.75 × 0.9966` is constant and always > 0.5, so sky coverage and linking efficiency
  have zero effect; `fraction_excluded` is exactly "fraction with mag < 20.5". The paper's 56%
  came from injecting synthetic P9s into actual ZTF epochs. The "paper_exclusion_fraction" test
  fabricates 564/1000 bright objects and asserts 56.4% — a tautology.
- **H4. p9-2021-oort-cloud: the key result (67% vs 88% f_ϖ) is hard-wired** —
  `injection_simulation.rs:96-176`. Eccentricity takes one arbitrary uniform kick over 4 Gyr
  (no Kozai cycles, no ϖ evolution); f_ϖ is computed from *initial* uniform angles (→ 0.50 by
  construction); the scattered-disk control draws Bernoulli with an invented
  `aligned_prob = 0.5 + 0.2·coupling`. Also a sign error: |ϖ−ϖ₉| < 90° labeled "anti-aligned".
  Fix: integrate the secular equations (p9-core machinery) and measure end-state Δϖ.
- **H5. p9-2021-orbit: no likelihood/MCMC; correlations ignored; split-normal sampler wrong** —
  `posterior.rs:27-156`. Marginals transcribed as independent split-normals (m–a, a–q strongly
  correlated in the real posterior); `AsymmetricGaussian::sample` picks sides 50/50 instead of
  σ_l/(σ_l+σ_u). At minimum fix side weights, sample (a, e) jointly, and label the crate a
  posterior-summary emulator.
- **H6. p9-2019-clustering: null omits observational bias — the paper's entire point** —
  `clustering_analysis.rs:67-100`. Angles drawn uniform → this computes "clustering vs isotropy",
  the null the paper's critics used. Consume the p9-2017-bias model when drawing nulls. Secondary:
  doc says f(i) ∝ sin i·exp(−i²/2σ²), code omits sin i; the paper holds each object's (a, e, i)
  fixed and randomizes only angles, the code resamples i globally.

### Medium

- **M1. 2014 FE72 elements wrong** — `kbo_sample.rs:165-173`: e = 0.99 gives q = 21.6 AU
  (Neptune-crossing, violates the q > 30 sample cut); real e ≈ 0.98, q ≈ 36. The test was
  loosened to q > 20 to mask it. Largest-Γ object in the sample → shifts the mean Poincaré vector.
- **M2. Unseeded `rand::random` in the ZTF efficiency MC** — `detection_efficiency.rs:58`.
  Also `self_calibration_correction` is dead physics (never applied) and inconsistent with the
  hard step at 20.5; `linking_threshold: 7` unused (no cadence model exists).
- **M3. No orbit→sky→footprint model in ZTF crate** — flat |β| < 60° cut on caller-supplied
  latitudes; ZTF is declination-limited (δ > −30°). Compute sky positions from the sampled
  elements (p9-core coords).
- **M4. Stability analytics never validated numerically** — `diffusion.rs:96-104` sets
  `d_measured = d_analytical` and the test asserts ratio ≈ 1 (tautology); `hansen_numerical`
  never compared to the approximations it certifies (would have caught H2).
- **M5. `measure_diffusion` is not a diffusion estimator** — `diffusion.rs:27-50`: single
  trajectory, single origin, can't distinguish drift from diffusion (its own test asserts D > 0
  for pure drift). Fit MSD(τ) over time origins / ensemble.
- **M6. p9-2019-review `evaluate_params` returns invented scaling relations** —
  `parameter_survey.rs:166-196` (admitted TODO). All Figure-16-style outputs are decorative;
  quarantine as placeholders or implement the semi-averaged secular integrator.
- **M7. Newton-only Kepler solvers, E₀ = M, e up to 0.99** — `reference_population.rs:119-129`,
  `hansen.rs:57-71` (+3 more copies repo-wide). Use E₀ = π for e > 0.8 or safeguarded Newton;
  consolidate in p9-core.
- **M8. 2021-orbit clustering significance**: textbook Rayleigh vs the paper's bias-aware MC;
  11 hard-coded unsourced ϖ values; test asserts only > 0.95 instead of pinning 0.996 —
  `statistical_measures.rs:107-118`.

### Low

- `standard_map_k = K²` asserted without derivation; ±2 AU "marginal" band arbitrary and
  non-scaling; albedo ~ U(0.3, 0.7) invented (review crate uses 0.40/0.75 — unify photometry in
  p9-core); `A_NEPTUNE`/`M_NEPTUNE_SOLAR` re-hard-coded (5.15e-5 vs core's 5.1513890e-5);
  `ioc_mass_fraction` 0.05√(ρ/ρ₀) invented with a circular test; scattered-disk histogram filled
  with hard-coded counts; `ossos_detection_power = 1/√n` placeholder; across all six crates no
  test reproduces a published number to stated precision.

---

## 5. 2022–2024 crates (des, neptune-crossing, oort-selfgrav, panstarrs)

### High

- **D1. DES: wrong distance law for reflected sunlight** — `recovery_analysis.rs:67-68`.
  Stellar `m = H + 5 log₁₀(d/10)` instead of reflected-light `m ≈ H + 5 log₁₀(rΔ)`
  (flux ∝ 1/r²Δ²). Gives m ≈ 12–15 at 300–1000 AU (naked-eye-adjacent), so completeness
  saturates at ~95% instead of the paper's 87%. The albedo-only H model also lacks radius
  dependence.
- **D2. DES footprint ~2× too large and inconsistent with its own `footprint_area`** —
  `survey_model.rs:34,63-72`. Box cut ≈ 9300 deg² vs the declared (unused) 5000. Exclusion is
  footprint-driven. Use a polygon/HEALPix mask or normalize the cut.
- **N1. Neptune-crossing: TNO table partially fabricated** — `observed_tnos.rs:43-237`.
  "2016 QV89" is a ~30 m NEA (a ≈ 1.7 AU), "2014 LU28" is a retrograde centaur — not the listed
  elements. The zeta statistic is computed from this table; everything downstream is meaningless.
  Internal-consistency tests pass because elements were generated to pass them. Rebuild from the
  paper's actual selection cross-checked against JPL SBDB, with an SBDB-tolerance test.
- **N2. Hypothesis test has no model behind it** — `hypothesis_test.rs:43-61,141-143`.
  CDF(q|r) = (q/r)^1.5 invented (the paper's CDF *is* the simulated perihelion distribution);
  discovery distances fabricated as r = q + 2 + 0.05a; no MC null for zeta → p-value.
- **N3. The "simulation" is a multiplicative fudge** — `simulation.rs:103-126`:
  `i_adjusted = 0.7·i; q_adjusted = 0.85·q` stands in for P9 secular chaos. Configs matching the
  paper (10⁴ particles, 4 Gyr, 5 M⊕/500 AU) are defined but no integration consumes them.
- **O1. Oort-selfgrav: no equations of motion; headline timescale dimensionally invalid** —
  `vzlk.rs:60-114`, `hamiltonian.rs:75-92`. `minimum_perihelion` ignores the Hamiltonian
  (returns the kinematic |J_z| bound, with a dead 500-iteration loop); `evolutionary_timescale_gyr`
  divides days by Gyr-days on a dH/de that is not a frequency (differentiate w.r.t. the canonical
  momentum G = √(GMa)·η). "Timescale ≫ 4.5 Gyr" is the paper's central conclusion; the test
  asserts only τ > 0. Fix: trace level sets of H at fixed (a, J_z), integrate Hamilton's
  equations, test H conservation along trajectories.
- **P1. Pan-STARRS: unseeded RNG in the detection MC** — `detection_pipeline.rs:60`
  (`rand::random`), unlike the seeded sibling crates.

### Medium

- **O2. Galactic tide absent from the oort-selfgrav Hamiltonian** (`hamiltonian.rs:55-59`) —
  a competing quadrupole at a ≈ 10³–10⁴ AU that breaks the axisymmetry `conservation.rs` relies
  on; include the secular tide term or document the validity regime.
- **O3. Likely factor-2 error in H_J2** — `hamiltonian.rs:48-52`: orbit-averaged interior-ring
  quadrupole is GM·(½Σm a²)·(3cos²i−1)/(4a³η³); the code has /8 with the ½ already inside J2eff.
  Tighten the q = 75 AU test from `< 0.5` to the paper's < 2%.
- **O4. Miyamoto–Nagai parameters need verification** (ã = 200 AU with b̃/ã = 5 is
  quasi-spherical with mass peaked inside the planetary region); uniform-M rectangle quadrature
  fragile at e ≥ 0.9 — add an n vs 2n convergence test or integrate in true anomaly.
- **O5. `kozai_constant` is the wrong invariant** — `conservation.rs:63-67` is neither the
  standard C_KL = e²(1 − 5/2 sin²i sin²ω) nor conserved by this Hamiltonian; replace with an
  H-conservation check along integrated trajectories.
- **D3. DES: no multi-epoch/motion/linking model** — single Bernoulli vs the paper's
  k-of-n-night linking over a 6-yr baseline; `n_nights: 575` unused.
- **D4. DES completeness parameters invented, single-band**; color models' `recovery_rate`s are
  stored constants whose only computational path runs through the broken D1 formula.
- **P2. Pan-STARRS conflates declination with ecliptic latitude** —
  `detection_pipeline.rs:53-54,86`: "dec > −30°" applied to β. The boundary is a sinusoid
  spanning β ≈ −53°…−7°; this misclassifies exactly the PS1/DES/ZTF unique-coverage regions.
- **P3. PS1 detection model omits cadence/linking** (`linking_threshold: 9` unused; hard step at
  21.5 vs DES's logistic); `compute_unique_detections` hard-codes DES depth 23.8 while the DES
  crate says 24.1 — two crates disagree about the same survey.
- **N4. Significance claim internally inconsistent** — `lib.rs:6-7` says p = 0.0034 ≈ "5σ";
  the code's own test says 2.7σ. Make doc, constant, and `sigma_rejection` agree with the paper.

### Low

- KS/AD statistics implemented but orphaned (no p-value mapping, nothing feeds them) —
  `hypothesis_test.rs:84-121`.
- Combined-exclusion is declarative addition (0.564 + 0.050 + 0.171 capped) with a stored 0.78 vs
  computed 0.785; the real combination is a per-orbit OR over the reference population.
  ZTF baseline appears as 0.564, 0.56, and 0.612-combined in three places — one shared table.
- `to_p9_params` hard-codes e = 0.3, i = 16°, ω = 150°, Ω = 100° without a source.

---

## 6. 2025 crates (clustering, iras-akari, perturbation)

### High

- **H1. Thermal model radius contradicts its docstring and flips the detectability conclusion** —
  `thermal_model.rs:41-43`. `R_EARTH·M^0.27` = 1.86 R⊕ at 10 M⊕ vs the documented ~3.5 R⊕
  (Fortney et al.). Flux ∝ R² → ~3.5× underestimate; the entire 7–17 M⊕ / 500–700 AU grid then
  falls below the IRAS 0.2 Jy limit, reporting the paper's premise as undetectable. The flux test
  is self-contradictory and vacuous (asserts F < 100 Jy); pin it to the paper's Table 1 values.
- **H2. Proper-motion model physically wrong and parallax-free** — `proper_motion.rs:42-46`.
  Uses total vis-viva speed / r as the sky rate; correct transverse rate is h/r² =
  √(GMa(1−e²))/r². Even a parabolic orbit maxes at √2×circular = 62.8′ at 500 AU, so heliocentric
  motion can never explain the paper's 69.6′ bound — the missing physics is geocentric parallax
  (±6.9′ at 500 AU). Infects `implied_distance` (which returns 371–525 AU instead of the paper's
  500–700). Model geocentric apparent positions at each survey epoch.
- **H3. Three mutually inconsistent forms of the same Neptune-resonance physics** —
  `p9-2025-perturbation/src/resonance.rs:41-56` (exp(−q²/2a_N²)) vs
  `p9-2021-stability/src/resonance_chain.rs:37-45` (exp(−q²/4a_N²)) vs
  `p9-2025-clustering/src/diffusion.rs:14-21` (a third form) — ~1.4× disagreement at q = 35 AU
  plus ~1.6× prefactor differences for the identical 2:j chain. Derive one sourced Hansen
  asymptotic + pendulum width in p9-core; add a cross-crate consistency test.
- **H4. Clustering crate can't produce its paper's results** — lib.rs claims 4 Gyr clone
  integrations but imports only `p9_core::constants`; no integration pipeline, no bias modeling,
  no significance test; results are hard-coded paper numbers and the "bimodal fit" is a PDF
  evaluator with no fitting. Wire clone generation → p9-core integrators → diffusion →
  clustering; add Rayleigh + bias-resampling MC.

### Medium

- **M1. Chance-alignment estimate off ~3 orders of magnitude** — `candidate_search.rs:145-155`:
  uniform-sky estimate gives ~11,300 pairs vs the paper's post-cut 13; no Poisson false-alarm
  probability for the single candidate (the key statistic). Condition on post-cut counts and
  latitude-dependent source density.
- **M2. Flux-ratio window [0.05, 3.0] unphysical** vs the blackbody's own 0.23–0.66 over
  30–50 K — derive bounds from `P9ThermalParams`.
- **M3. Positional uncertainties (20″/30″) carried but never used in cross-matching.**
- **M4. Modified Chirikov overlap is order-dependent** — `chirikov.rs:14-19`:
  `overlap(r1,r2) ≠ overlap(r2,r1)`; results depend on chain iteration direction.
- **M5. Farey machinery is a stub** (integer count, not mediant-tree enumeration — the fine comb
  is the paper's point); promised 4:j chain and l = 4 coefficients absent;
  `spherical_harmonics.rs` contains two hard-coded coefficient pairs, no harmonics;
  `total_disturbing_function` omits all Hansen/eccentricity factors; octupole Hansen indexing
  by `enumerate()` silently wrong for tables not starting at j = 0.
- **M6. Stability boundary `max()` of three regime curves** makes the fine comb (weak, bounded
  chaos) the binding instability boundary at a = 100–300 AU — implement piecewise per the
  paper's figure; verify ln vs log10 in the −63.1 + 19.4·log(a) fit.
- **M7. Hansen quadrature has no convergence control** (needed n grows ~ j(1−e)^{−3/2}; 256
  points at e > 0.9 with large j aliases); Newton Kepler from E₀ = M fragile at e → 0.999.
- **M8. Diffusion estimator biased** — max over t of (Δa)²/t overestimates (sup of a random walk
  grows as t log log t); MAD helper actually computes mean absolute deviation.
- **M9. Clone generation: no element covariances, angles never perturbed, silent class changes
  after clamps, 7 hard-coded TNOs vs the paper's ~2× larger sample.**

### Low

- Duplications → p9-core: `hansen_numerical` + Kepler solver (verbatim twins), Vincenty
  separation 2×, circular-PM inversion 3×, Neptune constants re-defined 4× (one diverging),
  giant-planet table inline.
- `is_detectable` hard-codes survey sensitivities duplicating the survey structs.
- MEGNO classifier exists but no MEGNO is ever computed.
- Pole misalignment computed as a 2D angle in (i cos Ω, i sin Ω) rather than spherical distance;
  the referenced Laplace plane never computed.
- Blackbody assumption fine vs the paper, but add an emissivity hook (H₂ CIA deviations at 40 K).

---

## 7. Cross-cutting statistics & reproducibility

- **Bias-aware nulls are the repo's biggest statistical gap**: p9-2017-bias's bias function has
  no ϖ dependence (H4 §3), and p9-2019-clustering / p9-2021-orbit test against uniformity — the
  null the published papers explicitly reject as inadequate. Fixing the bias function once and
  consuming it from the clustering crates fixes three crates.
- **Unseeded RNG** in p9-2021-ztf (`detection_efficiency.rs:58`) and p9-2024-panstarrs
  (`detection_pipeline.rs:60`); everywhere else seeds are threaded but never recorded in
  serialized outputs.
- **Rayleigh machinery** is re-implemented ≥ 6× with an approximation that's poor at n = 6;
  one `p9_core::analysis::circular` module (circular mean/std, R̄, Rayleigh p with small-n
  correction, Kuiper test) fixes all call sites.
- **Tautological tests** are endemic: tests that fabricate inputs to equal the published number
  (ZTF 56.4%, diffusion ratio 1.0, oort 5%), assert `is_finite()`, or verify a power law against
  a function defined as that power law. The single highest-leverage testing change is **one
  paper-number regression test per crate** plus **physical-invariant tests** (energy/L for
  integrators, H/H_z conservation for secular trajectories, element round-trips including
  hyperbolic and e≈0/i≈0 — the hyperbolic round trip alone would have caught the p9-core sign
  error).

## 8. Consolidation into p9-core (deduplication map)

| What | Current copies | Target |
|------|----------------|--------|
| Kepler solver | p9-core types.rs + kepler_step.rs, p9-2021-orbit, p9-2021-stability, p9-2024-oort-selfgrav, p9-2025-perturbation | one safeguarded solver in p9-core |
| Circular statistics / Rayleigh | ≥6 inline copies across 2016/2017/2018/2019/2021 crates | `p9_core::analysis::circular` |
| Hansen coefficients | p9-2021-stability, p9-2025-perturbation (verbatim twins) | `p9_core::analysis::hansen` with convergence control |
| Resonance width / Chirikov overlap | 3 inconsistent forms (2021-stability, 2025-perturbation, 2025-clustering) | one sourced implementation + consistency test |
| ETNO observational table | p9-2017-bias (vetted), p9-2018-resonance (corrupted), p9-2019-clustering, p9-2021-orbit | one SBDB-checked table with provenance |
| J2-effective giants | p9-core (canonical), p9-2017-dynamics, p9-2024-oort-selfgrav | p9-core, with planet-set parameter + GM boost |
| Photometry (H, apparent mag, albedo) | p9-2016-constraints, p9-2019-review, p9-2021-orbit, p9-2022-des, p9-2024-panstarrs | `p9_core::analysis::photometry` |
| Survey exclusion fractions | 0.564 / 0.56 / 0.612 in three crates | one shared table |
| Neptune constants | 4 redefinitions, one diverging (5.15e-5) | p9-core constants |
| Simulation snapshot driver | 3 copies with inconsistent snapshot triggers | one driver in p9-core |
| Resonant-angle definition | 2 inconsistent conventions (2017-dynamics, 2018-resonance) | one definition in p9-core |
