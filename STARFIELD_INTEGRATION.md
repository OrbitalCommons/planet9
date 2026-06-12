# Starfield ↔ Planet9 Integration Notes

Cross-map between this repo and `starfield` v0.13 (~/repos/starfield, the Rust
skyfield port) plus its `starfield-datasources` workspace. Two directions:
where starfield tooling would improve planet9's modeling, and where planet9
machinery would be generally useful upstream in starfield.

---

## Part 1 — Planet9 improvements available from starfield

Ordered by leverage.

### 1.1 Coordinate frames: replace every hand-rolled transform

The single biggest win. Planet9 currently carries **four local coordinate
conversions**, all flagged in REPRODUCTION_NOTES.md as consolidation debt:

| Planet9 location | What it hand-rolls | Starfield replacement |
|---|---|---|
| `p9-2022-des/src/sky.rs` (shared with p9-2024-panstarrs) | ecliptic→equatorial for declination cuts | `framelib::ECLIPTIC_J2000` / `ICRS` frame matrices |
| `p9-2021-ztf/src/sky.rs` | elements→ecliptic λ/β→declination | same |
| `p9-2025-iras-akari/src/orbital_constraints.rs` | local equatorial↔ecliptic | same |
| `p9-core/src/forces/galactic_tide.rs` | galactic-pole orientation (~60.2° literal) | `framelib::GALACTIC` matrix (SPICE-derived) |

Starfield's transforms come with IAU 2006 precession (`precessionlib`) and
full IAU 2000A nutation (1,365 terms, `nutationlib`), validated against
skyfield/astropy to ~1e-8 per matrix element. The survey-footprint findings
(P2: "declination/ecliptic-latitude conflation") exist precisely because
planet9 lacked this; `coordinates::{Equatorial, Ecliptic, Galactic}` with
`From` impls make those bugs unrepresentable.

Also relevant: `p9-2025-clustering`'s orbital-pole / Laplace-plane geometry
could use `Cartesian3::angular_distance` and the frame machinery instead of
its local spherical-distance helper.

### 1.2 Earth position & geocentric apparent positions (IRAS/AKARI, DES, PS1)

`p9-2025-iras-akari` models geocentric parallax with a **documented
circular-Earth approximation**; `p9-2022-des/sky.rs` approximates parallactic
motion similarly. Starfield gives the exact chain:

- `Loader::open("de421.bsp")` → `SpiceKernel` (auto-downloaded, mmap-backed)
- `kernel.at("earth", &t)` → barycentric Earth at the true survey epochs
- `Position::observe(...)` → `.apparent(...)` adds light-time iteration,
  relativistic aberration, and gravitational deflection — all currently
  ignored in planet9's survey models. At 500–700 AU, light-time alone is
  ~3–4 days; for the 1983→2006 IRAS/AKARI epoch pair that's a real (if small)
  systematic on the proper-motion window.

The epochs themselves should become `time::Time` values (TT/TDB-aware)
instead of raw year floats — `Timescale::utc((1983, ...))` etc.

### 1.3 Giant-planet states from real ephemerides

`p9-core/src/initial_conditions/planets.rs` hard-codes J2000 state vectors.
`planetlib::Ephemeris::get_state(Body::Neptune, &t)` (DE440/441) gives exact
states at **any** epoch — better initial conditions, exact Laplace-plane
references, and the ability to start integrations at survey/discovery epochs
rather than J2000 only.

### 1.4 Kepler solver and element conversions

- p9-core's `kepler_drift` universal-variable iteration has a **known
  convergence edge case** (a≈800 AU, e≈0.95, long dt — see
  REPRODUCTION_NOTES.md "Known numerical issues"). Starfield's `keplerlib`
  carries a SPICE `prop2b` port (universal variables + Stumpff) and an
  eccentric-anomaly solver (elliptic *and* hyperbolic, arXiv:2108.03215),
  both parity-tested against skyfield. Replacing or cross-validating
  p9-core's drift against it would close that issue.
- `elementslib::OsculatingElements` is a robust element extractor handling
  elliptic/parabolic/hyperbolic with explicit `inf` semantics — exactly the
  edge cases (e≈0, i≈0, retrograde, hyperbolic) where p9-core's
  `cartesian_to_elements` had bugs. Worth adopting or at least
  property-testing against.

### 1.5 Live small-body data instead of hand-vetted tables

The fabricated-TNO-table finding (N1: a near-Earth asteroid listed as a
110 AU TNO) happened because planet9 transcribes elements by hand.
Starfield's `SbdbClient::query(SbdbQueryParams)` can express the
p9-2024-neptune-crossing selection directly (a > 100 AU, q < 30 AU,
i < 40°, numbered/multi-opposition) and return `SmallBodyOrbit` records with
epochs, arcs, and condition codes. `starfield-mpc` parses MPCORB for the
same offline. The right architecture: `p9_core::data::etno` keeps the
frozen, provenance-commented snapshot used by regression tests, plus an
optional refresh path through SbdbClient that diffs against the snapshot —
catching both transcription errors and element drift.

`HorizonsClient::generate_spk_kernel` also covers the discovery-circumstance
gap (N2's documented r ≈ 1.15q approximation): real discovery-epoch
distances for each sample object.

### 1.6 Photometry cross-validation

- `starfield-mpc::hg_apparent_magnitude(h, g, phase)` is the proper IAU H–G
  phase law; p9-core's `apparent_magnitude` omits the phase function
  (negligible at opposition for P9, but the DES/PS1 models integrate over
  epochs where it isn't exactly zero).
- `magnitudelib::planetary_magnitude` (Mallama & Hilton 2018) covers
  Neptune — the natural validation anchor for
  `photometry::planet_apparent_magnitude`, which is already tested against
  Neptune's H ≈ −6.87 with our own formula. An independent cross-check test
  would pin both.

### 1.7 Event finding for survey models

`searchlib::{find_discrete, find_maxima}` plus
`almanac::oppositions_conjunctions` would let the DES/PS1/ZTF models compute
**when** a synthetic P9 is at opposition / inside a footprint, instead of
quantizing onto night grids by hand. Lower priority — the current
k-of-n-night models are adequate — but it's the principled version of the
`night_usability` calibration knob.

### 1.8 Shared constants

`starfield::constants` (AU_KM = 149,597,870.700, GM_SUN in km³/s², C_AUDAY,
J2000) and p9-core `constants.rs` should agree to full precision; today they
are independently transcribed. If p9-core takes starfield as a dependency,
re-export rather than redefine.

### Suggested adoption order

1. Add `starfield` as a p9-core dependency; replace the three local `sky.rs`
   conversions + galactic-tide orientation with `framelib` (deletes code,
   fixes nothing visible — pure de-risking).
2. Cross-validate / replace `kepler_drift` and `cartesian_to_elements`
   against `keplerlib`/`elementslib` (closes the known numerical issue).
3. Earth-from-ephemeris in IRAS/AKARI parallax + survey epoch handling via
   `time::Time`.
4. SBDB refresh path for `data::etno` with snapshot diffing.
5. Giant-planet initial conditions from `planetlib` at arbitrary epochs.

---

## Part 2 — Planet9 machinery generally useful upstream in starfield

Starfield's survey confirms: **no N-body or perturbed propagation exists
anywhere in it** (two-body Kepler, SGP4, SPK lookup only), and nothing in the
statistics/dynamics space. Planet9 built exactly that. Candidates, ordered by
generality:

### 2.1 N-body integrator stack → `starfield-nbody` (or `starfield::nbodylib`)

The headline contribution. p9-core now has:
- symplectic Wisdom–Holman in democratic-heliocentric coordinates (DLL98),
  with energy/angular-momentum diagnostics and dt-scaling tests
- Bulirsch–Stoer with step-size control
- Chambers (1999) hybrid switching with smooth changeover and
  minimum-separation encounter prediction
- the `ExtraForce` composable perturbation hook (J2 secular field, galactic
  tide, custom closures)

This fills starfield's documented gap and is already validated by 600+ tests.
As a companion crate (`starfield-nbody` in the datasources-style pattern) it
could consume starfield's `Time`, frames, and ephemerides for initial
conditions — the integration is symmetric with Part 1.3.

### 2.2 Secular / resonance dynamics toolkit

`p9_core::analysis::{secular, resonance, hansen}`:
- quadrupole + octupole secular Hamiltonians and numerical Gauss-ring double
  averaging (valid in the non-hierarchical regime where series diverge)
- Hansen coefficients with convergence-controlled quadrature
- canonical mean-motion resonance widths, Chirikov overlap, critical
  perihelion, resonant-angle/libration detection on time series

No equivalent exists in starfield or, to my knowledge, anywhere in the Rust
ecosystem. Useful to anyone doing long-term small-body dynamics (Kozai/vZLK,
resonance sticking, Oort-cloud work).

### 2.3 Circular statistics

`p9_core::analysis::circular` — circular mean/std, mean resultant length,
Rayleigh test with the small-n series correction, Kuiper test. Generic
directional statistics with astronomy-grade tests; natural fit for a
starfield statistics module (skyfield has no equivalent either — this would
be a differentiator).

### 2.4 Generic small-body photometry

p9-core's `photometry` (mass-radius relation for Neptunian bodies,
H from radius+albedo, reflected-light apparent magnitude) complements
starfield's `magnitudelib`, which only covers the 8 known planets via
Mallama–Hilton. A "hypothetical body" API (mass/radius/albedo → V at r, Δ)
is useful for occultation prediction, survey planning, and population
synthesis. Combine with `starfield-mpc`'s H–G phase law for the full picture.

### 2.5 Survey detection/completeness modeling

The per-orbit survey machinery built for ZTF/DES/PS1 — footprint membership,
per-epoch magnitude-dependent efficiency, k-of-n-night linking via exact
Poisson-binomial tails, per-orbit OR combination across surveys — generalizes
to any "would survey X have seen object Y" question. Fits the
`starfield-datasources` pattern (alongside `starfield-rubin`) as a
`starfield-surveysim` crate. The Poisson-binomial tail helper alone
(`p9-2022-des`) is a reusable numerical primitive.

### 2.6 Synthetic small-body populations

`p9_core::initial_conditions` (scattered-disk generators with documented
distributions, deterministic-N resampling, seeded) complements starfield's
synthetic *star* catalogs — same idea, solar-system flavored. Useful for
injection/recovery testing against any of the datasource crates.

### 2.7 Smaller pieces

- MSD-based diffusion estimator with multi-origin averaging
  (`p9-2025-clustering/src/diffusion.rs`) — time-series analysis primitive.
- The pyo3 parity-testing pattern flows the other way: starfield's
  `python-tests` harness (validating against skyfield/astropy) is a model
  planet9 should eventually adopt for validating the integrators against
  rebound.

---

## Non-overlaps (no action)

- Starfield's star catalogs (Gaia/Hipparcos), SGP4/satellite work, eclipses,
  almanac, image/star-finder tooling have no planet9 counterpart or need.
- Planet9's paper-specific crates (everything outside p9-core) are
  reproductions, not library material; only p9-core contents are upstream
  candidates.
