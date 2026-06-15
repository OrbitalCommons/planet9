# Reuse Consolidation Plan — lift cross-crate shared code into `p9-core`

Audit of every place a paper-crate depends on **another paper-crate** (not
`p9-core`). Goal: move the genuinely shared primitive into `p9-core` (with its
unit tests), then repoint **both** the original producer and every consumer at
the core copy, and drop the inter-crate `path` dependency.

24 inter-crate edges found, grouped below by the proposed `p9-core` home.
Producer = the crate that currently owns the code; Consumers = crates that
import it. "Lift" = the item(s) to move into core.

---

## G1. Reference Planet-Nine population  → `p9_core::data::reference_population`
The single biggest shared dependency (6 consumers).

- **Producer:** `p9-2021-orbit::reference_population` — `generate_reference_population`, `heliocentric_distance`.
- **Consumers:** `p9-2018-wise-search`, `p9-2021-ztf`, `p9-2022-des`, `p9-2024-panstarrs`, `p9-2025-ps1-holman`, `p9-2023-lsst-strategy`.
- **Lift:** the seeded synthetic P9 (mass, a, e, i, ω, Ω) population generator + the helio-distance helper.
- **Tests:** move the generator's unit tests to core; keep each survey crate's exclusion-fraction tests (they consume the population).

## G2. Thermal / blackbody IR–mm model  → `p9_core::analysis::thermal`
Duplicated across the entire IR/mm/thermal family.

- **Producers:** `p9-2018-wise-search::thermal_model` (`P9Thermal`, `INTERNAL_TEMP_FLOOR_K`, `W1_WAVELENGTH_M`, `W1_ZERO_POINT_JY`), `::survey_model` (`WiseSurvey`, `W1_DEPTH_REFERENCE`), `::detectability::max_detectable_distance`, `::sky::galactic_latitude_deg`; `p9-2025-iras-akari::thermal_model` (`P9ThermalParams`).
- **Consumers (path-dep):** `p9-2016-wise-coadd`, `p9-2021-act-mm` (← wise-search); `p9-2022-iras-candidate` (← iras-akari).
- **Also duplicate (no path-dep, but reimplement the same Planck/effective-temperature math):** `p9-2016-fortney-thermal`, `p9-2016-cowan-thermal`, `p9-2016-linder-evolution`, `p9-2025-simons-forecast`, `p9-2025-akari-refutation`.
- **Lift:** Planck `B_ν(T)` flux-density, effective-temperature energy balance (internal floor + insolation), band flux / Vega-zero-point→magnitude, `max_detectable_distance` bisection. One `thermal` module retires ~6 private copies.
- **Tests:** consolidate the Planck/T_eff/zero-point tests into core; crates keep their band-specific reference-value pins.

## G3. Brown & Batygin orbit geometry + favored-ν ephemeris  → `p9_core::types` (geometry) + `p9_core::data::ephemeris_constraint`
- **Producer:** `p9-2016-cassini-ranging` — `brown_batygin_orbit`, `geometry::{p9_position_at_true_anomaly, P9Geometry}`, `perturbation::*`, `published::{FAVORED_INTERVAL_DEG, TRUE_ANOMALY_TOLERANCE_DEG}`.
- **Consumers:** `p9-2016-holman-payne`, `p9-survey`.
- **Lift:** `p9_position_at_true_anomaly` (thin wrapper over `elements_to_cartesian` at a given ν) + `P9Geometry` → core types; the favored true-anomaly interval/tolerance constants → a small shared `ephemeris_constraint` data module. Leave the Cassini *range-residual* physics in the cassini crate.

## G4. Perturber-induced planet perihelion precession  → `p9_core::analysis::secular`
- **Producer:** `p9-2026-iorio-precession` — `Perturber`, the Δϖ̇ ∝ GM_p/a_p³ secular-precession formula, `quadrupole_prefactor_consistency`.
- **Consumers:** `p9-2016-iorio-precession` (re-exports it wholesale).
- **Lift:** the planet-precession function into `analysis::secular` (it already hosts `coplanar_quadrupole`); both Iorio crates apply the *bounds*.

## G5. Giant-planet angular-momentum / wire model  → `p9_core::initial_conditions`
- **Producer:** `p9-2016-obliquity::solar_model` — `GIANT_PLANETS` (mass/element/L wire set) + angular-momentum helpers.
- **Consumers:** `p9-2016-obliquity-gomes`, `p9-2016-obliquity-lai`.
- **Lift:** the giant-planet wire set + total-angular-momentum helper. Note `p9_core::initial_conditions::planets` already hard-codes giant-planet states — **merge these two giant-planet sources into one** so they can't drift.

## G6. Survey detection model: footprints, ZTF/DES depth, detection prob, sky position  → `p9_core::analysis::surveys` (+ `p9_core::coords`)
- **Producers:** `p9-2021-ztf` (`detection_efficiency::detection_probability_for_orbit`, `survey_model::ZtfSurvey`); `p9-2022-des` (`sky::{apparent_position_deg, apparent_position_with_earth_deg, phase_angle_with_earth}`, `color_models`, `survey_model::{DesBand, DesSurvey, poisson_binomial_tail}`); `p9-survey` (`schema::Footprint`, `telescope::catalog`).
- **Consumers:** `p9-2024-panstarrs` (← ztf + des), `p9-2025-ps1-holman` (← des::sky), `p9-2023-lsst-strategy` (← survey Footprint/telescope).
- **Lift:** the apparent-sky-position-from-elements helpers and `phase_angle_with_earth` → `coords`; the generic logistic detection-probability and `poisson_binomial_tail` → `analysis::surveys`; a shared `Footprint` type. (The `analysis::surveys` table already centralises depths.)

## G7. Reflex / annual parallax  → `p9_core::coords::observer`
- **Producer:** `p9-survey::parallax` — `annual_parallax_arcsec`, `parallax_arcsec`, `AU_KM`.
- **Consumers:** `p9-2025-parallax-search`.
- **Lift:** baseline→parallax helpers (`AU_KM` duplicates a literal already near `constants`).

## G8. Matched-filter stacking + orbit-search metric  → `p9_core::analysis::stacking`
- **Producer:** `p9-2025-stacking::{matched_filter, orbit_metric}`.
- **Consumers:** `p9-2020-tess-shiftstack`.
- **Lift:** √N stacked-SNR + 1.25·log₁₀N depth gain + the rate-error metric / trial-orbit count (generic survey-detection math, not paper-specific).

## G9. Von Mises special functions  → `p9_core::analysis::circular`
- **Producer:** `p9-2025-clustering::clustering` — `bessel_i0`, `kappa_from_r_bar`.
- **Consumers:** `p9-2026-apsidal-clustering`.
- **Lift:** `bessel_i0` and the A⁻¹(R̄) von-Mises κ estimator into `analysis::circular` (which already owns the circular statistics).

## G10. IRAS/AKARI candidate-pair geometry  → `p9_core::coords` (candidate-pair helpers)
- **Producer:** `p9-2025-iras-akari` — `proper_motion::transverse_rate_arcmin_yr`, `orbital_constraints::*`, `survey_model::{IrasSurvey, AkariFisSurvey}`.
- **Consumers:** `p9-2025-akari-refutation`.
- **Lift:** the two-epoch transverse-rate / parallax-budget helpers (claim and refutation crates share one geometry). Survey-epoch tables can stay crate-side.

## G11. Stable-KBO table + scattered-disk simulator  → `p9_core::data` + `p9_core::initial_conditions`
- **Producer:** `p9-2016-evidence` — `kbo_elements::{stable_kbos, longitude_of_perihelion}`, `scattered_disk_sim::{run_scattered_disk, DiskSimConfig}`.
- **Consumers:** `p9-2016-constraints`.
- **Lift:** the stable-KBO element table (overlaps `data::etno` — reconcile) and the scattered-disk generator (overlaps `initial_conditions::scattered_disk` — reconcile).

## G12. Siraj 2024 reference constants  → `p9_core::data` (minor)
- **Producer:** `p9-2024-siraj-orbit::reference`.
- **Consumer:** `p9-survey`.
- **Lift:** small set of labelled constants; low priority (a handful of `const`s).

---

## Suggested execution order (highest leverage first)
1. **G1** reference population (6 consumers) and **G2** thermal model (≈8 crates) — the big wins.
2. **G6** survey/sky detection helpers (3 consumers + the survey table).
3. **G9** von Mises, **G8** stacking, **G7** parallax, **G4** precession — clean single-edge lifts.
4. **G3** orbit-geometry/ephemeris, **G10** IRAS/AKARI geometry, **G5** giant-planet model (also de-dups an existing core copy), **G11** KBO/scattered-disk (reconcile with existing core modules).
5. **G12** Siraj constants — last, trivial.

Each step: move code + tests into `p9-core`, repoint producer + consumers, delete the now-unused inter-crate `path` dep, run `cargo build/test/fmt` for the whole workspace, one PR per group.
