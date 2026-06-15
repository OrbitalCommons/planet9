//! Recovery analysis for Planet Nine in the Dark Energy Survey.
//!
//! Injects a seeded reference population of synthetic Planet Nine orbits
//! (Brown & Batygin 2021 posterior, via `p9_2021_orbit`) through the DES
//! detection model: equatorial footprint membership of the *apparent*
//! per-night positions, band magnitudes from the surface color model, and
//! the paper's multi-night linking requirement. Paper results (Belyakov,
//! Bernardinelli & Brown 2022): of 100,000 synthetic orbits, 11,709 cross
//! the wide DES footprint, of which 10,187 (87.0%) are recovered; DES
//! uniquely excludes an additional 5.0% of parameter space beyond ZTF.

use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use p9_core::analysis::surveys::combined_unique_exclusion;
use p9_core::coords::observer::{EarthProvider, EarthState};
use p9_core::data::reference_population::{generate_reference_population, heliocentric_distance};
use p9_core::types::{OrbitalElements, P9Params};

use crate::color_models::{self, ColorModel};
use crate::survey_model::{DesBand, DesSurvey};

/// Result of the DES recovery simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecoveryResult {
    /// Number of synthetic orbits crossing the DES footprint
    pub n_crossing: usize,
    /// Number successfully recovered
    pub n_recovered: usize,
    /// Recovery fraction
    pub recovery_frac: f64,
}

impl RecoveryResult {
    /// Paper's reported result: 10187 / 11709 = 87.0%.
    pub fn paper_result() -> Self {
        Self {
            n_crossing: 11709,
            n_recovered: 10187,
            recovery_frac: 0.870,
        }
    }
}

/// Simulate DES recovery for a synthetic Planet Nine population with the
/// fiducial (BB21) color model.
pub fn compute_recovery(n_synthetic: usize, seed: u64) -> RecoveryResult {
    compute_recovery_for_model(&color_models::fiducial(), n_synthetic, seed)
}

/// Simulate DES recovery under a given surface color model.
///
/// For each synthetic orbit drawn from the Brown & Batygin (2021) posterior
/// emulator:
/// 1. the apparent equatorial position is evaluated at every DES observing
///    epoch (parallax + orbital drift included) for footprint membership;
/// 2. band magnitudes follow from the model's albedo/colors and the
///    p9-core reflected-light photometry at the orbit's heliocentric
///    distance;
/// 3. recovery is a Bernoulli draw at the exact Poisson-binomial
///    probability of detections on >= `min_nights` distinct nights.
pub fn compute_recovery_for_model(
    model: &ColorModel,
    n_synthetic: usize,
    seed: u64,
) -> RecoveryResult {
    run_recovery(model, n_synthetic, seed, None)
}

/// Ephemeris-path variant of `compute_recovery_for_model`: per-night
/// apparent positions from real Earth states (DE421 via `EphemerisEarth`,
/// or any other explicit `EarthProvider`). The night grid's Earth states
/// are computed once and shared across all orbits.
pub fn compute_recovery_for_model_with_earth(
    model: &ColorModel,
    n_synthetic: usize,
    seed: u64,
    earth: &mut impl EarthProvider,
) -> RecoveryResult {
    let nights = DesSurvey::default().observing_epochs_with_earth(earth);
    run_recovery(model, n_synthetic, seed, Some(&nights))
}

fn run_recovery(
    model: &ColorModel,
    n_synthetic: usize,
    seed: u64,
    nights: Option<&[(DesBand, f64, EarthState)]>,
) -> RecoveryResult {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let survey = DesSurvey::default();

    let population = generate_reference_population(n_synthetic, &mut rng);

    let mut n_crossing = 0usize;
    let mut n_recovered = 0usize;

    for obj in &population {
        let elements = OrbitalElements {
            a: obj.a,
            e: obj.e,
            i: obj.i,
            omega: obj.omega,
            omega_big: obj.omega_big,
            mean_anomaly: obj.mean_anomaly,
        };

        let crosses = match nights {
            Some(nights) => survey.crosses_footprint_with_earth(&elements, nights),
            None => survey.crosses_footprint(&elements),
        };
        if !crosses {
            continue;
        }
        n_crossing += 1;

        let r_au = heliocentric_distance(&P9Params {
            mass_earth: obj.mass,
            a: obj.a,
            e: obj.e,
            i: obj.i,
            omega: obj.omega,
            omega_big: obj.omega_big,
            mean_anomaly: obj.mean_anomaly,
        });
        let mags = model.band_magnitudes(obj.mass, r_au);

        let p_detect = match nights {
            Some(nights) => {
                survey.detection_probability_for_orbit_with_earth(&elements, &mags, nights)
            }
            None => survey.detection_probability_for_orbit(&elements, &mags),
        };
        if rng.gen::<f64>() < p_detect {
            n_recovered += 1;
        }
    }

    let recovery_frac = if n_crossing > 0 {
        n_recovered as f64 / n_crossing as f64
    } else {
        0.0
    };

    RecoveryResult {
        n_crossing,
        n_recovered,
        recovery_frac,
    }
}

/// Fraction of P9 parameter space uniquely excluded by DES beyond ZTF
/// (published value from the shared p9-core survey table: 5.0%,
/// Belyakov et al. 2022 via Belyakov et al. 2024).
pub fn unique_exclusion() -> f64 {
    p9_core::analysis::surveys::EXCLUSION_FRACTIONS
        .iter()
        .find(|s| s.survey == "DES")
        .expect("DES in shared exclusion table")
        .unique_fraction
}

/// Combined ZTF + DES exclusion of the Planet Nine parameter space, from
/// the shared p9-core unique fractions (0.564 + 0.050 = 0.614; the paper
/// rounds to 61.2%).
pub fn combined_ztf_des() -> f64 {
    let ztf = p9_core::analysis::surveys::EXCLUSION_FRACTIONS
        .iter()
        .find(|s| s.survey == "ZTF")
        .expect("ZTF in shared exclusion table")
        .unique_fraction;
    combined_unique_exclusion(&[ztf, unique_exclusion()])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paper_result_values() {
        let r = RecoveryResult::paper_result();
        assert_eq!(r.n_crossing, 11709);
        assert_eq!(r.n_recovered, 10187);
        assert!((r.recovery_frac - 0.87).abs() < 0.01);
    }

    /// X1/D1 regression: the computed recovery of footprint-crossing
    /// synthetic Planet Nines must reproduce the paper's 87.0% (the old
    /// stellar distance law saturated completeness at ~95%).
    #[test]
    fn computed_recovery_matches_paper() {
        let result = compute_recovery(6000, 2022);
        assert!(
            result.n_crossing > 200,
            "expected a usable crossing sample, got {}",
            result.n_crossing
        );
        assert!(
            (result.recovery_frac - 0.870).abs() < 0.04,
            "recovery = {:.3} vs paper 0.870",
            result.recovery_frac
        );
    }

    /// X1/D1 regression on the ephemeris path: the 87.0% pin must also
    /// hold with real DE421 Earth states (light-time + aberration); the
    /// real-Earth parallax differs from the analytic circular model by
    /// < 0.1 deg at P9 distances, so recovery moves well inside the
    /// shared 0.04 tolerance. Kernel-gated: skips without the cache file.
    #[test]
    fn computed_recovery_matches_paper_with_ephemeris_earth() {
        let mut earth = match p9_core::coords::observer::EphemerisEarth::try_new() {
            Ok(e) => e,
            Err(_) => {
                eprintln!("skipping: no ephemeris kernel in starfield cache");
                return;
            }
        };
        let result = compute_recovery_for_model_with_earth(
            &color_models::fiducial(),
            6000,
            2022,
            &mut earth,
        );
        assert!(
            result.n_crossing > 200,
            "expected a usable crossing sample, got {}",
            result.n_crossing
        );
        assert!(
            (result.recovery_frac - 0.870).abs() < 0.04,
            "ephemeris-path recovery = {:.3} vs paper 0.870",
            result.recovery_frac
        );
        // The analytic fallback agrees closely (same seed, same
        // population; only the Earth model differs).
        let analytic = compute_recovery(6000, 2022);
        assert!(
            (result.recovery_frac - analytic.recovery_frac).abs() < 0.02,
            "ephemeris {:.3} vs analytic {:.3}",
            result.recovery_frac,
            analytic.recovery_frac
        );
    }

    /// The crossing fraction is footprint-driven: ~11.7% of the reference
    /// population (11,709 / 100,000) crosses the 5000 deg^2 footprint.
    #[test]
    fn crossing_fraction_near_paper() {
        let result = compute_recovery(6000, 2022);
        let frac = result.n_crossing as f64 / 6000.0;
        assert!(
            (frac - 0.117).abs() < 0.06,
            "crossing fraction = {frac:.3} vs paper ~0.117"
        );
    }

    #[test]
    fn recovery_ordered_by_surface_brightness() {
        // Computed recoveries must order with detectability: the bright CH4
        // surface (albedo 0.75) beats the dark ultra-red Super-KBO surface
        // (albedo 0.1, faintest in the deep g band), as in Table 1
        // (88.1% vs 76.8%).
        let bright =
            compute_recovery_for_model(&color_models::methane_40k(), 3000, 7).recovery_frac;
        let dark = compute_recovery_for_model(&color_models::super_kbo(), 3000, 7).recovery_frac;
        assert!(
            bright > dark,
            "CH4 recovery {bright:.3} should exceed Super-KBO {dark:.3}"
        );
    }

    #[test]
    fn computed_recoveries_near_table1() {
        // Coarse per-model regression against Table 1 (the single
        // calibrated night-usability parameter is shared across models, so
        // inter-model differences are genuinely computed from colors).
        for model in color_models::all_models() {
            let r = compute_recovery_for_model(&model, 3000, 11).recovery_frac;
            assert!(
                (r - model.paper_recovery).abs() < 0.12,
                "{}: computed {r:.3} vs Table 1 {:.3}",
                model.name,
                model.paper_recovery
            );
        }
    }

    #[test]
    fn unique_exclusion_matches_shared_table() {
        assert!((unique_exclusion() - 0.05).abs() < 1e-10);
    }

    #[test]
    fn combined_exclusion_matches_paper() {
        // Shared-table combination 0.614 vs the paper's published 61.2%.
        assert!((combined_ztf_des() - 0.612).abs() < 0.005);
    }
}
