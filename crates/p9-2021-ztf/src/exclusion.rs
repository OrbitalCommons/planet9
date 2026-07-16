//! Parameter-space exclusion analysis for Planet Nine.
//!
//! Given a reference population of synthetic Planet Nine realizations drawn
//! from the posterior (Brown & Batygin 2021), each orbit's probability of
//! having been detected by ZTF is computed from its sky position
//! (declination footprint), its apparent magnitude (logistic recovery
//! efficiency) and the linking efficiency. Since no Planet Nine was found,
//! the *expected exclusion fraction* is the mean of these per-orbit
//! detection probabilities.
//!
//! (The previous implementation collapsed to a pure brightness cut: a
//! constant 0.75 × 0.9966 "probability" exceeded its 0.5 threshold for every
//! object brighter than 20.5, so footprint and linking had zero effect, and
//! the headline test fabricated 564 bright objects out of 1000 to "match"
//! the paper's 56.4%.)

use serde::{Deserialize, Serialize};

use p9_core::types::OrbitalElements;

use crate::detection_efficiency::detection_probability_for_orbit;
use crate::survey_model::ZtfSurvey;

pub use p9_core::analysis::surveys::{combined_unique_exclusion, EXCLUSION_FRACTIONS};

/// Result of the exclusion analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusionResult {
    /// Expected fraction of the prior parameter space excluded (paper: 56.4%)
    pub fraction_excluded: f64,
    /// Total number of prior draws evaluated
    pub n_total: u32,
    /// Expected number of draws that would have been detected
    pub expected_excluded: f64,
}

/// Compute the expected fraction of the Planet Nine prior excluded by the
/// ZTF non-detection, from per-orbit detection probabilities accumulated
/// deterministically over the (seeded, caller-generated) reference
/// population of `(elements, apparent magnitude)` pairs.
pub fn compute_exclusion(
    survey: &ZtfSurvey,
    population: &[(OrbitalElements, f64)],
) -> ExclusionResult {
    let n_total = population.len() as u32;
    let expected_excluded: f64 = population
        .iter()
        .map(|(elements, mag)| detection_probability_for_orbit(survey, elements, *mag))
        .sum();

    let fraction_excluded = if n_total > 0 {
        expected_excluded / n_total as f64
    } else {
        0.0
    };

    ExclusionResult {
        fraction_excluded,
        n_total,
        expected_excluded,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;
    use p9_core::data::reference_population::generate_reference_population;
    use rand::SeedableRng;

    /// H3 regression against the published number: a seeded reference
    /// population from the p9-2021-orbit posterior emulator, piped through
    /// the per-orbit survey model, must reproduce the ZTF exclusion fraction
    /// (0.564, shared table) to within a few points.
    #[test]
    fn paper_exclusion_fraction_from_reference_population() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(2021);
        let population: Vec<(OrbitalElements, f64)> = generate_reference_population(4000, &mut rng)
            .into_iter()
            .map(|p| {
                (
                    OrbitalElements {
                        a: p.a,
                        e: p.e,
                        i: p.i,
                        omega: p.omega,
                        omega_big: p.omega_big,
                        mean_anomaly: p.mean_anomaly,
                    },
                    p.v_magnitude,
                )
            })
            .collect();

        let survey = ZtfSurvey::default();
        let result = compute_exclusion(&survey, &population);

        let published = EXCLUSION_FRACTIONS[0].unique_fraction; // ZTF: 0.564
                                                                // Computed 0.528 with the BB21-faithful reference population
                                                                // ((m/3)R⊕ radius, U(0.2, 0.75) per-object albedo) and the
                                                                // object-level P(>= 7 detections) completeness whose 95% point is
                                                                // pinned at V = 20.5 as BB21 define it (was 0.477 when 20.5 was
                                                                // misread as a 50% logistic midpoint). The remaining -0.036 is the
                                                                // documented residual: per-field cadence/weather structure we don't
                                                                // model.
        assert!(
            (result.fraction_excluded - published).abs() < 0.05,
            "exclusion {:.3} vs published {published}",
            result.fraction_excluded
        );
    }

    /// The footprint must matter: removing the declination cut changes the
    /// exclusion (the old model was invariant to it).
    #[test]
    fn footprint_affects_exclusion() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let population: Vec<(OrbitalElements, f64)> = generate_reference_population(2000, &mut rng)
            .into_iter()
            .map(|p| {
                (
                    OrbitalElements {
                        a: p.a,
                        e: p.e,
                        i: p.i,
                        omega: p.omega,
                        omega_big: p.omega_big,
                        mean_anomaly: p.mean_anomaly,
                    },
                    p.v_magnitude,
                )
            })
            .collect();

        let ztf = ZtfSurvey::default();
        let all_sky = ZtfSurvey {
            dec_limit_deg: -90.0,
            ..ZtfSurvey::default()
        };
        let f_ztf = compute_exclusion(&ztf, &population).fraction_excluded;
        let f_all = compute_exclusion(&all_sky, &population).fraction_excluded;
        assert!(
            f_all > f_ztf + 0.02,
            "all-sky {f_all:.3} should exceed dec-limited {f_ztf:.3}"
        );
    }

    #[test]
    fn all_faint_nothing_excluded() {
        let survey = ZtfSurvey::default();
        let population: Vec<(OrbitalElements, f64)> = (0..200)
            .map(|k| {
                (
                    OrbitalElements {
                        a: 400.0,
                        e: 0.2,
                        i: 15.0 * DEG2RAD,
                        omega: 0.0,
                        omega_big: 0.0,
                        mean_anomaly: k as f64 * DEG2RAD,
                    },
                    24.0,
                )
            })
            .collect();
        let result = compute_exclusion(&survey, &population);
        assert!(result.fraction_excluded < 1e-4);
    }

    /// Combination across surveys uses the shared p9-core table.
    #[test]
    fn combined_exclusion_from_shared_table() {
        let fracs: Vec<f64> = EXCLUSION_FRACTIONS
            .iter()
            .map(|s| s.unique_fraction)
            .collect();
        let combined = combined_unique_exclusion(&fracs);
        assert!((combined - 0.785).abs() < 1e-12);
    }
}
