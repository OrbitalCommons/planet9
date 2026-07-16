//! Parameter-space exclusion for the Holman et al. (2025) PS1 search.
//!
//! Mirrors the ZTF/panstarrs exclusion pattern: a seeded reference P9
//! population (the `p9-2021-orbit` posterior emulator) is piped through the
//! shift-and-stack survey model. Each orbit's declination places it in or
//! out of the 3-pi footprint, and its apparent r magnitude sets its
//! recovery efficiency. Since Holman et al. find no Planet Nine, the
//! expected exclusion fraction is the mean per-orbit detection probability.
//!
//! Independent of and complementary to `p9-2024-panstarrs`: that crate
//! models the catalog-linking search; this one the pixel-level
//! shift-and-stack search of the same PS1 3-pi data.

use rand::{Rng, SeedableRng};
use serde::{Deserialize, Serialize};

use p9_core::coords::sky::apparent_position_deg;
use p9_core::data::reference_population::generate_reference_population;
use p9_core::types::OrbitalElements;
use p9_core::units::{degrees, Angle};

use crate::survey_model::Ps1StackSurvey;
use p9_core::analysis::photometry::SOLAR_V_MINUS_R;

/// A reference-population member reduced to the observables the survey model
/// needs: an equatorial declination and an apparent r magnitude.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct DetectionTarget {
    pub dec_deg: f64,
    pub r_magnitude: f64,
}

impl DetectionTarget {
    /// Equatorial declination as a dimension-checked [`Angle`] (degree storage).
    pub fn declination(&self) -> Angle {
        degrees(self.dec_deg)
    }
}

/// Result of the exclusion analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusionResult {
    /// Expected fraction of the reference population this PS1 search detects
    /// (and so excludes, given the non-detection).
    pub fraction_excluded: f64,
    /// Number of reference draws evaluated.
    pub n_total: u32,
    /// Expected number of draws detected (sum of per-orbit probabilities).
    pub expected_excluded: f64,
}

/// Equatorial declination (deg) of an orbit at its stored epoch, via the
/// shared `p9_core::coords::sky` ecliptic->equatorial apparent-position
/// conversion.
pub fn declination_deg(elements: &OrbitalElements) -> f64 {
    let (_, dec, _) = apparent_position_deg(elements, 0.0);
    dec
}

/// Generate a seeded reference population reduced to detection targets,
/// drawn from the `p9-2021-orbit` posterior emulator.
pub fn reference_targets(n: usize, seed: u64) -> Vec<DetectionTarget> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    generate_reference_population(n, &mut rng)
        .into_iter()
        .map(|p| {
            let elements = OrbitalElements {
                a: p.a,
                e: p.e,
                i: p.i,
                omega: p.omega,
                omega_big: p.omega_big,
                mean_anomaly: p.mean_anomaly,
            };
            DetectionTarget {
                dec_deg: declination_deg(&elements),
                // The reference population carries V; the PS1 depth is r.
                // Convert with the labeled neutral-reflector color.
                r_magnitude: p.v_magnitude - SOLAR_V_MINUS_R,
            }
        })
        .collect()
}

/// Expected exclusion fraction: mean per-target detection probability under
/// the survey model.
pub fn compute_exclusion(survey: &Ps1StackSurvey, targets: &[DetectionTarget]) -> ExclusionResult {
    let n_total = targets.len() as u32;
    let expected_excluded: f64 = targets
        .iter()
        .map(|t| survey.detection_probability(t.r_magnitude, t.dec_deg))
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

/// Seeded Monte-Carlo realization of the detected count (one Bernoulli per
/// target at its detection probability).
pub fn monte_carlo_detected(
    survey: &Ps1StackSurvey,
    targets: &[DetectionTarget],
    seed: u64,
) -> u32 {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    targets
        .iter()
        .filter(|t| rng.gen::<f64>() < survey.detection_probability(t.r_magnitude, t.dec_deg))
        .count() as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;

    #[test]
    fn declination_accessor_matches_f64_field() {
        use approx::assert_relative_eq;
        use uom::si::angle::degree;

        let t = DetectionTarget {
            dec_deg: -12.5,
            r_magnitude: 22.0,
        };
        assert_relative_eq!(t.declination().get::<degree>(), t.dec_deg, epsilon = 1e-12);
    }

    #[test]
    fn exclusion_fraction_in_open_unit_interval() {
        let targets = reference_targets(4000, 2025);
        let survey = Ps1StackSurvey::default();
        let result = compute_exclusion(&survey, &targets);
        // PINNED: a genuine fraction in (0, 1), neither trivially 0 nor 1.
        assert!(
            result.fraction_excluded > 0.0 && result.fraction_excluded < 1.0,
            "fraction = {}",
            result.fraction_excluded
        );
        // Order-of-magnitude anchor: a few tens of percent.
        assert!(
            result.fraction_excluded > 0.1 && result.fraction_excluded < 0.9,
            "fraction {} not in plausible range around {}",
            result.fraction_excluded,
            published::EXCLUDED_FRACTION_ORDER
        );
    }

    #[test]
    fn exclusion_increases_with_assumed_mass() {
        // Re-evaluate the same orbits but assign each the apparent magnitude
        // of a P9 of fixed assumed mass at the orbit's heliocentric distance.
        // More massive => brighter => larger excluded fraction.
        use crate::photometry::apparent_r_magnitude;
        use p9_core::data::reference_population::heliocentric_distance;
        use p9_core::types::P9Params;
        use rand::SeedableRng;

        let mut rng = rand::rngs::StdRng::seed_from_u64(2025);
        let pop = generate_reference_population(4000, &mut rng);

        let survey = Ps1StackSurvey::default();
        let exclusion_at_mass = |mass: f64| -> f64 {
            let targets: Vec<DetectionTarget> = pop
                .iter()
                .map(|p| {
                    let elements = OrbitalElements {
                        a: p.a,
                        e: p.e,
                        i: p.i,
                        omega: p.omega,
                        omega_big: p.omega_big,
                        mean_anomaly: p.mean_anomaly,
                    };
                    let params = P9Params {
                        mass_earth: mass,
                        a: p.a,
                        e: p.e,
                        i: p.i,
                        omega: p.omega,
                        omega_big: p.omega_big,
                        mean_anomaly: p.mean_anomaly,
                    };
                    let r = heliocentric_distance(&params);
                    DetectionTarget {
                        dec_deg: declination_deg(&elements),
                        r_magnitude: apparent_r_magnitude(mass, r),
                    }
                })
                .collect();
            compute_exclusion(&survey, &targets).fraction_excluded
        };

        let f5 = exclusion_at_mass(5.0);
        let f10 = exclusion_at_mass(10.0);
        let f20 = exclusion_at_mass(20.0);
        // PINNED: exclusion fraction increases with assumed mass.
        assert!(f10 > f5, "f10 {f10:.3} should exceed f5 {f5:.3}");
        assert!(f20 > f10, "f20 {f20:.3} should exceed f10 {f10:.3}");
    }

    #[test]
    fn footprint_reduces_exclusion_vs_all_sky() {
        let targets = reference_targets(4000, 7);
        let ps1 = Ps1StackSurvey::default();
        let all_sky = Ps1StackSurvey {
            dec_limit_deg: -90.0,
            ..Ps1StackSurvey::default()
        };
        let f_ps1 = compute_exclusion(&ps1, &targets).fraction_excluded;
        let f_all = compute_exclusion(&all_sky, &targets).fraction_excluded;
        // PINNED: the dec > -30 footprint excludes less than all-sky.
        assert!(
            f_all > f_ps1,
            "all-sky {f_all:.3} should exceed dec-limited {f_ps1:.3}"
        );
    }

    #[test]
    fn monte_carlo_seeded_reproducible() {
        let targets = reference_targets(2000, 11);
        let survey = Ps1StackSurvey::default();
        let a = monte_carlo_detected(&survey, &targets, 99);
        let b = monte_carlo_detected(&survey, &targets, 99);
        assert_eq!(a, b, "same seed must reproduce");
        // And it tracks the expected count.
        let expected = compute_exclusion(&survey, &targets).expected_excluded;
        assert!((a as f64 - expected).abs() < 0.1 * targets.len() as f64);
    }

    #[test]
    fn all_faint_nothing_excluded() {
        let survey = Ps1StackSurvey::default();
        let targets: Vec<DetectionTarget> = (0..200)
            .map(|k| DetectionTarget {
                dec_deg: 10.0 + (k % 30) as f64,
                r_magnitude: 26.0,
            })
            .collect();
        let result = compute_exclusion(&survey, &targets);
        assert!(result.fraction_excluded < 1e-4);
    }
}
