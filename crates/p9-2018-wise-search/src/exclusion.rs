//! Parameter-space exclusion from the WISE W1 non-detection.
//!
//! Meisner et al. (2018) found no Planet Nine in the WISE shift-and-stack
//! search. The fraction of a reference Planet Nine population that WISE
//! *would* have detected is therefore the fraction of parameter space the
//! non-detection excludes. Each realization's per-object detection
//! probability is the product of
//!   * the W1 magnitude-efficiency at its thermal+reflected W1 magnitude
//!     (`thermal_model` + `survey_model`), and
//!   * the galactic-plane footprint mask (`sky`/`survey_model`).
//! The expected exclusion fraction is the mean of these probabilities over a
//! seeded population sampled across the Batygin & Brown mass/distance/orbit
//! ranges.

use serde::{Deserialize, Serialize};

use crate::detectability::detection_probability;
use crate::sky::galactic_latitude_deg;
use crate::survey_model::WiseSurvey;
use crate::thermal_model::P9Thermal;
use p9_core::types::OrbitalElements;
use p9_core::units::{au, earth_masses, Length, Mass};

/// One reference Planet Nine realization: orbital elements, mass, and the
/// heliocentric distance at its current mean anomaly.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ReferenceP9 {
    pub elements: OrbitalElements,
    pub mass_earth: f64,
    pub distance_au: f64,
}

impl ReferenceP9 {
    /// Mass as a typed [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }

    /// Heliocentric distance as a typed [`Length`].
    pub fn distance(&self) -> Length {
        au(self.distance_au)
    }
}

/// Result of the WISE exclusion analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusionResult {
    /// Expected fraction of the reference population WISE would have detected.
    pub fraction_excluded: f64,
    /// Number of realizations evaluated.
    pub n_total: u32,
    /// Expected number detected (sum of per-object probabilities).
    pub expected_excluded: f64,
}

/// Compute the expected WISE W1 exclusion fraction over a reference
/// population, accumulating per-object detection probabilities
/// deterministically (no RNG inside).
pub fn compute_exclusion(survey: &WiseSurvey, population: &[ReferenceP9]) -> ExclusionResult {
    let n_total = population.len() as u32;
    let expected_excluded: f64 = population
        .iter()
        .map(|p| {
            let p9 = P9Thermal::new(p.mass_earth, p.distance_au);
            let b = galactic_latitude_deg(&p.elements);
            detection_probability(survey, &p9, b)
        })
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
    use crate::thermal_model::w1_magnitude;
    use p9_core::data::reference_population::{
        generate_reference_population, heliocentric_distance,
    };
    use p9_core::types::P9Params;
    use rand::SeedableRng;

    /// Build a WISE reference population from the shared p9-2021-orbit
    /// posterior emulator, recomputing heliocentric distance from the sampled
    /// elements (the emulator's `v_magnitude` is reflected-optical and not
    /// used here — WISE works in W1).
    fn wise_population(n: usize, seed: u64) -> Vec<ReferenceP9> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        generate_reference_population(n, &mut rng)
            .into_iter()
            .map(|p| {
                let params = P9Params {
                    mass_earth: p.mass,
                    a: p.a,
                    e: p.e,
                    i: p.i,
                    omega: p.omega,
                    omega_big: p.omega_big,
                    mean_anomaly: p.mean_anomaly,
                };
                ReferenceP9 {
                    elements: OrbitalElements {
                        a: p.a,
                        e: p.e,
                        i: p.i,
                        omega: p.omega,
                        omega_big: p.omega_big,
                        mean_anomaly: p.mean_anomaly,
                    },
                    mass_earth: p.mass,
                    distance_au: heliocentric_distance(&params),
                }
            })
            .collect()
    }

    #[test]
    fn reference_typed_accessors_match_f64() {
        use approx::assert_relative_eq;
        use uom::si::length::astronomical_unit;
        use uom::si::mass::kilogram;
        let p = ReferenceP9 {
            elements: OrbitalElements {
                a: 500.0,
                e: 0.2,
                i: 0.3,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            mass_earth: 6.2,
            distance_au: 480.0,
        };
        assert_relative_eq!(
            p.distance().get::<astronomical_unit>(),
            480.0,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.mass().get::<kilogram>(),
            6.2 * p9_core::units::EARTH_MASS_KG,
            epsilon = 1e12
        );
    }

    #[test]
    fn exclusion_fraction_is_a_probability() {
        let pop = wise_population(3000, 2018);
        let survey = WiseSurvey::default();
        let result = compute_exclusion(&survey, &pop);
        assert!(
            result.fraction_excluded > 0.0 && result.fraction_excluded < 1.0,
            "fraction = {}",
            result.fraction_excluded
        );
    }

    #[test]
    fn exclusion_increases_with_assumed_mass() {
        // Shift the whole population to a heavier (brighter in W1) mass and
        // the exclusion fraction must increase.
        let base = wise_population(3000, 7);
        let survey = WiseSurvey::default();

        let light: Vec<ReferenceP9> = base
            .iter()
            .map(|p| ReferenceP9 {
                mass_earth: 5.0,
                ..*p
            })
            .collect();
        let heavy: Vec<ReferenceP9> = base
            .iter()
            .map(|p| ReferenceP9 {
                mass_earth: 18.0,
                ..*p
            })
            .collect();

        let f_light = compute_exclusion(&survey, &light).fraction_excluded;
        let f_heavy = compute_exclusion(&survey, &heavy).fraction_excluded;
        assert!(
            f_heavy > f_light + 0.005,
            "heavy {f_heavy:.4} should exceed light {f_light:.4}"
        );
    }

    #[test]
    fn footprint_matters() {
        // Removing the galactic mask (all-sky) raises the exclusion fraction.
        let pop = wise_population(3000, 11);
        let masked = WiseSurvey::default();
        let all_sky = WiseSurvey {
            galactic_mask_deg: 0.0,
            ..WiseSurvey::default()
        };
        let f_masked = compute_exclusion(&masked, &pop).fraction_excluded;
        let f_all = compute_exclusion(&all_sky, &pop).fraction_excluded;
        // The galactic mask removes ~17% of the sky, so all-sky exceeds the
        // masked footprint by roughly that factor of the (small) exclusion.
        assert!(
            f_all > f_masked * 1.05,
            "all-sky {f_all:.4} should exceed masked {f_masked:.4}"
        );
    }

    #[test]
    fn nearby_massive_excluded_distant_lowmass_not() {
        let survey = WiseSurvey::default();
        // Nearby + massive: brighter than the W1 depth.
        assert!(w1_magnitude(18.0, 180.0) < survey.w1_depth);
        // Distant + low mass: fainter than the depth (undetectable, not excluded).
        assert!(w1_magnitude(5.0, 1200.0) > survey.w1_depth);
    }

    #[test]
    fn all_distant_lowmass_population_barely_excluded() {
        // A population of cold, light, very distant planets: WISE detects
        // essentially none of them.
        let survey = WiseSurvey::default();
        let pop: Vec<ReferenceP9> = (0..200)
            .map(|k| ReferenceP9 {
                elements: OrbitalElements {
                    a: 1500.0,
                    e: 0.1,
                    i: 0.3,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: k as f64 * 0.03,
                },
                mass_earth: 4.0,
                distance_au: 1400.0,
            })
            .collect();
        let result = compute_exclusion(&survey, &pop);
        assert!(
            result.fraction_excluded < 0.01,
            "f = {}",
            result.fraction_excluded
        );
    }
}
