//! Detection efficiency analysis for the ZTF Planet Nine search.
//!
//! The per-orbit detection probability is the product of
//! * the magnitude-dependent logistic recovery efficiency (the
//!   self-calibration curve, wired into `ZtfSurvey::detection_efficiency` —
//!   previously dead physics next to a hard step at 20.5),
//! * the declination footprint (ZTF is declination-limited, δ > −30°), and
//! * the measured tracklet-linking efficiency (99.66%).
//!
//! Expected counts are accumulated *deterministically* (sums of
//! probabilities), so no RNG is involved; the previous version drew unseeded
//! `rand::random` per object. A seeded Monte Carlo realization is available
//! for callers that need integer counts.

use rand::Rng;
use serde::{Deserialize, Serialize};

use p9_core::types::OrbitalElements;

use crate::sky::declination_deg;
use crate::survey_model::ZtfSurvey;

/// Result of a detection-efficiency computation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectionResult {
    /// Number of synthetic P9 orbits injected
    pub n_injected: u32,
    /// Expected number passing the magnitude efficiency alone
    pub expected_detected: f64,
    /// Expected number detected, in-footprint, and linked
    pub expected_linked: f64,
    /// Overall expected detection efficiency (expected_linked / n_injected)
    pub efficiency: f64,
}

/// Per-orbit detection probability: ε(m) × footprint(δ) × linking.
pub fn detection_probability(survey: &ZtfSurvey, magnitude: f64, dec_deg: f64) -> f64 {
    let eps = survey.detection_efficiency(magnitude);
    let footprint = if survey.in_footprint(dec_deg) {
        1.0
    } else {
        0.0
    };
    eps * footprint * survey.linking_efficiency()
}

/// Detection probability for an orbit, computing the sky position (and hence
/// the declination cut) from the orbital elements.
pub fn detection_probability_for_orbit(
    survey: &ZtfSurvey,
    elements: &OrbitalElements,
    magnitude: f64,
) -> f64 {
    detection_probability(survey, magnitude, declination_deg(elements))
}

/// Deterministic expected detection efficiency for a set of synthetic
/// Planet Nine orbits with predicted magnitudes.
pub fn compute_detection_efficiency(
    survey: &ZtfSurvey,
    orbits: &[(OrbitalElements, f64)],
) -> DetectionResult {
    let n_injected = orbits.len() as u32;
    let mut expected_detected = 0.0;
    let mut expected_linked = 0.0;

    for (elements, mag) in orbits {
        expected_detected += survey.detection_efficiency(*mag);
        expected_linked += detection_probability_for_orbit(survey, elements, *mag);
    }

    let efficiency = if n_injected > 0 {
        expected_linked / n_injected as f64
    } else {
        0.0
    };

    DetectionResult {
        n_injected,
        expected_detected,
        expected_linked,
        efficiency,
    }
}

/// Seeded Monte Carlo realization of the detection counts (one Bernoulli per
/// orbit at its detection probability). Thread your own seeded `Rng`.
pub fn realize_detections<R: Rng>(
    survey: &ZtfSurvey,
    orbits: &[(OrbitalElements, f64)],
    rng: &mut R,
) -> u32 {
    orbits
        .iter()
        .filter(|(elements, mag)| {
            rng.gen::<f64>() < detection_probability_for_orbit(survey, elements, *mag)
        })
        .count() as u32
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;
    use rand::SeedableRng;

    fn orbit(i_deg: f64, m_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: 400.0,
            e: 0.2,
            i: i_deg * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: m_deg * DEG2RAD,
        }
    }

    #[test]
    fn bright_northern_orbits_detected() {
        let survey = ZtfSurvey::default();
        let orbits: Vec<(OrbitalElements, f64)> = (0..50)
            .map(|k| (orbit(5.0, k as f64 * 2.0), 18.0))
            .collect();
        let result = compute_detection_efficiency(&survey, &orbits);
        assert_eq!(result.n_injected, 50);
        // Low-inclination orbits never reach δ < −30°, magnitudes are bright.
        assert!(
            result.efficiency > 0.9,
            "efficiency = {}",
            result.efficiency
        );
    }

    #[test]
    fn faint_objects_missed() {
        let survey = ZtfSurvey::default();
        let orbits: Vec<(OrbitalElements, f64)> = (0..50)
            .map(|k| (orbit(5.0, k as f64 * 2.0), 23.0))
            .collect();
        let result = compute_detection_efficiency(&survey, &orbits);
        assert!(
            result.efficiency < 0.01,
            "efficiency = {}",
            result.efficiency
        );
    }

    #[test]
    fn southern_orbits_outside_footprint() {
        let survey = ZtfSurvey::default();
        // High-inclination orbit positioned far south of the ecliptic.
        let mut south = None;
        for m in 0..360 {
            let o = orbit(40.0, m as f64);
            if declination_deg(&o) < -35.0 {
                south = Some(o);
                break;
            }
        }
        let south = south.expect("found a southern sky position");
        assert_eq!(detection_probability_for_orbit(&survey, &south, 18.0), 0.0);
    }

    #[test]
    fn deterministic_and_mc_agree() {
        let survey = ZtfSurvey::default();
        let orbits: Vec<(OrbitalElements, f64)> = (0..400)
            .map(|k| (orbit(15.0, k as f64 * 0.9), 19.5 + (k % 30) as f64 * 0.1))
            .collect();
        let expected = compute_detection_efficiency(&survey, &orbits).expected_linked;

        let mut rng = rand::rngs::StdRng::seed_from_u64(99);
        let n_runs = 200;
        let mean_mc: f64 = (0..n_runs)
            .map(|_| realize_detections(&survey, &orbits, &mut rng) as f64)
            .sum::<f64>()
            / n_runs as f64;
        assert!(
            (mean_mc - expected).abs() < 0.05 * expected.max(1.0) + 3.0,
            "MC mean {mean_mc:.1} vs expected {expected:.1}"
        );
    }
}
