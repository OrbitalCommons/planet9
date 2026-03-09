//! Detection efficiency analysis for the ZTF Planet Nine search.
//!
//! Computes the fraction of synthetic Planet Nine orbits that would be
//! detected and linked in the ZTF survey, accounting for sky coverage,
//! depth limits, and tracklet-linking requirements.

use serde::{Deserialize, Serialize};

use crate::survey_model::ZtfSurvey;

/// Result of a detection-efficiency Monte Carlo run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DetectionResult {
    /// Number of synthetic P9 orbits injected
    pub n_injected: u32,
    /// Number passing the single-exposure depth cut
    pub n_detected: u32,
    /// Number successfully linked into tracklets
    pub n_linked: u32,
    /// Overall detection efficiency (n_linked / n_injected)
    pub efficiency: f64,
}

/// Compute detection efficiency from a set of synthetic Planet Nine
/// apparent magnitudes and ecliptic latitudes.
///
/// Each entry in `magnitudes` is the predicted V-band magnitude of a
/// synthetic P9 orbit at the survey epoch. Each entry in `ecliptic_lat_deg`
/// is the ecliptic latitude in degrees (used for sky-coverage cuts).
///
/// Returns the fraction that would be both detected and linked.
pub fn compute_detection_efficiency(
    survey: &ZtfSurvey,
    magnitudes: &[f64],
    ecliptic_lat_deg: &[f64],
) -> DetectionResult {
    assert_eq!(magnitudes.len(), ecliptic_lat_deg.len());

    let n_injected = magnitudes.len() as u32;
    let mut n_detected = 0u32;
    let mut n_linked = 0u32;

    for (mag, lat) in magnitudes.iter().zip(ecliptic_lat_deg.iter()) {
        if !survey.is_detectable(*mag) {
            continue;
        }
        n_detected += 1;

        // ZTF covers |b_ecl| < ~60 deg well; galactic plane avoidance
        // removes a band around |b_gal| < 7 deg. We model this as a
        // flat sky-coverage probability applied to each injected orbit.
        let in_footprint = lat.abs() < 60.0;
        if !in_footprint {
            continue;
        }

        // Apply linking efficiency
        let linked = rand::random::<f64>() < survey.linking_efficiency();
        if linked {
            n_linked += 1;
        }
    }

    let efficiency = if n_injected > 0 {
        n_linked as f64 / n_injected as f64
    } else {
        0.0
    };

    DetectionResult {
        n_injected,
        n_detected,
        n_linked,
        efficiency,
    }
}

/// Self-calibration correction factor derived from known asteroid recoveries.
///
/// Brown & Batygin (2021) inject known asteroids into the pipeline and
/// measure the recovery rate as a function of magnitude. The correction
/// factor accounts for systematic losses not captured by the simple
/// depth-limit model (e.g., trailing losses, crowded fields).
///
/// Returns the multiplicative correction to apply to the raw efficiency.
pub fn self_calibration_correction(magnitude: f64) -> f64 {
    // Near the survey limit, trailing losses and confusion reduce recovery.
    // Model as a smooth sigmoid roll-off centered at V = 20.0.
    let midpoint = 20.0;
    let steepness = 2.0;
    1.0 / (1.0 + ((magnitude - midpoint) * steepness).exp())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_bright_objects_detected() {
        let survey = ZtfSurvey::default();
        let mags = vec![18.0; 100];
        let lats = vec![10.0; 100];
        let result = compute_detection_efficiency(&survey, &mags, &lats);
        assert_eq!(result.n_injected, 100);
        assert_eq!(result.n_detected, 100);
        // With 99.66% linking, nearly all should link
        assert!(result.n_linked >= 90);
    }

    #[test]
    fn all_faint_objects_missed() {
        let survey = ZtfSurvey::default();
        let mags = vec![23.0; 50];
        let lats = vec![10.0; 50];
        let result = compute_detection_efficiency(&survey, &mags, &lats);
        assert_eq!(result.n_injected, 50);
        assert_eq!(result.n_detected, 0);
        assert_eq!(result.n_linked, 0);
        assert!((result.efficiency - 0.0).abs() < 1e-10);
    }

    #[test]
    fn calibration_bright_near_unity() {
        let correction = self_calibration_correction(17.0);
        assert!(
            correction > 0.99,
            "Bright objects: correction = {correction}"
        );
    }

    #[test]
    fn calibration_faint_near_zero() {
        let correction = self_calibration_correction(24.0);
        assert!(
            correction < 0.01,
            "Faint objects: correction = {correction}"
        );
    }
}
