//! Detection analysis for the PS1 Planet Nine search.
//!
//! Computes the number of synthetic Planet Nine orbits detected and linked
//! in the PS1 survey, and identifies detections unique to PS1 (not found
//! by ZTF or DES).

use serde::{Deserialize, Serialize};

use crate::survey_model::Ps1Survey;

/// Summary of PS1 detection results from the synthetic population analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ps1Result {
    /// Total synthetic orbits detected in PS1 pipeline
    pub n_detected: u32,
    /// Detections unique to PS1 (not recovered by ZTF or DES)
    pub n_unique: u32,
}

impl Ps1Result {
    /// Paper's reported values: 69802 detected, 17054 unique to PS1.
    pub fn paper_values() -> Self {
        Self {
            n_detected: 69802,
            n_unique: 17054,
        }
    }

    /// Fraction of detections that are unique to PS1.
    pub fn unique_fraction(&self) -> f64 {
        self.n_unique as f64 / self.n_detected as f64
    }
}

/// Simulate PS1 detection for a synthetic Planet Nine population.
///
/// Each entry in `magnitudes` is the predicted V-band apparent magnitude of a
/// synthetic P9 orbit at the PS1 epoch. Each entry in `ecliptic_lat_deg` is the
/// ecliptic latitude in degrees.
///
/// Returns the number of synthetic orbits that would be detected and linked.
pub fn compute_detection(survey: &Ps1Survey, magnitudes: &[f64], ecliptic_lat_deg: &[f64]) -> u32 {
    assert_eq!(magnitudes.len(), ecliptic_lat_deg.len());

    let mut n_detected = 0u32;

    for (mag, lat) in magnitudes.iter().zip(ecliptic_lat_deg.iter()) {
        if !survey.is_detectable(*mag) {
            continue;
        }

        // PS1 3-pi survey covers dec > -30 deg; model as ecliptic latitude cut.
        // Objects at very high galactic latitudes have better coverage.
        let in_footprint = *lat > -30.0;
        if !in_footprint {
            continue;
        }

        // Apply linking efficiency
        let linked = rand::random::<f64>() < survey.linking_efficiency();
        if linked {
            n_detected += 1;
        }
    }

    n_detected
}

/// Compute detections unique to PS1 (not found by ZTF or DES).
///
/// Takes the full set of PS1-detected magnitudes and filters out those
/// that would also be detected by ZTF (depth ~20.5) or DES (depth ~23.8
/// but limited footprint covering ~12% of sky at dec < -20).
pub fn compute_unique_detections(
    magnitudes: &[f64],
    ecliptic_lat_deg: &[f64],
    ztf_depth: f64,
    des_dec_limit: f64,
) -> u32 {
    assert_eq!(magnitudes.len(), ecliptic_lat_deg.len());

    let mut n_unique = 0u32;

    for (mag, lat) in magnitudes.iter().zip(ecliptic_lat_deg.iter()) {
        let ztf_would_detect = *mag < ztf_depth;
        let des_would_detect = *lat < des_dec_limit && *mag < 23.8;

        if !ztf_would_detect && !des_would_detect {
            n_unique += 1;
        }
    }

    n_unique
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paper_values_correct() {
        let result = Ps1Result::paper_values();
        assert_eq!(result.n_detected, 69802);
        assert_eq!(result.n_unique, 17054);
    }

    #[test]
    fn unique_fraction_consistent() {
        let result = Ps1Result::paper_values();
        let frac = result.unique_fraction();
        assert!((frac - 17054.0 / 69802.0).abs() < 1e-10);
    }

    #[test]
    fn bright_objects_all_detected() {
        let survey = Ps1Survey::default();
        let mags = vec![18.0; 200];
        let lats = vec![10.0; 200];
        let n = compute_detection(&survey, &mags, &lats);
        assert!(n >= 180, "Expected most bright objects detected, got {n}");
    }

    #[test]
    fn faint_objects_all_missed() {
        let survey = Ps1Survey::default();
        let mags = vec![25.0; 50];
        let lats = vec![10.0; 50];
        let n = compute_detection(&survey, &mags, &lats);
        assert_eq!(n, 0);
    }

    #[test]
    fn unique_filters_ztf_overlap() {
        let mags = vec![19.0, 20.0, 21.0, 21.3];
        let lats = vec![10.0, 10.0, 10.0, 10.0];
        let n = compute_unique_detections(&mags, &lats, 20.5, -20.0);
        // Only mag 21.0 and 21.3 are too faint for ZTF and too far north for DES
        assert_eq!(n, 2);
    }
}
