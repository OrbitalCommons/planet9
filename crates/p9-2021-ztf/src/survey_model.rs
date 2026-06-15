//! ZTF survey characteristics for the Planet Nine search.
//!
//! Models the Zwicky Transient Facility's footprint (a *declination-limited*
//! survey from Palomar, δ ≳ −30°), depth, and tracklet-linking efficiency as
//! used in Brown & Batygin (2021). The depth limit comes from the shared
//! survey table in `p9_core::analysis::surveys`.

use p9_core::analysis::surveys::limiting_magnitude;
use serde::{Deserialize, Serialize};

/// ZTF survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZtfSurvey {
    /// 50%-efficiency magnitude of the recovery pipeline (~20.5, ZTF r-band)
    pub depth_limit: f64,
    /// Logistic steepness of the efficiency roll-off near the depth limit
    /// (mag⁻¹); calibrated against known-asteroid recoveries in the paper,
    /// which show a sharp (few-tenths of a magnitude) transition.
    pub efficiency_steepness: f64,
    /// Southern declination limit of the footprint (degrees); ZTF observes
    /// from Palomar and covers δ > −30°.
    pub dec_limit_deg: f64,
}

impl Default for ZtfSurvey {
    fn default() -> Self {
        Self {
            depth_limit: limiting_magnitude("ZTF").expect("ZTF in shared survey table"),
            efficiency_steepness: 4.0,
            dec_limit_deg: -30.0,
        }
    }
}

impl ZtfSurvey {
    /// Magnitude-dependent single-epoch detection efficiency: a logistic
    /// roll-off centered on the depth limit (the self-calibration curve from
    /// known-asteroid injections), replacing the previous hard step at 20.5.
    pub fn detection_efficiency(&self, apparent_magnitude: f64) -> f64 {
        1.0 / (1.0 + ((apparent_magnitude - self.depth_limit) * self.efficiency_steepness).exp())
    }

    /// Whether a sky position falls in the ZTF footprint (declination cut).
    pub fn in_footprint(&self, dec_deg: f64) -> bool {
        dec_deg > self.dec_limit_deg
    }

    /// Tracklet-linking efficiency measured from known asteroid recoveries.
    ///
    /// Brown & Batygin (2021) report 99.66% linking efficiency for objects
    /// with sufficient detections, calibrated against known solar system objects.
    pub fn linking_efficiency(&self) -> f64 {
        0.9966
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_shared_survey_table() {
        let survey = ZtfSurvey::default();
        assert!((survey.depth_limit - 20.5).abs() < 1e-10);
        assert!((survey.dec_limit_deg - (-30.0)).abs() < 1e-10);
    }

    #[test]
    fn efficiency_is_logistic_not_step() {
        let survey = ZtfSurvey::default();
        // Near unity well above the limit, ~0.5 at the limit, ~0 below.
        assert!(survey.detection_efficiency(18.0) > 0.99);
        assert!((survey.detection_efficiency(20.5) - 0.5).abs() < 1e-10);
        assert!(survey.detection_efficiency(23.0) < 0.01);
        // Smooth: intermediate efficiency just past the limit.
        let e = survey.detection_efficiency(20.8);
        assert!(e > 0.05 && e < 0.5, "ε(20.8) = {e:.3}");
    }

    #[test]
    fn footprint_is_declination_limited() {
        let survey = ZtfSurvey::default();
        assert!(survey.in_footprint(0.0));
        assert!(survey.in_footprint(60.0));
        assert!(!survey.in_footprint(-45.0));
    }

    #[test]
    fn linking_efficiency_matches_paper() {
        let survey = ZtfSurvey::default();
        assert!((survey.linking_efficiency() - 0.9966).abs() < 1e-10);
    }
}
