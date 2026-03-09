//! ZTF survey characteristics for the Planet Nine search.
//!
//! Models the Zwicky Transient Facility's sky coverage, depth limits,
//! and tracklet-linking requirements as used in Brown & Batygin (2021).

use serde::{Deserialize, Serialize};

/// ZTF survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZtfSurvey {
    /// Single-exposure depth limit in V-band magnitude (~20.5 for ZTF r-band)
    pub depth_limit: f64,
    /// Fraction of sky covered by ZTF public survey (~0.75 of accessible sky)
    pub sky_coverage_frac: f64,
    /// Minimum number of detections required to form a tracklet and link
    pub linking_threshold: u32,
}

impl Default for ZtfSurvey {
    fn default() -> Self {
        Self {
            depth_limit: 20.5,
            sky_coverage_frac: 0.75,
            linking_threshold: 7,
        }
    }
}

impl ZtfSurvey {
    /// Returns true if an object at the given apparent magnitude is bright enough
    /// to be detected in a single ZTF exposure.
    pub fn is_detectable(&self, apparent_magnitude: f64) -> bool {
        apparent_magnitude < self.depth_limit
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
    fn default_matches_paper_values() {
        let survey = ZtfSurvey::default();
        assert!((survey.depth_limit - 20.5).abs() < 1e-10);
        assert!((survey.sky_coverage_frac - 0.75).abs() < 1e-10);
        assert_eq!(survey.linking_threshold, 7);
    }

    #[test]
    fn bright_object_is_detectable() {
        let survey = ZtfSurvey::default();
        assert!(survey.is_detectable(19.0));
        assert!(survey.is_detectable(20.4));
    }

    #[test]
    fn faint_object_not_detectable() {
        let survey = ZtfSurvey::default();
        assert!(!survey.is_detectable(21.0));
        assert!(!survey.is_detectable(25.0));
    }

    #[test]
    fn linking_efficiency_matches_paper() {
        let survey = ZtfSurvey::default();
        assert!((survey.linking_efficiency() - 0.9966).abs() < 1e-10);
    }
}
