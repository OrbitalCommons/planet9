//! PS1 survey characteristics for the Planet Nine search.
//!
//! Models the Pan-STARRS1 sky coverage, depth limits, quality cuts,
//! and tracklet-linking requirements as used in Brown, Holman & Batygin (2024).

use serde::{Deserialize, Serialize};

/// Quality-cut thresholds applied to PS1 detections.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityCuts {
    /// Minimum PSF quality factor (psfQfPerfect threshold)
    pub psf_qf_perfect_min: f64,
    /// Maximum apparent magnitude retained after quality filtering
    pub mag_max: f64,
}

/// PS1 survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ps1Survey {
    /// Single-exposure depth limit in V-band magnitude (~21.5 for PS1)
    pub depth_limit: f64,
    /// Fraction of sky covered by the PS1 3-pi survey (~0.75)
    pub coverage_frac: f64,
    /// Minimum number of detections required to link a moving object (n >= 9)
    pub linking_threshold: u32,
    /// Quality cuts applied to individual detections
    pub quality_cuts: QualityCuts,
}

impl Default for Ps1Survey {
    fn default() -> Self {
        Self {
            depth_limit: 21.5,
            coverage_frac: 0.75,
            linking_threshold: 9,
            quality_cuts: QualityCuts {
                psf_qf_perfect_min: 0.99,
                mag_max: 22.5,
            },
        }
    }
}

impl Ps1Survey {
    /// Returns true if an object at the given apparent magnitude is bright enough
    /// to be detected in a single PS1 exposure.
    pub fn is_detectable(&self, apparent_magnitude: f64) -> bool {
        apparent_magnitude < self.depth_limit
    }

    /// Apply quality cuts to a detection. Returns true if the detection passes.
    pub fn quality_cuts(&self, psf_qf_perfect: f64, magnitude: f64) -> bool {
        psf_qf_perfect > self.quality_cuts.psf_qf_perfect_min
            && magnitude < self.quality_cuts.mag_max
    }

    /// Tracklet-linking efficiency measured from known asteroid recoveries.
    ///
    /// Brown, Holman & Batygin (2024) report 99.2% linking efficiency for
    /// objects with sufficient detections in the PS1 cadence.
    pub fn linking_efficiency(&self) -> f64 {
        0.992
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_paper_values() {
        let survey = Ps1Survey::default();
        assert!((survey.depth_limit - 21.5).abs() < 1e-10);
        assert!((survey.coverage_frac - 0.75).abs() < 1e-10);
        assert_eq!(survey.linking_threshold, 9);
    }

    #[test]
    fn quality_cuts_match_paper() {
        let survey = Ps1Survey::default();
        assert!((survey.quality_cuts.psf_qf_perfect_min - 0.99).abs() < 1e-10);
        assert!((survey.quality_cuts.mag_max - 22.5).abs() < 1e-10);
    }

    #[test]
    fn bright_object_is_detectable() {
        let survey = Ps1Survey::default();
        assert!(survey.is_detectable(19.0));
        assert!(survey.is_detectable(21.4));
    }

    #[test]
    fn faint_object_not_detectable() {
        let survey = Ps1Survey::default();
        assert!(!survey.is_detectable(22.0));
        assert!(!survey.is_detectable(25.0));
    }

    #[test]
    fn quality_cuts_filter_correctly() {
        let survey = Ps1Survey::default();
        assert!(survey.quality_cuts(0.995, 21.0));
        assert!(!survey.quality_cuts(0.98, 21.0));
        assert!(!survey.quality_cuts(0.995, 23.0));
        assert!(!survey.quality_cuts(0.5, 24.0));
    }

    #[test]
    fn linking_efficiency_matches_paper() {
        let survey = Ps1Survey::default();
        assert!((survey.linking_efficiency() - 0.992).abs() < 1e-10);
    }
}
