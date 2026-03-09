//! DES survey characteristics for the Planet Nine search.
//!
//! Models the Dark Energy Survey's footprint, depth, and detection
//! completeness as used in Belyakov, Bernardinelli & Brown (2022).

use serde::{Deserialize, Serialize};

/// Photometric band used in DES observations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum DesBand {
    G,
    R,
    I,
    Z,
    Y,
}

/// DES survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DesSurvey {
    /// Effective footprint area in square degrees
    pub footprint_area: f64,
    /// Total number of observing nights across 6 years
    pub n_nights: u32,
    /// Single-epoch 10-sigma depth in g-band (mag)
    pub depth_g: f64,
    /// Bands used in the survey
    pub bands: Vec<DesBand>,
}

impl Default for DesSurvey {
    fn default() -> Self {
        Self {
            footprint_area: 5000.0,
            n_nights: 575,
            depth_g: 24.1,
            bands: vec![DesBand::G, DesBand::R, DesBand::I, DesBand::Z, DesBand::Y],
        }
    }
}

impl DesSurvey {
    /// Detection completeness as a function of apparent magnitude.
    ///
    /// Uses a logistic (logit) function:
    ///   p(m) = c / (1 + exp(k * (m - m50)))
    ///
    /// where c is the asymptotic completeness, k is the steepness, and
    /// m50 is the magnitude at 50% completeness. Parameters are calibrated
    /// to match DES transient detection pipeline performance.
    pub fn completeness(&self, apparent_mag: f64) -> f64 {
        let c = 0.95;
        let k = 3.0;
        let m50 = self.depth_g - 0.5;
        c / (1.0 + (k * (apparent_mag - m50)).exp())
    }

    /// Simplified DES footprint check.
    ///
    /// DES covers approximately 5000 deg^2 of the southern sky.
    /// This simplified model requires dec < 0 and applies a rough
    /// RA cut to approximate the actual footprint boundaries.
    pub fn is_in_footprint(&self, ra_deg: f64, dec_deg: f64) -> bool {
        if dec_deg > 0.0 {
            return false;
        }
        if dec_deg < -65.0 {
            return false;
        }
        // DES wide-field footprint spans roughly RA 0-120 and 300-360
        // with some additional coverage near the Galactic caps
        (ra_deg < 120.0) || (ra_deg > 300.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_paper_values() {
        let survey = DesSurvey::default();
        assert!((survey.footprint_area - 5000.0).abs() < 1e-10);
        assert_eq!(survey.n_nights, 575);
        assert!((survey.depth_g - 24.1).abs() < 1e-10);
        assert_eq!(survey.bands.len(), 5);
    }

    #[test]
    fn bright_objects_high_completeness() {
        let survey = DesSurvey::default();
        let p = survey.completeness(20.0);
        assert!(
            p > 0.90,
            "Bright objects should have >90% completeness, got {p}"
        );
    }

    #[test]
    fn faint_objects_low_completeness() {
        let survey = DesSurvey::default();
        let p = survey.completeness(26.0);
        assert!(
            p < 0.10,
            "Faint objects should have <10% completeness, got {p}"
        );
    }

    #[test]
    fn completeness_monotonically_decreasing() {
        let survey = DesSurvey::default();
        let mut prev = survey.completeness(18.0);
        for mag_tenth in 185..=270 {
            let mag = mag_tenth as f64 / 10.0;
            let p = survey.completeness(mag);
            assert!(
                p <= prev + 1e-12,
                "Completeness should decrease with magnitude"
            );
            prev = p;
        }
    }

    #[test]
    fn footprint_rejects_northern_sky() {
        let survey = DesSurvey::default();
        assert!(!survey.is_in_footprint(30.0, 10.0));
        assert!(!survey.is_in_footprint(30.0, 45.0));
    }

    #[test]
    fn footprint_accepts_southern_sky() {
        let survey = DesSurvey::default();
        assert!(survey.is_in_footprint(30.0, -30.0));
        assert!(survey.is_in_footprint(350.0, -20.0));
    }

    #[test]
    fn footprint_rejects_deep_south() {
        let survey = DesSurvey::default();
        assert!(!survey.is_in_footprint(30.0, -70.0));
    }
}
