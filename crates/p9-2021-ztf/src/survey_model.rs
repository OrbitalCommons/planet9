//! ZTF survey characteristics for the Planet Nine search.
//!
//! Models the Zwicky Transient Facility's footprint (a *declination-limited*
//! survey from Palomar, δ ≳ −30°), depth, and tracklet-linking efficiency as
//! used in Brown & Batygin (2021). The depth limit comes from the shared
//! survey table in `p9_core::analysis::surveys`.

use p9_core::analysis::surveys::{limiting_magnitude, NORTHERN_SURVEY_DEC_LIMIT_DEG};
use p9_core::units::{degrees, Angle};
use serde::{Deserialize, Serialize};

/// ZTF survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ZtfSurvey {
    /// V-band 95% object-completeness depth: BB21 define V = 20.5 as the
    /// magnitude to which the survey recovers >= [`Self::min_detections`]
    /// epochs for 95% of objects — NOT a 50% logistic midpoint (the earlier
    /// reading here, which under-detected near the limit).
    pub depth_limit: f64,
    /// Single-epoch 50%-efficiency magnitude (V): the DES paper's
    /// independent reading of ZTF, m50 = 20.75 g ~ 20.48 V. Per-epoch
    /// efficiency is a logistic centered here.
    pub single_epoch_m50: f64,
    /// Logistic steepness of the per-epoch efficiency roll-off (mag⁻¹);
    /// calibrated against known-asteroid recoveries in the paper, which show
    /// a sharp (few-tenths of a magnitude) transition.
    pub efficiency_steepness: f64,
    /// Detections required for the tracklet linker to consider an object
    /// found (BB21 require >= 7).
    pub min_detections: u32,
    /// Effective number of independent epochs a sky position gets over the
    /// survey — the one calibrated quantity, chosen as the smallest N whose
    /// object-level completeness at [`Self::depth_limit`] reaches the
    /// defined 95% (N = 22 gives 96.0%).
    pub effective_epochs: u32,
    /// Southern declination limit of the footprint (degrees); ZTF observes
    /// from Palomar and covers δ > −30°.
    pub dec_limit_deg: f64,
}

impl Default for ZtfSurvey {
    fn default() -> Self {
        Self {
            depth_limit: limiting_magnitude("ZTF").expect("ZTF in shared survey table"),
            single_epoch_m50: 20.48,
            efficiency_steepness: 4.0,
            min_detections: 7,
            effective_epochs: 22,
            dec_limit_deg: NORTHERN_SURVEY_DEC_LIMIT_DEG,
        }
    }
}

impl ZtfSurvey {
    /// Single-epoch detection efficiency: a logistic roll-off centered on
    /// [`Self::single_epoch_m50`] (the self-calibration curve from
    /// known-asteroid injections).
    pub fn single_epoch_efficiency(&self, apparent_magnitude: f64) -> f64 {
        1.0 / (1.0
            + ((apparent_magnitude - self.single_epoch_m50) * self.efficiency_steepness).exp())
    }

    /// Object-level detection efficiency: the probability of accumulating at
    /// least [`Self::min_detections`] single-epoch detections over
    /// [`Self::effective_epochs`] independent epochs — the quantity BB21's
    /// V = 20.5 / 95% completeness statement is about. Binomial survival:
    /// P(X >= 7), X ~ Bin(N, ε(m)).
    pub fn detection_efficiency(&self, apparent_magnitude: f64) -> f64 {
        let eps = self.single_epoch_efficiency(apparent_magnitude);
        let n = self.effective_epochs;
        let mut p_lt = 0.0_f64; // P(X < min_detections)
        let mut coeff = 1.0_f64; // C(n, k)
        for k in 0..self.min_detections {
            if k > 0 {
                coeff *= (n - k + 1) as f64 / k as f64;
            }
            p_lt += coeff * eps.powi(k as i32) * (1.0 - eps).powi((n - k) as i32);
        }
        (1.0 - p_lt).clamp(0.0, 1.0)
    }

    /// Whether a sky position falls in the ZTF footprint (declination cut).
    pub fn in_footprint(&self, dec_deg: f64) -> bool {
        dec_deg > self.dec_limit_deg
    }

    /// Southern declination limit of the footprint as a dimension-checked
    /// [`Angle`].
    pub fn dec_limit(&self) -> Angle {
        degrees(self.dec_limit_deg)
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
    use approx::assert_relative_eq;
    use uom::si::angle::degree;

    #[test]
    fn default_matches_shared_survey_table() {
        let survey = ZtfSurvey::default();
        assert!((survey.depth_limit - 20.5).abs() < 1e-10);
        assert!((survey.dec_limit_deg - (-30.0)).abs() < 1e-10);
    }

    #[test]
    fn typed_dec_limit_matches_f64() {
        let survey = ZtfSurvey::default();
        assert_relative_eq!(
            survey.dec_limit().get::<degree>(),
            survey.dec_limit_deg,
            epsilon = 1e-12
        );
    }

    #[test]
    fn object_completeness_matches_bb21_definition() {
        let survey = ZtfSurvey::default();
        // The defining property: 95% object-level completeness (>= 7
        // detections) at the published V = 20.5 depth.
        let p95 = survey.detection_efficiency(survey.depth_limit);
        assert!((0.95..0.98).contains(&p95), "P(20.5) = {p95:.4}");
        // Saturated bright, dead faint, sharp few-tenths-mag transition
        // (the BB21 self-calibration behavior).
        assert!(survey.detection_efficiency(18.0) > 0.999);
        assert!(survey.detection_efficiency(21.0) < 0.05);
        assert!(survey.detection_efficiency(23.0) < 1e-6);
        // Monotone non-increasing with magnitude.
        let mut last = 1.0;
        for k in 0..=60 {
            let m = 18.0 + 0.1 * k as f64;
            let e = survey.detection_efficiency(m);
            assert!(e <= last + 1e-12, "non-monotone at m = {m}");
            last = e;
        }
        // Per-epoch curve still sits where the DES reading puts it.
        let e50 = survey.single_epoch_efficiency(survey.single_epoch_m50);
        assert!((e50 - 0.5).abs() < 1e-10);
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
