//! WISE/NEOWISE survey characteristics for the Planet Nine search of
//! Meisner et al. (2018).
//!
//! WISE is a near-all-sky thermal-IR survey. The shift-and-stack co-adds used
//! by Meisner et al. reach a W1 single-frame-stack depth of ≈ 16.5 (Vega),
//! taken from the shared `p9_core::analysis::surveys` table ("WISE W1"). The
//! search covers essentially the whole sky in inertial coordinates; the
//! dominant losses are confusion near the galactic plane (high source density
//! and cirrus) which is masked. We model the footprint as all-sky minus a
//! galactic-latitude exclusion band, evaluated with `p9_core::coords::sky`.

use serde::{Deserialize, Serialize};

use p9_core::analysis::surveys::limiting_magnitude;

/// Published WISE W1 single-exposure / shift-and-stack depth (Vega mag),
/// kept as a labelled reference constant; the live value comes from the
/// shared survey table and the two are pinned equal in the tests.
pub const W1_DEPTH_REFERENCE: f64 = 16.5;

/// Galactic-latitude half-width (degrees) of the masked confusion band about
/// the galactic plane. Meisner et al. (2018) note the galactic plane is the
/// region of poorest sensitivity in the WISE shift-and-stack search; |b| < 10°
/// is a representative mask that removes ~17% of the sky (sin 10° ≈ 0.174).
pub const GALACTIC_MASK_DEG: f64 = 10.0;

/// WISE survey model for the Planet Nine W1 search.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct WiseSurvey {
    /// W1 (3.4 µm) limiting magnitude (Vega), from the shared survey table.
    pub w1_depth: f64,
    /// Logistic steepness of the magnitude-efficiency roll-off (mag⁻¹).
    pub efficiency_steepness: f64,
    /// Half-width of the masked galactic-plane band in galactic latitude
    /// (degrees); positions with |b| below this are excluded from the search.
    pub galactic_mask_deg: f64,
}

impl Default for WiseSurvey {
    fn default() -> Self {
        Self {
            w1_depth: limiting_magnitude("WISE W1").expect("WISE W1 in shared survey table"),
            efficiency_steepness: 4.0,
            galactic_mask_deg: GALACTIC_MASK_DEG,
        }
    }
}

impl WiseSurvey {
    /// Magnitude-dependent recovery efficiency: a logistic roll-off centered
    /// on the W1 depth (ε = 0.5 at the limit, → 1 well above it in flux /
    /// below it in magnitude, → 0 for faint sources).
    pub fn detection_efficiency(&self, w1_mag: f64) -> f64 {
        1.0 / (1.0 + ((w1_mag - self.w1_depth) * self.efficiency_steepness).exp())
    }

    /// Whether a galactic latitude (degrees) lies outside the masked plane.
    pub fn in_footprint(&self, galactic_lat_deg: f64) -> bool {
        galactic_lat_deg.abs() >= self.galactic_mask_deg
    }

    /// Fraction of the celestial sphere covered: all-sky minus the masked
    /// galactic band |b| < galactic_mask. The unmasked solid-angle fraction
    /// is cos-weighted: f = 1 − sin(b_mask).
    pub fn sky_coverage_fraction(&self) -> f64 {
        1.0 - self.galactic_mask_deg.to_radians().sin()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_w1_depth_matches_shared_table() {
        let survey = WiseSurvey::default();
        assert!((survey.w1_depth - 16.5).abs() < 1e-10);
        assert!((survey.w1_depth - W1_DEPTH_REFERENCE).abs() < 1e-10);
    }

    #[test]
    fn efficiency_is_logistic() {
        let survey = WiseSurvey::default();
        assert!(survey.detection_efficiency(14.0) > 0.99);
        assert!((survey.detection_efficiency(16.5) - 0.5).abs() < 1e-10);
        assert!(survey.detection_efficiency(19.0) < 0.01);
    }

    #[test]
    fn galactic_plane_is_masked() {
        let survey = WiseSurvey::default();
        assert!(survey.in_footprint(60.0));
        assert!(survey.in_footprint(10.0));
        assert!(!survey.in_footprint(0.0));
        assert!(!survey.in_footprint(-5.0));
    }

    #[test]
    fn coverage_is_near_all_sky() {
        let survey = WiseSurvey::default();
        let f = survey.sky_coverage_fraction();
        // |b| < 10° removes ~17% of the sky; coverage ~83%.
        assert!(f > 0.80 && f < 0.85, "coverage = {f:.3}");
    }

    #[test]
    fn wider_mask_lowers_coverage() {
        let narrow = WiseSurvey::default();
        let wide = WiseSurvey {
            galactic_mask_deg: 20.0,
            ..WiseSurvey::default()
        };
        assert!(wide.sky_coverage_fraction() < narrow.sky_coverage_fraction());
    }
}
