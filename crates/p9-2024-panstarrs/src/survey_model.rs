//! PS1 survey characteristics for the Planet Nine search of
//! Brown, Holman & Batygin (2024), AJ (arXiv:2401.17977).
//!
//! Models the Pan-STARRS1 3-pi coverage (all sky north of dec = -30 deg),
//! the multi-epoch cadence (the DR2 catalog provides ~12 visits per filter
//! over 2009-2015), the >= 9-detection linking requirement, the quality
//! cuts, and the measured 99.2% linking efficiency. The single-exposure
//! depth comes from the shared `p9_core::analysis::surveys` table.

use p9_core::analysis::surveys::{limiting_magnitude, poisson_binomial_tail};
use p9_core::coords::observer::{Time, Timescale};
use serde::{Deserialize, Serialize};

/// Quality-cut thresholds applied to PS1 detections.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QualityCuts {
    /// Minimum PSF quality factor (psfQfPerfect > 0.99, paper Sec. 2)
    pub psf_qf_perfect_min: f64,
    /// Maximum apparent magnitude retained after quality filtering
    /// (objects fainter than V = 22.5 are rejected, paper Sec. 2)
    pub mag_max: f64,
}

/// PS1 survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ps1Survey {
    /// 50%-completeness V magnitude of the end-to-end search (21.5;
    /// abstract of Brown, Holman & Batygin 2024), shared via the p9-core
    /// survey table ("PS1 3pi").
    pub depth_limit: f64,
    /// Southern declination limit of the 3-pi survey footprint (degrees).
    pub dec_limit_deg: f64,
    /// Number of usable observing epochs per sky position. PS1 DR2 imaged
    /// each position ~12 times in each of five filters; only a subset
    /// reaches the search's effective depth. 17 is *calibrated*, not
    /// counted: with the per-epoch 50% point at `depth_limit`, requiring
    /// >= 9 of 17 epochs puts the end-to-end 50%-completeness exactly at
    /// V = 21.5 by binomial symmetry (asserted in tests).
    pub n_epochs: u32,
    /// Minimum number of detections required to link a moving object
    /// (9; "requiring 9 detections... brings the number of false positives
    /// to an acceptably low number", paper Sec. 3).
    pub linking_threshold: u32,
    /// Logistic steepness of the per-epoch efficiency roll-off (mag^-1).
    pub efficiency_steepness: f64,
    /// Quality cuts applied to individual detections
    pub quality_cuts: QualityCuts,
}

impl Default for Ps1Survey {
    fn default() -> Self {
        Self {
            depth_limit: limiting_magnitude("PS1 3pi").expect("PS1 in shared survey table"),
            dec_limit_deg: -30.0,
            n_epochs: 17,
            linking_threshold: 9,
            efficiency_steepness: 3.0,
            quality_cuts: QualityCuts {
                psf_qf_perfect_min: 0.99,
                mag_max: 22.5,
            },
        }
    }
}

impl Ps1Survey {
    /// Per-epoch detection efficiency at apparent magnitude `m`: logistic
    /// roll-off centered on the depth limit (P3 fix: previously a hard step
    /// at 21.5 next to DES's logistic, and the cadence was ignored).
    pub fn epoch_efficiency(&self, apparent_magnitude: f64) -> f64 {
        1.0 / (1.0 + ((apparent_magnitude - self.depth_limit) * self.efficiency_steepness).exp())
    }

    /// Whether a sky position falls in the 3-pi footprint (equatorial
    /// declination cut; callers must pass *declination*, not ecliptic
    /// latitude — see the sky mapping in `detection_pipeline`).
    pub fn in_footprint(&self, dec_deg: f64) -> bool {
        dec_deg > self.dec_limit_deg
    }

    /// Apply quality cuts to a detection. Returns true if the detection passes.
    pub fn passes_quality_cuts(&self, psf_qf_perfect: f64, magnitude: f64) -> bool {
        psf_qf_perfect > self.quality_cuts.psf_qf_perfect_min
            && magnitude < self.quality_cuts.mag_max
    }

    /// Tracklet-linking efficiency measured from the reference-population
    /// recovery: 99.2% of qualifying members are correctly linked
    /// (Brown, Holman & Batygin 2024).
    pub fn linking_efficiency(&self) -> f64 {
        0.992
    }

    /// Start of PS1 3-pi survey operations, 2009-06-02 (the DR2 epochs
    /// used by the Planet Nine search span 2009-2014). Anchor for
    /// `epoch_offsets_days`.
    pub fn survey_start(&self, ts: &Timescale) -> Time {
        ts.utc((2009, 6, 2))
    }

    /// End of the 3-pi imaging used by DR2 (2014-03-31).
    pub fn survey_end(&self, ts: &Timescale) -> Time {
        ts.utc((2014, 3, 31))
    }

    /// Day offsets of the `n_epochs` observing-epoch grid from
    /// `survey_start`: evenly spaced across the DR2 baseline. The
    /// detection model itself is a k-of-n binomial (per-epoch positions
    /// are not resolved); the grid exists so callers that need apparent
    /// positions (e.g. the ephemeris path in `detection_pipeline`) sample
    /// real epochs instead of raw year floats.
    pub fn epoch_offsets_days(&self) -> Vec<f64> {
        let ts = Timescale::default();
        let span = self.survey_end(&ts).tt() - self.survey_start(&ts).tt();
        let n = self.n_epochs.max(2);
        (0..n).map(|i| span * i as f64 / (n - 1) as f64).collect()
    }

    /// End-to-end probability that an object of apparent magnitude `m` at
    /// declination `dec_deg` is detected and linked:
    ///
    ///   footprint(dec) x P(Binomial(n_epochs, eps(m)) >= linking_threshold)
    ///     x linking_efficiency,
    ///
    /// with the magnitude quality cut (V < 22.5) applied. This is the
    /// cadence/linking model (P3): a faint object can be seen on a few
    /// epochs yet fail to accumulate the 9 detections linking requires.
    pub fn detection_probability(&self, apparent_magnitude: f64, dec_deg: f64) -> f64 {
        if !self.in_footprint(dec_deg) || apparent_magnitude >= self.quality_cuts.mag_max {
            return 0.0;
        }
        let eps = self.epoch_efficiency(apparent_magnitude);
        let probs = vec![eps; self.n_epochs as usize];
        poisson_binomial_tail(&probs, self.linking_threshold) * self.linking_efficiency()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_paper_values() {
        let survey = Ps1Survey::default();
        // Depth from the shared p9-core table.
        assert!((survey.depth_limit - 21.5).abs() < 1e-10);
        assert!((survey.dec_limit_deg - (-30.0)).abs() < 1e-10);
        assert_eq!(survey.linking_threshold, 9);
    }

    #[test]
    fn epoch_grid_spans_dr2_baseline() {
        let survey = Ps1Survey::default();
        let ts = Timescale::default();
        let offsets = survey.epoch_offsets_days();
        assert_eq!(offsets.len(), survey.n_epochs as usize);
        assert_eq!(offsets[0], 0.0);
        let span = survey.survey_end(&ts).tt() - survey.survey_start(&ts).tt();
        // 2009-06-02 .. 2014-03-31 is ~4.8 yr.
        assert!((span / 365.25 - 4.83).abs() < 0.02, "span = {span} days");
        assert!((offsets.last().unwrap() - span).abs() < 1e-9);
        assert!(offsets.windows(2).all(|w| w[1] > w[0]));
    }

    #[test]
    fn quality_cuts_match_paper() {
        let survey = Ps1Survey::default();
        assert!((survey.quality_cuts.psf_qf_perfect_min - 0.99).abs() < 1e-10);
        assert!((survey.quality_cuts.mag_max - 22.5).abs() < 1e-10);
    }

    #[test]
    fn epoch_efficiency_is_logistic_not_step() {
        let survey = Ps1Survey::default();
        assert!(survey.epoch_efficiency(19.0) > 0.99);
        assert!((survey.epoch_efficiency(21.5) - 0.5).abs() < 1e-10);
        assert!(survey.epoch_efficiency(23.5) < 0.01);
        let e = survey.epoch_efficiency(21.8);
        assert!(e > 0.05 && e < 0.5, "eps(21.8) = {e:.3}");
    }

    #[test]
    fn end_to_end_50_percent_at_paper_depth() {
        // Calibration regression: the k-of-n linking model's 50% point must
        // sit at the paper's V = 21.5 (P(Bin(17, 1/2) >= 9) = 1/2 exactly,
        // times the 99.2% linking efficiency).
        let survey = Ps1Survey::default();
        let p = survey.detection_probability(21.5, 10.0);
        assert!(
            (p - 0.5 * survey.linking_efficiency()).abs() < 1e-10,
            "P(detect | V=21.5) = {p}"
        );
    }

    #[test]
    fn bright_objects_nearly_always_linked() {
        let survey = Ps1Survey::default();
        let p = survey.detection_probability(19.0, 10.0);
        assert!(p > 0.98, "P = {p}");
    }

    #[test]
    fn faint_objects_fail_linking() {
        let survey = Ps1Survey::default();
        assert!(survey.detection_probability(22.4, 10.0) < 0.01);
        // The V < 22.5 quality cut zeroes anything fainter.
        assert_eq!(survey.detection_probability(22.6, 10.0), 0.0);
    }

    #[test]
    fn southern_sky_outside_footprint() {
        let survey = Ps1Survey::default();
        assert_eq!(survey.detection_probability(19.0, -45.0), 0.0);
        assert!(survey.detection_probability(19.0, -25.0) > 0.9);
    }

    #[test]
    fn quality_cuts_filter_correctly() {
        let survey = Ps1Survey::default();
        assert!(survey.passes_quality_cuts(0.995, 21.0));
        assert!(!survey.passes_quality_cuts(0.98, 21.0));
        assert!(!survey.passes_quality_cuts(0.995, 23.0));
        assert!(!survey.passes_quality_cuts(0.5, 24.0));
    }

    #[test]
    fn linking_efficiency_matches_paper() {
        let survey = Ps1Survey::default();
        assert!((survey.linking_efficiency() - 0.992).abs() < 1e-10);
    }
}
