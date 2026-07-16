//! The LSST observing-strategy model and its tunable knobs.
//!
//! Schwamb et al. (2023), "Tuning the Legacy Survey of Space and Time (LSST)
//! Observing Strategy for Solar System Science" (arXiv:2303.02355), argue that
//! the cadence and footprint choices of the Rubin/LSST main survey strongly
//! affect the discovery of distant, slow-moving solar-system bodies — Planet
//! Nine among them. Three levers dominate the discoverable fraction of slow
//! distant movers:
//!
//!   1. **Depth per visit / co-added depth.** A single `r`-band visit reaches
//!      `r ≈ 24.5` at SNR 5; the ten-year co-add of the Wide-Fast-Deep
//!      footprint reaches `r ≈ 27`. Deeper imaging recovers fainter (more
//!      distant / lower-mass) Planet Nine solutions.
//!   2. **Number of visits used for linking.** A slow mover only enters the
//!      Solar System Processing (SSP) catalog once it has been detected on
//!      enough separate nights to form linkable tracklets and a confirmed
//!      orbit. Requiring more independent detections raises the per-object
//!      completeness threshold and lowers the discoverable fraction.
//!   3. **Footprint / declination extent.** LSST images the southern sky
//!      (`δ ≲ +12°`) and avoids the Galactic plane. Shrinking the footprint
//!      (e.g. a smaller northern dec limit) removes orbits whose current sky
//!      position lies outside it.
//!
//! These published LSST numbers are kept as labelled reference constants in
//! [`published`]; the [`LsstStrategy`] struct exposes them as the knobs we
//! vary. Nothing here is a hard-coded "answer": the discoverable fraction is
//! computed from a real reference Planet Nine population in
//! [`crate::discoverable`].

use p9_core::analysis::surveys::poisson_binomial_tail;
use p9_core::units::{degrees, Angle};
use serde::{Deserialize, Serialize};

/// Published LSST single-visit and co-added depths and footprint parameters,
/// used only as labelled reference constants (never tuned to a target answer).
pub mod published {
    /// Median `r`-band single-visit 5σ point-source depth (LSST System &
    /// Survey design; Bianco et al. 2022 / Schwamb et al. 2023 baseline).
    pub const SINGLE_VISIT_DEPTH_R: f64 = 24.5;

    /// Ten-year Wide-Fast-Deep co-added `r`-band 5σ depth.
    pub const STACK_DEPTH_R: f64 = 27.0;

    /// Northern declination limit of the main survey footprint (degrees).
    /// Matches the `p9-survey` Rubin/LSST preset.
    pub const FOOTPRINT_DEC_MAX_DEG: f64 = 12.0;

    /// Southern declination limit of the main survey footprint (degrees).
    pub const FOOTPRINT_DEC_MIN_DEG: f64 = -75.0;

    /// Galactic-plane avoidance: require `|b|` above this (degrees).
    pub const GALACTIC_LAT_MIN_DEG: f64 = 10.0;

    /// Fraction of the in-band sky actually imaged to depth (WFD coverage of
    /// the declination band after dithering / rolling cadence).
    pub const COVERAGE_FRACTION: f64 = 0.85;

    /// Number of independent-night detections SSP requires to link a slow
    /// mover into a confirmed tracklet/orbit (the baseline linking threshold;
    /// Schwamb et al. discuss the sensitivity of discovery to this choice).
    pub const VISITS_FOR_LINKING: u32 = 3;

    /// Total number of independent nights a given sky position is visited
    /// over the survey, in the band used for distant-mover linking. A WFD
    /// field accumulates of order this many `r` + adjacent-band epochs that
    /// are usable for slow-mover linking over ten years.
    pub const VISITS_PER_FIELD: u32 = 30;
}

/// One LSST observing-strategy configuration. The fields are the tunable
/// levers from Schwamb et al. (2023); [`LsstStrategy::baseline`] seeds them
/// from the [`published`] reference constants.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LsstStrategy {
    /// Limiting magnitude for a *single* detection (per-visit depth).
    pub single_visit_depth: f64,
    /// Limiting magnitude of the ten-year co-add (shift-and-stack / SSP
    /// deep-drilling reach for a slow mover).
    pub stack_depth: f64,
    /// Logistic steepness of the recovery roll-off near the depth limit
    /// (mag⁻¹). Same self-calibration form used across the workspace survey
    /// crates (`p9-2021-ztf`): a few-tenths-of-a-magnitude transition.
    pub efficiency_steepness: f64,
    /// Southern declination limit of the footprint (degrees).
    pub dec_min_deg: f64,
    /// Northern declination limit of the footprint (degrees).
    pub dec_max_deg: f64,
    /// Galactic-plane avoidance: require `|b| ≥` this (degrees).
    pub galactic_lat_min_deg: f64,
    /// Fraction of the in-band sky imaged to depth.
    pub coverage_fraction: f64,
    /// Independent-night detections required to link a slow mover.
    pub visits_for_linking: u32,
    /// Independent nights a field is visited over the survey.
    pub visits_per_field: u32,
}

impl LsstStrategy {
    /// Southern declination limit as a typed [`Angle`].
    pub fn dec_min(&self) -> Angle {
        degrees(self.dec_min_deg)
    }
    /// Northern declination limit as a typed [`Angle`].
    pub fn dec_max(&self) -> Angle {
        degrees(self.dec_max_deg)
    }
    /// Galactic-plane avoidance limit as a typed [`Angle`].
    pub fn galactic_lat_min(&self) -> Angle {
        degrees(self.galactic_lat_min_deg)
    }

    /// The baseline LSST strategy from the [`published`] reference constants.
    /// Detection uses the single-visit depth (a slow mover must be detectable
    /// on individual nights to be *linked*).
    pub fn baseline() -> Self {
        Self {
            single_visit_depth: published::SINGLE_VISIT_DEPTH_R,
            stack_depth: published::STACK_DEPTH_R,
            efficiency_steepness: 4.0,
            dec_min_deg: published::FOOTPRINT_DEC_MIN_DEG,
            dec_max_deg: published::FOOTPRINT_DEC_MAX_DEG,
            galactic_lat_min_deg: published::GALACTIC_LAT_MIN_DEG,
            coverage_fraction: published::COVERAGE_FRACTION,
            visits_for_linking: published::VISITS_FOR_LINKING,
            visits_per_field: published::VISITS_PER_FIELD,
        }
    }

    /// A deep-stack variant: discovery uses the ten-year co-added depth
    /// instead of the single-visit depth (the shift-and-stack / SSP deep
    /// search). All other levers unchanged from `self`.
    pub fn with_stack_depth(mut self) -> Self {
        self.single_visit_depth = self.stack_depth;
        self
    }

    /// The same strategy with a different number of visits required to link.
    pub fn with_visits_for_linking(mut self, n: u32) -> Self {
        self.visits_for_linking = n;
        self
    }

    /// The same strategy with a different single-visit detection depth (a
    /// real cadence lever: exposure time, snaps-per-visit, and how the
    /// available dark time is split across filters all move the per-visit
    /// `r`-band depth that a slow mover must beat to be linkable).
    pub fn with_single_visit_depth(mut self, depth: f64) -> Self {
        self.single_visit_depth = depth;
        self
    }

    /// The same strategy with a different northern declination limit.
    pub fn with_dec_max(mut self, dec_max_deg: f64) -> Self {
        self.dec_max_deg = dec_max_deg;
        self
    }

    /// Single-night recovery efficiency `ε(m)`: a logistic roll-off centred on
    /// the detection depth (the self-calibration curve from injected-source
    /// recovery), identical in form to the ZTF/DES survey models.
    pub fn single_visit_efficiency(&self, apparent_magnitude: f64) -> f64 {
        let x = (apparent_magnitude - self.single_visit_depth) * self.efficiency_steepness;
        1.0 / (1.0 + x.exp())
    }

    /// Probability that a slow mover at apparent magnitude `m` and the
    /// declination band test (handled by the caller) is **linked**: it must
    /// be recovered on at least `visits_for_linking` of `visits_per_field`
    /// independent nights, each an independent Bernoulli trial at the
    /// single-visit efficiency `ε(m)`. This is the survival function of a
    /// Binomial(`visits_per_field`, `ε`) at the linking threshold.
    ///
    /// Raising `visits_for_linking` raises the threshold and so can only
    /// lower the linking probability; deepening (raising the effective
    /// detection depth, e.g. via the stack) raises `ε` and so raises it.
    pub fn linking_probability(&self, apparent_magnitude: f64) -> f64 {
        let eps = self.single_visit_efficiency(apparent_magnitude);
        poisson_binomial_tail(
            &vec![eps; self.visits_per_field as usize],
            self.visits_for_linking,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn typed_declination_accessors_match_f64() {
        use uom::si::angle::degree;
        let s = LsstStrategy::baseline();
        assert_relative_eq!(s.dec_min().get::<degree>(), s.dec_min_deg, epsilon = 1e-12);
        assert_relative_eq!(s.dec_max().get::<degree>(), s.dec_max_deg, epsilon = 1e-12);
        assert_relative_eq!(
            s.galactic_lat_min().get::<degree>(),
            s.galactic_lat_min_deg,
            epsilon = 1e-12
        );
    }

    #[test]
    fn baseline_matches_published_reference_constants() {
        let s = LsstStrategy::baseline();
        assert_eq!(s.single_visit_depth, published::SINGLE_VISIT_DEPTH_R);
        assert_eq!(s.stack_depth, published::STACK_DEPTH_R);
        assert_eq!(s.dec_max_deg, published::FOOTPRINT_DEC_MAX_DEG);
        assert_eq!(s.dec_min_deg, published::FOOTPRINT_DEC_MIN_DEG);
        assert_eq!(s.coverage_fraction, published::COVERAGE_FRACTION);
        assert_eq!(s.visits_for_linking, published::VISITS_FOR_LINKING);
    }

    #[test]
    fn single_visit_depth_matches_p9_survey_rubin_preset() {
        // The footprint band must agree with the reused p9-survey Rubin preset.
        let rubin = p9_survey::telescope::catalog()
            .into_iter()
            .find(|t| t.name.contains("Rubin"))
            .expect("Rubin preset present");
        let s = LsstStrategy::baseline();
        assert_eq!(s.dec_min_deg, rubin.footprint.dec_min_deg);
        assert_eq!(s.dec_max_deg, rubin.footprint.dec_max_deg);
        assert_eq!(s.galactic_lat_min_deg, rubin.footprint.galactic_lat_min_deg);
        assert_eq!(s.coverage_fraction, rubin.footprint.coverage_fraction);
        // The Rubin preset's stacked limiting mag sits between our single-visit
        // and ten-year-stack depths.
        assert!(rubin.limiting_mag >= s.single_visit_depth);
        assert!(rubin.limiting_mag <= s.stack_depth);
    }

    #[test]
    fn efficiency_is_logistic_half_at_depth() {
        let s = LsstStrategy::baseline();
        assert_relative_eq!(
            s.single_visit_efficiency(s.single_visit_depth),
            0.5,
            epsilon = 1e-12
        );
        assert!(s.single_visit_efficiency(s.single_visit_depth - 2.0) > 0.99);
        assert!(s.single_visit_efficiency(s.single_visit_depth + 2.0) < 0.01);
    }

    #[test]
    fn binomial_survival_matches_known_values() {
        // P[X >= 1], X ~ Bin(n, p): 1 - (1-p)^n.
        assert_relative_eq!(
            poisson_binomial_tail(&[0.1; 30], 1),
            1.0 - 0.9_f64.powi(30),
            epsilon = 1e-12
        );
        // Symmetric coin: P[X >= 3], X ~ Bin(5, 0.5) = 16/32 = 0.5.
        assert_relative_eq!(poisson_binomial_tail(&[0.5; 5], 3), 0.5, epsilon = 1e-12);
        // P[X >= 0] = 1 always; P[X >= n+1] = 0.
        assert_eq!(poisson_binomial_tail(&[0.3; 10], 0), 1.0);
        assert_eq!(poisson_binomial_tail(&[0.3; 10], 11), 0.0);
        // Edge probabilities.
        assert_eq!(poisson_binomial_tail(&[0.0; 10], 4), 0.0);
        assert_eq!(poisson_binomial_tail(&[1.0; 10], 4), 1.0);
    }

    #[test]
    fn linking_probability_drops_when_more_visits_required() {
        // A faint object near the single-visit depth: requiring more
        // independent detections monotonically reduces the linking chance.
        let mag = published::SINGLE_VISIT_DEPTH_R; // ε = 0.5 per visit
        let mut prev = 1.0;
        for n in 1..=10 {
            let p = LsstStrategy::baseline()
                .with_visits_for_linking(n)
                .linking_probability(mag);
            assert!(p <= prev + 1e-12, "linking prob rose at n = {n}");
            prev = p;
        }
    }

    #[test]
    fn deeper_depth_raises_linking_probability() {
        // A V ≈ 25.5 object is below the single-visit depth but above the
        // ten-year stack depth: the stack links it, the single visit barely.
        let mag = 25.5;
        let baseline = LsstStrategy::baseline().linking_probability(mag);
        let deep = LsstStrategy::baseline()
            .with_stack_depth()
            .linking_probability(mag);
        assert!(deep > baseline, "deep {deep:.4} !> baseline {baseline:.4}");
    }
}
