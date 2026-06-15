//! Survey model for the Holman et al. (2025) shift-and-stack PS1 search.
//!
//! The single-epoch PS1 3-pi r-band depth (21.5) comes from the shared
//! `p9_core::analysis::surveys` table. The distinguishing feature of this
//! search, relative to the catalog-linking search in `p9-2024-panstarrs`,
//! is the *shift-and-stack* depth gain: co-adding the N multi-epoch
//! exposures along a trial motion vector lowers the noise by ~sqrt(N),
//! pushing the effective 50%-completeness magnitude fainter by
//! `stack_depth_gain` mag. The effective depth is therefore
//!
//!   m_eff = m_single + stack_depth_gain.
//!
//! Detection probability is a logistic roll-off centered on `m_eff`,
//! gated by the declination > -30 deg 3-pi footprint.

use p9_core::analysis::surveys::limiting_magnitude;
use serde::{Deserialize, Serialize};

use crate::published;

/// PS1 3-pi shift-and-stack survey model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Ps1StackSurvey {
    /// Single-epoch 50%-completeness r magnitude (21.5), from the shared
    /// `p9-core` "PS1 3pi" survey-depth table.
    pub single_epoch_depth: f64,
    /// Magnitudes of effective depth gained by the shift-and-stack co-add
    /// (~sqrt(N) noise reduction over N epochs). Holman et al. (2025) reach
    /// ~1 mag deeper than a single exposure.
    pub stack_depth_gain: f64,
    /// Southern declination limit of the 3-pi footprint (degrees).
    pub dec_limit_deg: f64,
    /// Logistic steepness of the efficiency roll-off (mag^-1).
    pub efficiency_steepness: f64,
}

impl Default for Ps1StackSurvey {
    fn default() -> Self {
        Self {
            single_epoch_depth: limiting_magnitude("PS1 3pi")
                .expect("PS1 3pi in shared survey table"),
            stack_depth_gain: published::SHIFT_STACK_DEPTH_GAIN_MAG,
            dec_limit_deg: published::PS1_DEC_LIMIT_DEG,
            efficiency_steepness: 3.0,
        }
    }
}

impl Ps1StackSurvey {
    /// Effective 50%-completeness depth after the shift-and-stack co-add:
    /// `single_epoch_depth + stack_depth_gain`.
    pub fn effective_depth(&self) -> f64 {
        self.single_epoch_depth + self.stack_depth_gain
    }

    /// Whether an equatorial declination falls in the 3-pi footprint.
    pub fn in_footprint(&self, dec_deg: f64) -> bool {
        dec_deg > self.dec_limit_deg
    }

    /// Logistic recovery efficiency at apparent r magnitude `m`, centered on
    /// the effective (stacked) depth. 0.5 exactly at `effective_depth()`.
    pub fn recovery_efficiency(&self, apparent_magnitude: f64) -> f64 {
        1.0 / (1.0
            + ((apparent_magnitude - self.effective_depth()) * self.efficiency_steepness).exp())
    }

    /// End-to-end detection probability for an object of apparent r
    /// magnitude `m` at equatorial declination `dec_deg`: the footprint cut
    /// times the logistic recovery efficiency.
    pub fn detection_probability(&self, apparent_magnitude: f64, dec_deg: f64) -> f64 {
        if !self.in_footprint(dec_deg) {
            return 0.0;
        }
        self.recovery_efficiency(apparent_magnitude)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depth_anchored_to_shared_table() {
        let s = Ps1StackSurvey::default();
        // PINNED: PS1 3-pi single-epoch r depth == 21.5 from p9-core.
        assert!((s.single_epoch_depth - 21.5).abs() < 1e-12);
        assert!((s.single_epoch_depth - limiting_magnitude("PS1 3pi").unwrap()).abs() < 1e-12);
    }

    #[test]
    fn footprint_is_dec_minus_30() {
        let s = Ps1StackSurvey::default();
        assert!((s.dec_limit_deg - (-30.0)).abs() < 1e-12);
        assert!(s.in_footprint(-29.0));
        assert!(!s.in_footprint(-31.0));
        assert_eq!(s.detection_probability(20.0, -45.0), 0.0);
    }

    #[test]
    fn shift_stack_deepens_the_search() {
        let s = Ps1StackSurvey::default();
        // The stack reaches fainter than the single-epoch depth.
        assert!(s.effective_depth() > s.single_epoch_depth);
        assert!((s.effective_depth() - 22.5).abs() < 1e-12);
    }

    #[test]
    fn efficiency_half_at_effective_depth() {
        let s = Ps1StackSurvey::default();
        let e = s.recovery_efficiency(s.effective_depth());
        assert!((e - 0.5).abs() < 1e-12, "eff(m_eff) = {e}");
        assert!(s.recovery_efficiency(19.0) > 0.99);
        assert!(s.recovery_efficiency(24.5) < 0.01);
    }

    #[test]
    fn stack_detects_fainter_than_single_epoch() {
        // At r = 22.2 (fainter than single-epoch 21.5) the stack still
        // recovers most objects; without the gain it would be sub-threshold.
        let stack = Ps1StackSurvey::default();
        let no_stack = Ps1StackSurvey {
            stack_depth_gain: 0.0,
            ..Ps1StackSurvey::default()
        };
        let p_stack = stack.detection_probability(22.2, 0.0);
        let p_single = no_stack.detection_probability(22.2, 0.0);
        assert!(p_stack > 0.6, "stack p = {p_stack}");
        assert!(p_single < 0.2, "single-epoch p = {p_single}");
        assert!(p_stack > p_single);
    }
}
