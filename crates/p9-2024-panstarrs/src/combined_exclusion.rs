//! Multi-survey combined exclusion analysis.
//!
//! Combines the exclusion fractions from ZTF, DES, and PS1 to compute the
//! total fraction of Planet Nine parameter space ruled out. The three-survey
//! combination excludes ~78% of the prior, leaving ~22% viable.

use serde::{Deserialize, Serialize};

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;

/// Updated Planet Nine parameter estimates from the combined analysis.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct UpdatedParameters {
    /// Semi-major axis in AU (median)
    pub a_median: f64,
    /// Semi-major axis upper uncertainty (+170 AU)
    pub a_upper: f64,
    /// Semi-major axis lower uncertainty (-120 AU)
    pub a_lower: f64,
    /// Mass in Earth masses (median)
    pub mass_earth_median: f64,
    /// Mass upper uncertainty (+2.6 ME)
    pub mass_upper: f64,
    /// Mass lower uncertainty (-1.7 ME)
    pub mass_lower: f64,
    /// Apparent V-band magnitude (median)
    pub v_mag_median: f64,
    /// V magnitude upper uncertainty (+1.1)
    pub v_mag_upper: f64,
    /// V magnitude lower uncertainty (-1.4)
    pub v_mag_lower: f64,
}

impl UpdatedParameters {
    /// Paper's updated parameter estimates after three-survey exclusion.
    pub fn paper_values() -> Self {
        Self {
            a_median: 500.0,
            a_upper: 170.0,
            a_lower: 120.0,
            mass_earth_median: 6.6,
            mass_upper: 2.6,
            mass_lower: 1.7,
            v_mag_median: 22.0,
            v_mag_upper: 1.1,
            v_mag_lower: 1.4,
        }
    }

    /// Convert median parameters to a P9Params for orbital computations.
    pub fn to_p9_params(&self) -> P9Params {
        P9Params {
            mass_earth: self.mass_earth_median,
            a: self.a_median,
            e: 0.3,
            i: 16.0 * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }
}

/// Combined exclusion fractions from the three-survey analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CombinedExclusion {
    /// Fraction excluded by ZTF alone (56.4%)
    pub ztf_frac: f64,
    /// Additional unique fraction excluded by DES (5.0%)
    pub des_unique: f64,
    /// Additional unique fraction excluded by PS1 (17.1%)
    pub ps1_unique: f64,
    /// Total combined exclusion (78%)
    pub combined: f64,
}

impl CombinedExclusion {
    /// Paper's reported exclusion fractions.
    pub fn paper_values() -> Self {
        Self {
            ztf_frac: 0.564,
            des_unique: 0.050,
            ps1_unique: 0.171,
            combined: 0.78,
        }
    }
}

/// Combine exclusion fractions from ZTF, DES, and PS1.
///
/// The surveys have overlapping sky coverage, so the unique contributions
/// from DES and PS1 are added to the ZTF baseline. The combined fraction
/// is capped at 1.0.
pub fn compute_combined(ztf_frac: f64, des_unique: f64, ps1_unique: f64) -> CombinedExclusion {
    let combined = (ztf_frac + des_unique + ps1_unique).min(1.0);
    CombinedExclusion {
        ztf_frac,
        des_unique,
        ps1_unique,
        combined,
    }
}

/// Describe the remaining viable parameter space after combined exclusion.
///
/// Returns the fraction of the prior that has not been excluded by any survey.
pub fn remaining_parameter_space(exclusion: &CombinedExclusion) -> f64 {
    1.0 - exclusion.combined
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn paper_exclusion_values() {
        let ex = CombinedExclusion::paper_values();
        assert_relative_eq!(ex.ztf_frac, 0.564, epsilon = 1e-10);
        assert_relative_eq!(ex.des_unique, 0.050, epsilon = 1e-10);
        assert_relative_eq!(ex.ps1_unique, 0.171, epsilon = 1e-10);
        assert_relative_eq!(ex.combined, 0.78, epsilon = 1e-10);
    }

    #[test]
    fn compute_combined_matches_paper() {
        let ex = compute_combined(0.564, 0.050, 0.171);
        assert_relative_eq!(ex.combined, 0.785, epsilon = 1e-10);
    }

    #[test]
    fn remaining_space_is_complement() {
        let ex = CombinedExclusion::paper_values();
        let remaining = remaining_parameter_space(&ex);
        assert_relative_eq!(remaining, 0.22, epsilon = 1e-10);
    }

    #[test]
    fn combined_capped_at_unity() {
        let ex = compute_combined(0.8, 0.2, 0.3);
        assert_relative_eq!(ex.combined, 1.0, epsilon = 1e-10);
    }

    #[test]
    fn updated_parameters_match_paper() {
        let params = UpdatedParameters::paper_values();
        assert_relative_eq!(params.a_median, 500.0, epsilon = 1e-10);
        assert_relative_eq!(params.a_upper, 170.0, epsilon = 1e-10);
        assert_relative_eq!(params.a_lower, 120.0, epsilon = 1e-10);
        assert_relative_eq!(params.mass_earth_median, 6.6, epsilon = 1e-10);
        assert_relative_eq!(params.mass_upper, 2.6, epsilon = 1e-10);
        assert_relative_eq!(params.mass_lower, 1.7, epsilon = 1e-10);
        assert_relative_eq!(params.v_mag_median, 22.0, epsilon = 1e-10);
        assert_relative_eq!(params.v_mag_upper, 1.1, epsilon = 1e-10);
        assert_relative_eq!(params.v_mag_lower, 1.4, epsilon = 1e-10);
    }

    #[test]
    fn updated_params_to_p9_params() {
        let params = UpdatedParameters::paper_values();
        let p9 = params.to_p9_params();
        assert_relative_eq!(p9.mass_earth, 6.6, epsilon = 1e-10);
        assert_relative_eq!(p9.a, 500.0, epsilon = 1e-10);
    }
}
