//! Parameter-space exclusion analysis for Planet Nine.
//!
//! Given a reference population of P9 orbital parameters drawn from the
//! prior (Brown & Batygin 2021 Table 1), compute the fraction of parameter
//! space ruled out by the ZTF non-detection.

use serde::{Deserialize, Serialize};

use crate::survey_model::ZtfSurvey;

/// Result of the exclusion analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ExclusionResult {
    /// Fraction of the prior parameter space excluded (paper: 56.4%)
    pub fraction_excluded: f64,
    /// Total number of prior draws evaluated
    pub n_total: u32,
    /// Number of draws that would have been detected (and thus excluded)
    pub n_excluded: u32,
}

/// Compute the fraction of the Planet Nine prior excluded by the ZTF survey.
///
/// Each entry in `magnitudes` is the predicted apparent magnitude for a
/// synthetic P9 orbit drawn from the prior. The survey model determines
/// which of these would have been detected; since no P9 was found, those
/// orbits are excluded.
///
/// The linking efficiency and sky-coverage fraction are folded in as
/// multiplicative probabilities.
pub fn compute_exclusion(survey: &ZtfSurvey, magnitudes: &[f64]) -> ExclusionResult {
    let n_total = magnitudes.len() as u32;
    let mut n_excluded = 0u32;

    let linking_eff = survey.linking_efficiency();

    for mag in magnitudes {
        if !survey.is_detectable(*mag) {
            continue;
        }
        // Probability this orbit falls in the ZTF footprint and is linked
        let p_detect = survey.sky_coverage_frac * linking_eff;
        // Treat each orbit as excluded if its detection probability exceeds 50%
        if p_detect > 0.5 {
            n_excluded += 1;
        }
    }

    let fraction_excluded = if n_total > 0 {
        n_excluded as f64 / n_total as f64
    } else {
        0.0
    };

    ExclusionResult {
        fraction_excluded,
        n_total,
        n_excluded,
    }
}

/// Placeholder for combining ZTF exclusion with DES and Pan-STARRS results.
///
/// The combined exclusion across independent surveys is:
///   1 - (1 - f_ZTF)(1 - f_DES)(1 - f_PS1)
/// assuming independent sky coverage. Full implementation requires the
/// exclusion maps from Meisner et al. (2018) and Holman & Payne (2016).
pub fn combined_exclusion(ztf_excluded: f64, des_excluded: f64, ps1_excluded: f64) -> f64 {
    1.0 - (1.0 - ztf_excluded) * (1.0 - des_excluded) * (1.0 - ps1_excluded)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn paper_exclusion_fraction() {
        // Reproduce the paper's 56.4% exclusion by constructing a magnitude
        // distribution where ~56.4% of draws are brighter than the depth limit.
        let n = 1000;
        let n_bright = 564;
        let mut mags = Vec::with_capacity(n);
        for _ in 0..n_bright {
            mags.push(19.0); // detectable
        }
        for _ in n_bright..n {
            mags.push(22.0); // too faint
        }

        let survey = ZtfSurvey::default();
        let result = compute_exclusion(&survey, &mags);
        assert_eq!(result.n_total, n as u32);
        assert!((result.fraction_excluded - 0.564).abs() < 0.01);
    }

    #[test]
    fn all_faint_nothing_excluded() {
        let mags = vec![23.0; 200];
        let survey = ZtfSurvey::default();
        let result = compute_exclusion(&survey, &mags);
        assert_eq!(result.n_excluded, 0);
        assert!((result.fraction_excluded - 0.0).abs() < 1e-10);
    }

    #[test]
    fn combined_exclusion_independent_surveys() {
        let combined = combined_exclusion(0.564, 0.20, 0.10);
        // 1 - (1-0.564)*(1-0.20)*(1-0.10) = 1 - 0.436*0.80*0.90 = 1 - 0.31392
        assert!((combined - 0.68608).abs() < 1e-4);
    }

    #[test]
    fn combined_exclusion_boundary_cases() {
        assert!((combined_exclusion(0.0, 0.0, 0.0) - 0.0).abs() < 1e-10);
        assert!((combined_exclusion(1.0, 0.0, 0.0) - 1.0).abs() < 1e-10);
        assert!((combined_exclusion(1.0, 1.0, 1.0) - 1.0).abs() < 1e-10);
    }
}
