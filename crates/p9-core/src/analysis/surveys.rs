//! Shared survey table: limiting magnitudes and published Planet Nine
//! exclusion fractions.
//!
//! Single source of truth — the ZTF baseline previously appeared as 0.564,
//! 0.56 and 0.612-combined in three crates, and two crates disagreed about
//! the DES depth (23.8 vs 24.1).

/// A survey's published contribution to the Planet Nine parameter-space
/// exclusion.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurveyExclusion {
    pub survey: &'static str,
    /// Fraction of the reference P9 population excluded *uniquely* by this
    /// survey (i.e. not already excluded by the surveys above it in the
    /// table).
    pub unique_fraction: f64,
    pub reference: &'static str,
}

/// Published unique exclusion fractions, in the combination order of
/// Belyakov, Bernardinelli & Brown (2024): ZTF first, then DES, then PS1.
pub const EXCLUSION_FRACTIONS: [SurveyExclusion; 3] = [
    SurveyExclusion {
        survey: "ZTF",
        unique_fraction: 0.564,
        reference: "Brown & Batygin (2022), ZTF non-detection",
    },
    SurveyExclusion {
        survey: "DES",
        unique_fraction: 0.050,
        reference: "Bernardinelli et al. (2022) via Belyakov et al. (2024)",
    },
    SurveyExclusion {
        survey: "PS1",
        unique_fraction: 0.171,
        reference: "Belyakov, Bernardinelli & Brown (2024)",
    },
];

/// Combine *unique* (already-disjoint) exclusion fractions: a simple sum,
/// capped at 1. Note the rigorous combination is a per-orbit OR over a common
/// reference population — use these published unique fractions only as the
/// papers report them.
pub fn combined_unique_exclusion(unique_fractions: &[f64]) -> f64 {
    unique_fractions.iter().sum::<f64>().min(1.0)
}

/// Combine exclusion fractions of *independent* surveys: 1 − Π(1 − f).
/// Only valid when the surveys' excluded regions are uncorrelated, which sky
/// surveys generally are not; prefer per-orbit combination.
pub fn combined_independent_exclusion(fractions: &[f64]) -> f64 {
    1.0 - fractions.iter().map(|f| 1.0 - f).product::<f64>()
}

/// A survey's approximate single-epoch limiting magnitude.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurveyDepth {
    pub survey: &'static str,
    /// Limiting magnitude (band noted in `band`)
    pub limiting_mag: f64,
    pub band: &'static str,
    pub reference: &'static str,
}

/// Approximate survey depths used across the detectability crates.
pub const SURVEY_DEPTHS: [SurveyDepth; 5] = [
    SurveyDepth {
        survey: "WISE W1",
        limiting_mag: 16.5,
        band: "W1",
        reference: "Wright et al. (2010)",
    },
    SurveyDepth {
        survey: "CRTS",
        limiting_mag: 19.5,
        band: "V",
        reference: "Drake et al. (2009)",
    },
    SurveyDepth {
        survey: "ZTF",
        limiting_mag: 20.5,
        band: "r",
        reference: "Brown & Batygin (2022)",
    },
    SurveyDepth {
        survey: "PS1 3pi",
        limiting_mag: 21.5,
        band: "r",
        reference: "Chambers et al. (2016)",
    },
    SurveyDepth {
        survey: "DES",
        limiting_mag: 23.8,
        band: "r",
        reference: "Bernardinelli et al. (2022)",
    },
];

/// Look up a survey's limiting magnitude by name.
pub fn limiting_magnitude(survey: &str) -> Option<f64> {
    SURVEY_DEPTHS
        .iter()
        .find(|d| d.survey.eq_ignore_ascii_case(survey))
        .map(|d| d.limiting_mag)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_combined_unique_matches_paper() {
        // Belyakov et al. (2024): 0.564 + 0.050 + 0.171 = 0.785 ≈ 78%.
        let fracs: Vec<f64> = EXCLUSION_FRACTIONS
            .iter()
            .map(|s| s.unique_fraction)
            .collect();
        let combined = combined_unique_exclusion(&fracs);
        assert_relative_eq!(combined, 0.785, epsilon = 1e-12);
    }

    #[test]
    fn test_combined_unique_caps_at_one() {
        assert_eq!(combined_unique_exclusion(&[0.7, 0.6]), 1.0);
    }

    #[test]
    fn test_independent_combination() {
        let c = combined_independent_exclusion(&[0.5, 0.5]);
        assert_relative_eq!(c, 0.75, epsilon = 1e-12);
        assert_eq!(combined_independent_exclusion(&[]), 0.0);
    }

    #[test]
    fn test_depth_lookup() {
        assert_eq!(limiting_magnitude("DES"), Some(23.8));
        assert_eq!(limiting_magnitude("ztf"), Some(20.5));
        assert_eq!(limiting_magnitude("LSST"), None);
    }

    #[test]
    fn test_fractions_are_probabilities() {
        for s in &EXCLUSION_FRACTIONS {
            assert!((0.0..=1.0).contains(&s.unique_fraction), "{}", s.survey);
        }
    }
}
