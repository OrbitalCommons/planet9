//! Three-class stability classification based on diffusion coefficients.
//!
//! D_a,crit = 10^{-3} AU^2/yr
//! Stable: <D_a> < 0.5 * D_crit AND sigma < <D_a>
//! Metastable: <D_a> < 5 * D_crit
//! Unstable: everything else

use serde::{Deserialize, Serialize};

use crate::diffusion::DiffusionResult;

/// Critical diffusion coefficient (AU^2/yr).
pub const D_CRIT: f64 = 1e-3;

/// Stability classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum StabilityClass {
    Stable,
    Metastable,
    Unstable,
}

/// Classify a TNO based on its diffusion measurement.
pub fn classify(result: &DiffusionResult) -> StabilityClass {
    if result.d_mean < 0.5 * D_CRIT && result.d_std < result.d_mean {
        StabilityClass::Stable
    } else if result.d_mean < 5.0 * D_CRIT {
        StabilityClass::Metastable
    } else {
        StabilityClass::Unstable
    }
}

/// Classify a batch of TNOs and return counts per class.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClassificationSummary {
    pub n_stable: usize,
    pub n_metastable: usize,
    pub n_unstable: usize,
}

pub fn classify_batch(results: &[DiffusionResult]) -> ClassificationSummary {
    let mut summary = ClassificationSummary {
        n_stable: 0,
        n_metastable: 0,
        n_unstable: 0,
    };

    for r in results {
        match classify(r) {
            StabilityClass::Stable => summary.n_stable += 1,
            StabilityClass::Metastable => summary.n_metastable += 1,
            StabilityClass::Unstable => summary.n_unstable += 1,
        }
    }

    summary
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_result(d_mean: f64, d_std: f64) -> DiffusionResult {
        DiffusionResult {
            name: "test",
            d_mean,
            d_std,
            d_analytical: 0.0,
        }
    }

    #[test]
    fn test_classify_stable() {
        let r = make_result(1e-4, 5e-5);
        assert_eq!(classify(&r), StabilityClass::Stable);
    }

    #[test]
    fn test_classify_metastable() {
        let r = make_result(2e-3, 1e-3);
        assert_eq!(classify(&r), StabilityClass::Metastable);
    }

    #[test]
    fn test_classify_unstable() {
        let r = make_result(0.1, 0.05);
        assert_eq!(classify(&r), StabilityClass::Unstable);
    }

    #[test]
    fn test_classify_batch() {
        let results = vec![
            make_result(1e-4, 5e-5), // stable
            make_result(2e-3, 1e-3), // metastable
            make_result(0.1, 0.05),  // unstable
            make_result(1e-5, 1e-6), // stable
        ];
        let summary = classify_batch(&results);
        assert_eq!(summary.n_stable, 2);
        assert_eq!(summary.n_metastable, 1);
        assert_eq!(summary.n_unstable, 1);
    }
}
