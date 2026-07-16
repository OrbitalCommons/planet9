//! Conversion of the log-likelihood ratio to a Gaussian-equivalent
//! significance in σ.
//!
//! Under the isotropic null H0, Wilks' theorem states that 2Λ is
//! asymptotically distributed as χ² with k = 2 degrees of freedom, where k is
//! the number of free parameters of the von Mises model (µ and κ). For k = 2
//! the χ² survival function is closed form:
//!
//!   p = P(χ²₂ ≥ 2Λ) = exp(−Λ).
//!
//! We treat p as a two-sided Gaussian tail probability and map it to a
//! Gaussian-equivalent significance via the inverse normal CDF:
//!
//!   σ = Φ⁻¹(1 − p/2),
//!
//! so that p = 0.05 ↔ 1.96σ and p = 0.0027 ↔ 3σ. The normal quantile uses
//! the Acklam (2003) rational approximation, accurate to ~1e-9 relative.

use p9_core::analysis::stats::normal_quantile;

/// χ²-with-2-dof survival function P(χ²₂ ≥ x) = exp(−x/2).
fn chi2_2dof_survival(x: f64) -> f64 {
    (-0.5 * x).exp()
}

/// Two-sided p-value for the apsidal-clustering log-likelihood ratio via
/// Wilks' theorem with 2 degrees of freedom: p = P(χ²₂ ≥ 2Λ) = exp(−Λ).
pub fn lambda_to_p_value(lambda: f64) -> f64 {
    let lambda = lambda.max(0.0);
    chi2_2dof_survival(2.0 * lambda).clamp(f64::MIN_POSITIVE, 1.0)
}

/// Convert a two-sided p-value to a Gaussian-equivalent significance in σ:
///
///   σ = Φ⁻¹(1 − p/2).
pub fn p_value_to_sigma(p: f64) -> f64 {
    let p = p.clamp(f64::MIN_POSITIVE, 1.0 - 1e-15);
    normal_quantile(1.0 - 0.5 * p)
}

/// Full pipeline: log-likelihood ratio Λ → Wilks p-value → Gaussian σ.
pub fn sigma(lambda: f64) -> f64 {
    p_value_to_sigma(lambda_to_p_value(lambda))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_known_points() {
        // Median.
        assert!((normal_quantile(0.5)).abs() < 1e-9);
        // 97.5th percentile ≈ 1.95996.
        assert!((normal_quantile(0.975) - 1.959_963_98).abs() < 1e-6);
        // 99.865th percentile ≈ 3.0 (one-sided 3σ).
        assert!((normal_quantile(0.998_650_1) - 3.0).abs() < 1e-3);
    }

    #[test]
    fn quantile_symmetry() {
        for p in [0.01, 0.1, 0.3, 0.45] {
            assert!((normal_quantile(p) + normal_quantile(1.0 - p)).abs() < 1e-8);
        }
    }

    #[test]
    fn p_to_sigma_known_values() {
        // p = 0.05 (two-sided) ↔ 1.96σ.
        let s = p_value_to_sigma(0.05);
        assert!((s - 1.959_963_98).abs() < 1e-4, "σ(0.05) = {s}");
        // p = 0.0027 (two-sided) ↔ 3σ.
        let s = p_value_to_sigma(0.0027);
        assert!((s - 3.0).abs() < 2e-3, "σ(0.0027) = {s}");
    }

    #[test]
    fn sigma_monotone_in_lambda() {
        assert!(sigma(1.0) < sigma(2.0));
        assert!(sigma(2.0) < sigma(5.0));
        // Λ = 0 → p = 1 → σ = 0.
        assert!(sigma(0.0) < 1e-6);
    }
}
