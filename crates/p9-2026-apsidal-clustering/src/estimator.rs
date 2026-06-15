//! Conditional-likelihood estimator for apsidal clustering.
//!
//! The estimator compares two models for the longitude of perihelion ϖ:
//!
//! - H1 (clustered): ϖ ~ von Mises(µ, κ), with (µ, κ) maximum-likelihood
//!   fitted from the sample — µ = circular mean, κ = A⁻¹(R̄).
//! - H0 (isotropic): ϖ uniform on the circle, density 1/2π.
//!
//! The test statistic is the log-likelihood ratio
//!
//!   Λ = ln L(H1) − ln L(H0) = Σ_i [ ln f(ϖ_i | µ, κ) − ln(1/2π) ]
//!     = Σ_i [ κ·cos(ϖ_i − µ) − ln I₀(κ) ]
//!
//! (the 1/2π normalizations cancel). Λ ≥ 0 always, since H0 is nested in H1
//! at κ = 0. Under H0, by Wilks' theorem 2Λ is asymptotically χ² with 2
//! degrees of freedom (the two free parameters µ, κ of the von Mises model);
//! [`crate::significance`] performs the χ² → p → σ conversion.

use p9_core::analysis::circular::{
    bessel_i0, circular_mean, kappa_from_r_bar, mean_resultant_length,
};
use serde::{Deserialize, Serialize};

/// Maximum-likelihood von Mises fit of a ϖ sample plus the resulting
/// log-likelihood-ratio statistic against the isotropic null.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ApsidalFit {
    /// Sample size.
    pub n: usize,
    /// Maximum-likelihood mean direction µ (rad).
    pub mu: f64,
    /// Maximum-likelihood concentration κ.
    pub kappa: f64,
    /// Mean resultant length R̄ of the sample.
    pub r_bar: f64,
    /// Log-likelihood ratio Λ = ln L(H1) − ln L(H0).
    pub lambda: f64,
}

/// Maximum-likelihood von Mises fit (µ = circular mean, κ = A⁻¹(R̄)).
pub fn von_mises_fit(varpi: &[f64]) -> (f64, f64) {
    let mu = circular_mean(varpi).unwrap_or(0.0);
    let r_bar = mean_resultant_length(varpi);
    (mu, kappa_from_r_bar(r_bar))
}

/// Log-likelihood ratio Λ = ln L(H1) − ln L(H0) for the apsidal-clustering
/// von Mises (H1) vs uniform (H0) models.
///
/// Both models share the 1/2π normalization, so it cancels term-by-term:
///
///   Λ = Σ_i [ κ·cos(ϖ_i − µ) − ln I₀(κ) ].
pub fn log_likelihood_ratio(varpi: &[f64], mu: f64, kappa: f64) -> f64 {
    let ln_i0 = bessel_i0(kappa).ln();
    varpi.iter().map(|&w| kappa * (w - mu).cos() - ln_i0).sum()
}

/// Fit the apsidal-clustering estimator to a ϖ sample and assemble the
/// [`ApsidalFit`] summary.
pub fn fit(varpi: &[f64]) -> ApsidalFit {
    let (mu, kappa) = von_mises_fit(varpi);
    let r_bar = mean_resultant_length(varpi);
    let lambda = log_likelihood_ratio(varpi, mu, kappa);
    ApsidalFit {
        n: varpi.len(),
        mu,
        kappa,
        r_bar,
        lambda,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn lambda_nonnegative_for_clustered() {
        // A clustered sample must prefer H1 (Λ ≥ 0).
        let varpi: Vec<f64> = (0..10).map(|k| 1.0 + 0.05 * k as f64).collect();
        let (mu, kappa) = von_mises_fit(&varpi);
        let lambda = log_likelihood_ratio(&varpi, mu, kappa);
        assert!(lambda > 0.0, "Λ = {lambda}");
    }

    #[test]
    fn lambda_near_zero_for_uniform() {
        // An evenly spread sample has R̄ ≈ 0, κ ≈ 0, Λ ≈ 0.
        let varpi: Vec<f64> = (0..12).map(|k| k as f64 * TAU / 12.0).collect();
        let fit = fit(&varpi);
        assert!(fit.kappa < 0.1, "κ = {}", fit.kappa);
        assert!(fit.lambda.abs() < 0.5, "Λ = {}", fit.lambda);
    }

    #[test]
    fn fit_recovers_mean_direction() {
        let varpi: Vec<f64> = (0..20).map(|k| 2.0 + 0.02 * (k as f64 - 10.0)).collect();
        let fit = fit(&varpi);
        assert!((fit.mu - 2.0).abs() < 0.05, "µ = {}", fit.mu);
        assert!(fit.kappa > 1.0, "κ = {}", fit.kappa);
    }
}
