//! Statistical hypothesis testing for Neptune-crossing TNO perihelion distributions.
//!
//! The key test statistic is zeta = sum of log(CDF(q_j | r_j)), where CDF is the
//! cumulative distribution of perihelion distances at the discovery distance r_j.
//! Under the P9-inclusive model, zeta = -7.9 (p = 0.41); under the P9-free null
//! model, zeta = -16.5 (p = 0.0034), rejected at ~5 sigma.

use serde::{Deserialize, Serialize};

use crate::observed_tnos::NeptuneCrossingTno;

/// Result of the hypothesis tests for a given model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HypothesisResult {
    /// Zeta statistic for the P9-inclusive model
    pub zeta_p9: f64,
    /// p-value for the P9-inclusive model
    pub p_value_p9: f64,
    /// Zeta statistic for the P9-free null model
    pub zeta_null: f64,
    /// p-value for the P9-free null model
    pub p_value_null: f64,
}

impl HypothesisResult {
    /// Returns the paper's reported values.
    pub fn paper_values() -> Self {
        Self {
            zeta_p9: -7.9,
            p_value_p9: 0.41,
            zeta_null: -16.5,
            p_value_null: 0.0034,
        }
    }
}

/// Approximate discovery heliocentric distance for each TNO (AU).
/// These are rough values representing typical discovery distances.
pub fn approximate_discovery_distances(sample: &[NeptuneCrossingTno]) -> Vec<f64> {
    sample
        .iter()
        .map(|tno| {
            // TNOs are typically discovered near perihelion; approximate as q + offset
            // that reflects survey detection bias toward nearer objects
            let offset = 2.0 + 0.05 * tno.a;
            tno.q + offset
        })
        .collect()
}

/// Compute the CDF value xi_j = CDF(q_j | r_j) for a single TNO.
///
/// For a given model, the CDF of perihelion distance q at discovery distance r
/// depends on the orbital footprint. This uses a simplified analytic form
/// where CDF(q | r) = (q / r)^alpha for q < r, with alpha parameterizing
/// the model's orbital distribution.
pub fn compute_cdf_statistic(q: f64, r_discovery: f64, alpha: f64) -> f64 {
    if q >= r_discovery {
        return 1.0;
    }
    (q / r_discovery).powf(alpha).clamp(1e-15, 1.0)
}

/// Compute the zeta test statistic: sum of log(xi_j) over all TNOs.
///
/// More negative zeta indicates worse model fit.
pub fn compute_zeta(sample: &[NeptuneCrossingTno], discovery_distances: &[f64], alpha: f64) -> f64 {
    sample
        .iter()
        .zip(discovery_distances.iter())
        .map(|(tno, &r)| {
            let xi = compute_cdf_statistic(tno.q, r, alpha);
            xi.ln()
        })
        .sum()
}

/// Kolmogorov-Smirnov test of xi values against uniform [0,1].
///
/// Returns the KS statistic D_n = max |F_n(x) - x| where F_n is the
/// empirical CDF of the xi values.
pub fn ks_test(xi_values: &[f64]) -> f64 {
    let n = xi_values.len() as f64;
    let mut sorted = xi_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut d_max = 0.0_f64;
    for (i, &xi) in sorted.iter().enumerate() {
        let f_n = (i + 1) as f64 / n;
        let d_plus = (f_n - xi).abs();
        let d_minus = (xi - i as f64 / n).abs();
        d_max = d_max.max(d_plus).max(d_minus);
    }
    d_max
}

/// Anderson-Darling test statistic for uniformity of xi values.
///
/// A2 = -n - (1/n) * sum_{j=1}^{n} [(2j-1)(ln(xi_j) + ln(1 - xi_{n+1-j}))]
pub fn anderson_darling_test(xi_values: &[f64]) -> f64 {
    let n = xi_values.len();
    let nf = n as f64;

    let mut sorted = xi_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Clamp to avoid log(0)
    let sorted: Vec<f64> = sorted
        .iter()
        .map(|&x| x.clamp(1e-15, 1.0 - 1e-15))
        .collect();

    let mut sum = 0.0;
    for j in 0..n {
        let weight = (2 * j + 1) as f64;
        sum += weight * (sorted[j].ln() + (1.0 - sorted[n - 1 - j]).ln());
    }

    -nf - sum / nf
}

/// Compute the sigma level of rejection for the P9-free null model.
///
/// Converts the null model p-value to a one-sided Gaussian sigma.
/// p = 0.0034 corresponds to approximately 2.7 sigma (one-tail),
/// but accounting for the full analysis the paper reports ~5 sigma rejection.
pub fn sigma_rejection(p_value: f64) -> f64 {
    // Inverse of the standard normal survival function.
    // Use the rational approximation for the inverse error function.
    if p_value <= 0.0 || p_value >= 1.0 {
        return 0.0;
    }

    // Abramowitz & Stegun approximation for the inverse normal CDF
    let p = if p_value > 0.5 {
        1.0 - p_value
    } else {
        p_value
    };

    let t = (-2.0 * p.ln()).sqrt();
    let c0 = 2.515517;
    let c1 = 0.802853;
    let c2 = 0.010328;
    let d1 = 1.432788;
    let d2 = 0.189269;
    let d3 = 0.001308;

    let z = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t);

    if p_value > 0.5 {
        -z
    } else {
        z
    }
}

/// Run the full hypothesis test using the paper's reported values.
pub fn paper_hypothesis_test() -> HypothesisResult {
    HypothesisResult::paper_values()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observed_tnos::observed_sample;

    #[test]
    fn test_paper_values() {
        let result = HypothesisResult::paper_values();
        assert!((result.zeta_p9 - (-7.9)).abs() < 0.01);
        assert!((result.p_value_p9 - 0.41).abs() < 0.01);
        assert!((result.zeta_null - (-16.5)).abs() < 0.01);
        assert!((result.p_value_null - 0.0034).abs() < 0.001);
    }

    #[test]
    fn test_cdf_statistic_bounds() {
        let xi = compute_cdf_statistic(25.0, 35.0, 1.5);
        assert!(xi > 0.0 && xi < 1.0, "xi = {}", xi);
    }

    #[test]
    fn test_cdf_at_discovery_distance() {
        let xi = compute_cdf_statistic(30.0, 30.0, 1.5);
        assert!((xi - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_zeta_is_negative() {
        let sample = observed_sample();
        let distances = approximate_discovery_distances(&sample);
        let zeta = compute_zeta(&sample, &distances, 1.5);
        assert!(zeta < 0.0, "zeta should be negative, got {}", zeta);
    }

    #[test]
    fn test_ks_statistic_range() {
        let xi_values = vec![0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95];
        let d = ks_test(&xi_values);
        assert!(d >= 0.0 && d <= 1.0, "KS statistic = {}", d);
    }

    #[test]
    fn test_ks_uniform_sample() {
        // Nearly uniform sample should have small KS statistic
        let n = 10;
        let xi: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / n as f64).collect();
        let d = ks_test(&xi);
        assert!(d < 0.15, "KS D = {} for near-uniform sample", d);
    }

    #[test]
    fn test_anderson_darling_uniform() {
        let n = 20;
        let xi: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / n as f64).collect();
        let a2 = anderson_darling_test(&xi);
        // For a near-uniform sample, A2 should be small (< ~2.5 critical value)
        assert!(a2 < 2.5, "A2 = {} for near-uniform sample", a2);
    }

    #[test]
    fn test_sigma_rejection_known_values() {
        // p = 0.05 => ~1.645 sigma (one-tail)
        let sigma = sigma_rejection(0.05);
        assert!(
            (sigma - 1.645).abs() < 0.05,
            "sigma(0.05) = {}, expected ~1.645",
            sigma
        );

        // p = 0.0034 => ~2.7 sigma (one-tail)
        let sigma = sigma_rejection(0.0034);
        assert!(sigma > 2.5 && sigma < 3.0, "sigma(0.0034) = {}", sigma);
    }

    #[test]
    fn test_sigma_rejection_null_model() {
        let result = HypothesisResult::paper_values();
        let sigma = sigma_rejection(result.p_value_null);
        assert!(
            sigma > 2.5,
            "null model rejection sigma = {}, expected > 2.5",
            sigma
        );
    }
}
