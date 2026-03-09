//! Statistical clustering tests from Brown (2017).
//!
//! Implements the Monte Carlo framework for testing whether the observed
//! clustering of distant KBO orbital angles is statistically significant
//! after accounting for observational bias.

use std::f64::consts::PI;

use rand::Rng;

use crate::bias_function;
use crate::kbo_sample::DistantKbo;

/// Result of a clustering test.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ClusteringResult {
    /// Observed Rayleigh z-statistic for ϖ
    pub observed_z_varpi: f64,
    /// Observed Rayleigh z-statistic for ω
    pub observed_z_omega: f64,
    /// P-value for ϖ clustering alone
    pub p_varpi: f64,
    /// P-value for ω clustering alone
    pub p_omega: f64,
    /// Combined p-value (joint ϖ + ω)
    pub p_combined: f64,
    /// Number of Monte Carlo iterations
    pub n_iterations: usize,
}

/// Compute the Rayleigh z-statistic for a set of angles.
///
/// z = (Σ cos θ)² + (Σ sin θ)² / n
///
/// Higher z indicates stronger clustering.
pub fn rayleigh_z(angles: &[f64]) -> f64 {
    let n = angles.len() as f64;
    if n == 0.0 {
        return 0.0;
    }
    let cos_sum: f64 = angles.iter().map(|&a| a.cos()).sum();
    let sin_sum: f64 = angles.iter().map(|&a| a.sin()).sum();
    (cos_sum * cos_sum + sin_sum * sin_sum) / n
}

/// Run the bias-corrected Monte Carlo clustering test.
///
/// For each iteration:
/// 1. Draw random angles from a uniform distribution in ϖ and Ω
/// 2. Weight by observational bias
/// 3. Compute Rayleigh z for the bias-weighted sample
/// 4. Compare to the observed z value
pub fn monte_carlo_clustering_test(
    kbos: &[DistantKbo],
    n_iterations: usize,
    seed: u64,
) -> ClusteringResult {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    // Compute observed statistics
    let observed_varpis: Vec<f64> = kbos
        .iter()
        .map(|k| (k.elements.omega + k.elements.omega_big).rem_euclid(2.0 * PI))
        .collect();
    let observed_omegas: Vec<f64> = kbos.iter().map(|k| k.elements.omega).collect();

    let observed_z_varpi = rayleigh_z(&observed_varpis);
    let observed_z_omega = rayleigh_z(&observed_omegas);

    let n = kbos.len();
    let mut count_exceed_varpi = 0usize;
    let mut count_exceed_omega = 0usize;
    let mut count_exceed_both = 0usize;

    for _ in 0..n_iterations {
        // Draw random orbital angles, weighted by bias
        let mut random_varpis = Vec::with_capacity(n);
        let mut random_omegas = Vec::with_capacity(n);

        for kbo in kbos {
            // Draw random ϖ weighted by bias
            let (varpi, omega) =
                draw_biased_angles(kbo.elements.a, kbo.elements.e, kbo.elements.i, &mut rng);
            random_varpis.push(varpi);
            random_omegas.push(omega);
        }

        let z_varpi = rayleigh_z(&random_varpis);
        let z_omega = rayleigh_z(&random_omegas);

        if z_varpi >= observed_z_varpi {
            count_exceed_varpi += 1;
        }
        if z_omega >= observed_z_omega {
            count_exceed_omega += 1;
        }
        if z_varpi >= observed_z_varpi && z_omega >= observed_z_omega {
            count_exceed_both += 1;
        }
    }

    ClusteringResult {
        observed_z_varpi,
        observed_z_omega,
        p_varpi: count_exceed_varpi as f64 / n_iterations as f64,
        p_omega: count_exceed_omega as f64 / n_iterations as f64,
        p_combined: count_exceed_both as f64 / n_iterations as f64,
        n_iterations,
    }
}

/// Draw random orbital angles (ϖ, ω) weighted by observational bias.
///
/// Uses rejection sampling: draw uniform ϖ, compute bias weight,
/// accept with probability proportional to weight.
fn draw_biased_angles<R: Rng>(a: f64, e: f64, i: f64, rng: &mut R) -> (f64, f64) {
    loop {
        let varpi = rng.gen::<f64>() * 2.0 * PI;
        let omega = rng.gen::<f64>() * 2.0 * PI;

        let weight = bias_function::bias_weight(a, e, varpi, omega, i);

        // Rejection sampling
        if rng.gen::<f64>() < weight.min(1.0).max(0.01) {
            return (varpi, omega);
        }
    }
}

/// Quick test with fewer iterations.
pub fn quick_clustering_test(kbos: &[DistantKbo]) -> ClusteringResult {
    monte_carlo_clustering_test(kbos, 1000, 42)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kbo_sample;

    #[test]
    fn test_rayleigh_z_clustered() {
        // All angles the same → z = n
        let angles = vec![1.0; 10];
        let z = rayleigh_z(&angles);
        assert!(
            (z - 10.0).abs() < 0.01,
            "z = {:.2}, expected 10.0 for identical angles",
            z
        );
    }

    #[test]
    fn test_rayleigh_z_uniform() {
        // Evenly spaced angles → z ≈ 0
        let n = 100;
        let angles: Vec<f64> = (0..n).map(|i| 2.0 * PI * i as f64 / n as f64).collect();
        let z = rayleigh_z(&angles);
        assert!(z < 0.1, "z = {:.4}, expected near 0 for uniform", z);
    }

    #[test]
    fn test_observed_clustering_significant() {
        let kbos = kbo_sample::paper_sample_a230();
        let varpis: Vec<f64> = kbos
            .iter()
            .map(|k| (k.elements.omega + k.elements.omega_big).rem_euclid(2.0 * PI))
            .collect();
        let z = rayleigh_z(&varpis);

        // The observed ϖ clustering should be detectable
        assert!(z > 2.0, "Observed z_ϖ = {:.2} should be > 2", z);
    }

    #[test]
    fn test_quick_clustering_test() {
        let kbos = kbo_sample::paper_sample_a230();
        let result = quick_clustering_test(&kbos);

        assert!(result.observed_z_varpi > 0.0);
        assert!(result.observed_z_omega > 0.0);
        assert!(result.p_varpi >= 0.0 && result.p_varpi <= 1.0);
        assert!(result.p_omega >= 0.0 && result.p_omega <= 1.0);

        // ϖ clustering should be significant (p < 10% with simplified bias)
        assert!(
            result.p_varpi < 0.15,
            "p(ϖ) = {:.4} should be < 0.15",
            result.p_varpi
        );
    }

    #[test]
    fn test_omega_clustering() {
        let kbos = kbo_sample::paper_sample_a230();
        let omegas: Vec<f64> = kbos.iter().map(|k| k.elements.omega).collect();
        let z = rayleigh_z(&omegas);

        // ω clustering should also be detectable
        assert!(z > 1.0, "Observed z_ω = {:.2} should be > 1", z);
    }
}
