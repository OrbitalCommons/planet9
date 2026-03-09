//! Clustering metrics for evaluating P9 parameter acceptability.
//!
//! Implements the statistical tests from Brown & Batygin (2016) Section 2:
//! - Perihelion confinement: fraction of survivors with anti-aligned Δϖ
//! - High-perihelion production: fraction of Sedna-like objects (q > 60 AU)
//! - Survival rate: number of particles remaining after 4 Gyr

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::types::OrbitalElements;

/// Evaluate clustering metrics for a set of final orbital elements.
///
/// `elements`: final orbital elements of surviving particles
/// `varpi_p9`: Planet Nine's longitude of perihelion (radians)
/// `a_min`, `a_max`: semi-major axis range for clustering analysis
///
/// Returns (clustering_fraction, high_perihelion_fraction, n_qualifying)
pub fn evaluate_clustering(
    elements: &[OrbitalElements],
    varpi_p9: f64,
    a_min: f64,
    a_max: f64,
) -> (f64, f64, usize) {
    let qualifying: Vec<&OrbitalElements> = elements
        .iter()
        .filter(|e| e.a >= a_min && e.a <= a_max && e.e < 1.0)
        .collect();

    let n = qualifying.len();
    if n == 0 {
        return (0.0, 0.0, 0);
    }

    let mut n_anti_aligned = 0;
    let mut n_high_q = 0;

    for elem in &qualifying {
        let varpi = elem.omega + elem.omega_big;
        let mut dv = varpi - varpi_p9;
        while dv > PI {
            dv -= TWO_PI;
        }
        while dv < -PI {
            dv += TWO_PI;
        }

        // Anti-aligned: |Δϖ| > π/2
        if dv.abs() > PI / 2.0 {
            n_anti_aligned += 1;
        }

        // High perihelion (Sedna-like): q > 60 AU
        let q = elem.a * (1.0 - elem.e);
        if q > 60.0 {
            n_high_q += 1;
        }
    }

    (
        n_anti_aligned as f64 / n as f64,
        n_high_q as f64 / n as f64,
        n,
    )
}

/// Rayleigh test for angular concentration.
///
/// The Rayleigh test statistic R² = (ΣsinΘ)² + (ΣcosΘ)² / n²
/// Under uniform distribution, 2n*R² ~ χ²(2).
///
/// Returns (R_bar, p_value) where R_bar is the mean resultant length
/// and p_value is the probability of this concentration occurring by chance.
pub fn rayleigh_test(angles: &[f64]) -> (f64, f64) {
    let n = angles.len() as f64;
    if n < 2.0 {
        return (0.0, 1.0);
    }

    let sin_sum: f64 = angles.iter().map(|&a| a.sin()).sum();
    let cos_sum: f64 = angles.iter().map(|&a| a.cos()).sum();

    let r_bar = (sin_sum * sin_sum + cos_sum * cos_sum).sqrt() / n;

    // Approximate p-value: p ≈ exp(-n * R_bar²)
    // This is exact for large n and R_bar not too close to 1
    let p_value = (-n * r_bar * r_bar).exp();

    (r_bar, p_value)
}

/// Compute confinement probability for a subsample.
///
/// From the paper: randomly draw `n_draw` objects from qualifying particles
/// (a ∈ [300,700], q < 80, i < 50°, survived > 3 Gyr) and check if they
/// show clustering similar to the observed 6 KBOs.
pub fn confinement_probability(
    elements: &[OrbitalElements],
    varpi_p9: f64,
    n_draw: usize,
    n_trials: usize,
    seed: u64,
) -> f64 {
    use rand::seq::SliceRandom;
    use rand::SeedableRng;

    let qualifying: Vec<&OrbitalElements> = elements
        .iter()
        .filter(|e| {
            e.a >= 300.0
                && e.a <= 700.0
                && e.perihelion() < 80.0
                && e.i < 50.0 * DEG2RAD
                && e.e < 1.0
        })
        .collect();

    if qualifying.len() < n_draw {
        return 0.0;
    }

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut n_confined = 0;

    for _ in 0..n_trials {
        let sample: Vec<&&OrbitalElements> = qualifying.choose_multiple(&mut rng, n_draw).collect();

        let varpis: Vec<f64> = sample
            .iter()
            .map(|e| {
                let varpi = e.omega + e.omega_big;
                let mut dv = varpi - varpi_p9;
                while dv > PI {
                    dv -= TWO_PI;
                }
                while dv < -PI {
                    dv += TWO_PI;
                }
                dv
            })
            .collect();

        let (r_bar, _) = rayleigh_test(&varpis);

        // Threshold for "confined": R_bar > 0.5
        if r_bar > 0.5 {
            n_confined += 1;
        }
    }

    n_confined as f64 / n_trials as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_clustering_evaluation() {
        // Create elements that are anti-aligned (ω + Ω ≈ π)
        let elements = vec![
            OrbitalElements {
                a: 400.0,
                e: 0.7,
                i: 0.0,
                omega: 3.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 500.0,
                e: 0.8,
                i: 0.0,
                omega: 3.2,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 350.0,
                e: 0.5,
                i: 0.0,
                omega: 0.5,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ];

        let varpi_p9 = 0.0;
        let (clust, high_q, n) = evaluate_clustering(&elements, varpi_p9, 300.0, 700.0);

        assert_eq!(n, 3);
        assert!(clust >= 0.0 && clust <= 1.0);
        assert!(high_q >= 0.0 && high_q <= 1.0);
    }

    #[test]
    fn test_rayleigh_uniform() {
        // Uniformly spaced angles should have low R_bar
        let angles: Vec<f64> = (0..20).map(|i| i as f64 / 20.0 * TWO_PI).collect();
        let (r_bar, p) = rayleigh_test(&angles);
        assert!(
            r_bar < 0.2,
            "Uniform angles should have low R_bar: {}",
            r_bar
        );
        assert!(p > 0.1, "Uniform angles should have high p-value: {}", p);
    }

    #[test]
    fn test_rayleigh_clustered() {
        // Tightly clustered angles should have high R_bar
        let angles: Vec<f64> = (0..20).map(|i| 1.0 + i as f64 * 0.01).collect();
        let (r_bar, p) = rayleigh_test(&angles);
        assert!(
            r_bar > 0.9,
            "Clustered angles should have high R_bar: {}",
            r_bar
        );
        assert!(p < 0.001, "Clustered angles should have low p-value: {}", p);
    }
}
