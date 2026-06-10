//! Semi-major axis diffusion analysis for scattered disk objects.
//!
//! Diffusion timescales and MEGNO-like chaos characterization.

use serde::{Deserialize, Serialize};

use crate::stability::diffusion_coefficient;

/// Result of a diffusion measurement for a single orbit.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiffusionMeasurement {
    /// Initial semi-major axis (AU)
    pub a_initial: f64,
    /// Initial perihelion distance (AU)
    pub q_initial: f64,
    /// Measured diffusion coefficient (AU^2/yr)
    pub d_measured: f64,
    /// Analytical prediction (AU^2/yr)
    pub d_analytical: f64,
}

/// Measure the diffusion coefficient from a time series of semi-major axis values.
///
/// D_a = <(a(t) - a(0))^2> / t
///
/// Uses windowed mean-squared displacement.
pub fn measure_diffusion(a_series: &[f64], dt_yr: f64) -> f64 {
    if a_series.len() < 2 {
        return 0.0;
    }

    let a0 = a_series[0];
    let n = a_series.len();
    let t_total = (n - 1) as f64 * dt_yr;

    if t_total < 1e-10 {
        return 0.0;
    }

    // Mean-squared displacement at final time
    let msd: f64 = a_series
        .iter()
        .skip(n / 2) // Use second half for better statistics
        .map(|&a| (a - a0).powi(2))
        .sum::<f64>()
        / (n / 2) as f64;

    msd / t_total
}

/// Compare measured vs analytical diffusion coefficients across a sample.
pub fn diffusion_comparison(measurements: &[DiffusionMeasurement]) -> DiffusionComparisonResult {
    if measurements.is_empty() {
        return DiffusionComparisonResult {
            n_samples: 0,
            mean_ratio: 0.0,
            std_ratio: 0.0,
        };
    }

    let ratios: Vec<f64> = measurements
        .iter()
        .filter(|m| m.d_analytical > 1e-30)
        .map(|m| m.d_measured / m.d_analytical)
        .collect();

    if ratios.is_empty() {
        return DiffusionComparisonResult {
            n_samples: 0,
            mean_ratio: 0.0,
            std_ratio: 0.0,
        };
    }

    let n = ratios.len() as f64;
    let mean = ratios.iter().sum::<f64>() / n;
    let var = ratios.iter().map(|&r| (r - mean).powi(2)).sum::<f64>() / n;

    DiffusionComparisonResult {
        n_samples: ratios.len(),
        mean_ratio: mean,
        std_ratio: var.sqrt(),
    }
}

/// Summary statistics for diffusion coefficient comparison.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiffusionComparisonResult {
    pub n_samples: usize,
    pub mean_ratio: f64,
    pub std_ratio: f64,
}

/// Generate a synthetic diffusion measurement from the analytical model.
pub fn synthetic_measurement(a: f64, q: f64) -> DiffusionMeasurement {
    let d_analytical = diffusion_coefficient(q);
    DiffusionMeasurement {
        a_initial: a,
        q_initial: q,
        d_measured: d_analytical, // Perfect agreement for analytical test
        d_analytical,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_measure_diffusion_constant() {
        // Constant semi-major axis => zero diffusion
        let series = vec![100.0; 100];
        let d = measure_diffusion(&series, 1.0);
        assert!(d.abs() < 1e-10, "D = {} for constant a", d);
    }

    #[test]
    fn test_measure_diffusion_linear_drift() {
        // Linear drift a(t) = 100 + 0.01*t
        let series: Vec<f64> = (0..100).map(|i| 100.0 + 0.01 * i as f64).collect();
        let d = measure_diffusion(&series, 1.0);
        assert!(d > 0.0, "D should be positive for drifting a");
    }

    #[test]
    fn test_synthetic_measurement() {
        let m = synthetic_measurement(200.0, 35.0);
        assert!((m.d_measured - m.d_analytical).abs() < 1e-20);
    }

    #[test]
    fn test_diffusion_comparison_perfect() {
        let measurements: Vec<DiffusionMeasurement> = (0..10)
            .map(|i| synthetic_measurement(100.0 + i as f64 * 50.0, 35.0))
            .collect();
        let result = diffusion_comparison(&measurements);
        assert_eq!(result.n_samples, 10);
        assert!(
            (result.mean_ratio - 1.0).abs() < 0.01,
            "mean ratio = {}",
            result.mean_ratio
        );
    }

    #[test]
    fn test_diffusion_comparison_empty() {
        let result = diffusion_comparison(&[]);
        assert_eq!(result.n_samples, 0);
    }
}
