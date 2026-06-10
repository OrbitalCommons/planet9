//! Semi-major axis diffusion measurement and analytical criterion.
//!
//! D_a ~ delta_a^2 / tau, where tau = 10-100 Myr depending on q.
//! Analytical: D_a ~ (8/5pi)(m_N/M_sun)*sqrt(GM_sun*a_N)*exp(-(q/a_N)^2/2)

use p9_core::constants::GM_SUN;
use serde::{Deserialize, Serialize};

/// Neptune parameters.
const A_NEPTUNE: f64 = 30.07;
const M_NEPTUNE_SOLAR: f64 = 5.15e-5;

/// Analytical diffusion coefficient (AU^2/yr) from Batygin et al. (2021).
pub fn analytical_diffusion(q: f64) -> f64 {
    let v_n = (GM_SUN * A_NEPTUNE).sqrt();
    let d_per_day = (8.0 / (5.0 * std::f64::consts::PI))
        * M_NEPTUNE_SOLAR
        * v_n
        * (-(q / A_NEPTUNE).powi(2) / 2.0).exp();
    d_per_day * 365.25
}

/// Measure diffusion from a time series of semi-major axis values.
///
/// Applies a uniform filter to remove short-term (~Myr) oscillations,
/// then computes D_a = max(<delta_a^2>) / tau.
pub fn measure_diffusion(a_series: &[f64], dt_yr: f64, window_size: usize) -> f64 {
    if a_series.len() < window_size + 1 {
        return 0.0;
    }

    // Smooth with uniform filter
    let smoothed = uniform_filter(a_series, window_size);

    // Compute maximum mean-squared displacement
    let a0 = smoothed[0];
    let mut max_d = 0.0_f64;

    for (i, &a) in smoothed.iter().enumerate().skip(1) {
        let t = i as f64 * dt_yr;
        if t < 1e-10 {
            continue;
        }
        let d = (a - a0).powi(2) / t;
        max_d = max_d.max(d);
    }

    max_d
}

/// Simple uniform (boxcar) filter for a 1D series.
fn uniform_filter(data: &[f64], window: usize) -> Vec<f64> {
    if window == 0 || data.len() < window {
        return data.to_vec();
    }

    let half_w = window / 2;
    let n = data.len();
    let mut filtered = Vec::with_capacity(n);

    for i in 0..n {
        let lo = if i >= half_w { i - half_w } else { 0 };
        let hi = (i + half_w + 1).min(n);
        let sum: f64 = data[lo..hi].iter().sum();
        filtered.push(sum / (hi - lo) as f64);
    }

    filtered
}

/// Diffusion measurement result for a single TNO (averaged over clones).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiffusionResult {
    pub name: &'static str,
    /// Mean diffusion coefficient across clones (AU^2/yr)
    pub d_mean: f64,
    /// Std of diffusion coefficient across clones
    pub d_std: f64,
    /// Analytical prediction
    pub d_analytical: f64,
}

/// Compute mean and std of diffusion measurements across clones.
pub fn aggregate_clone_diffusion(name: &'static str, d_values: &[f64], q: f64) -> DiffusionResult {
    if d_values.is_empty() {
        return DiffusionResult {
            name,
            d_mean: 0.0,
            d_std: 0.0,
            d_analytical: analytical_diffusion(q),
        };
    }

    // Remove outliers (> 3 sigma from median)
    let mut sorted = d_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = sorted[sorted.len() / 2];
    let mad: f64 = sorted.iter().map(|&d| (d - median).abs()).sum::<f64>() / sorted.len() as f64;

    let filtered: Vec<f64> = d_values
        .iter()
        .filter(|&&d| (d - median).abs() < 3.0 * mad.max(1e-10))
        .copied()
        .collect();

    let n = filtered.len() as f64;
    let mean = filtered.iter().sum::<f64>() / n;
    let var = filtered.iter().map(|&d| (d - mean).powi(2)).sum::<f64>() / n;

    DiffusionResult {
        name,
        d_mean: mean,
        d_std: var.sqrt(),
        d_analytical: analytical_diffusion(q),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analytical_diffusion_positive() {
        let d = analytical_diffusion(35.0);
        assert!(d > 0.0);
    }

    #[test]
    fn test_analytical_diffusion_decreases_with_q() {
        let d30 = analytical_diffusion(30.0);
        let d40 = analytical_diffusion(40.0);
        assert!(d30 > d40);
    }

    #[test]
    fn test_measure_diffusion_constant() {
        let series = vec![100.0; 100];
        let d = measure_diffusion(&series, 1e6, 5);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_uniform_filter_identity() {
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let filtered = uniform_filter(&data, 1);
        for (a, b) in data.iter().zip(filtered.iter()) {
            assert!((a - b).abs() < 1e-10);
        }
    }

    #[test]
    fn test_aggregate_clone_diffusion() {
        let d_values = vec![1e-4, 1.1e-4, 0.9e-4, 1.05e-4, 0.95e-4];
        let result = aggregate_clone_diffusion("test", &d_values, 35.0);
        assert!(result.d_mean > 0.0);
        assert!(result.d_std < result.d_mean);
    }
}
