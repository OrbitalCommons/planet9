//! Semi-major axis diffusion measurement and analytical criterion.
//!
//! Measured: D_a from a least-squares fit of the mean-squared
//! displacement MSD(τ) = ⟨(a(t+τ) − a(t))²⟩_t over all time origins
//! (ensemble-of-origins average), D = MSD/τ slope. The previous estimator
//! took max over t of (Δa)²/t from a single origin, which overestimates —
//! the running supremum of a random walk grows like t·log log t.
//!
//! Analytical: D_a ~ (8/5π)(m_N/M_sun)·√(GM_sun·a_N)·exp(−q²/(2a_N²))
//! (Batygin et al. 2021). The exponential is the canonical
//! exp(−q²/2a_N²) of `p9_core::analysis::resonance` — D_a ∝ K, the
//! Chirikov overlap parameter, at fixed a (see the consistency test).

use p9_core::constants::{A_NEPTUNE_AU, GM_SUN, MASS_NEPTUNE_SOLAR, YEAR_DAYS};
use serde::{Deserialize, Serialize};

/// Analytical diffusion coefficient (AU^2/yr) from Batygin et al. (2021).
pub fn analytical_diffusion(q: f64) -> f64 {
    let v_n = (GM_SUN * A_NEPTUNE_AU).sqrt();
    let d_per_day = (8.0 / (5.0 * std::f64::consts::PI))
        * MASS_NEPTUNE_SOLAR
        * v_n
        * (-(q / A_NEPTUNE_AU).powi(2) / 2.0).exp();
    d_per_day * YEAR_DAYS
}

/// Measure the diffusion coefficient (AU²/yr) from a time series of
/// semi-major axis values sampled every `dt_yr` years.
///
/// A boxcar filter of `window_size` samples first removes short-period
/// oscillations; the MSD is then averaged over all time origins for lags
/// up to a quarter of the series, and D is the least-squares slope of
/// MSD(τ) = D·τ through the origin.
pub fn measure_diffusion(a_series: &[f64], dt_yr: f64, window_size: usize) -> f64 {
    if a_series.len() < 4 {
        return 0.0;
    }

    let smoothed = uniform_filter(a_series, window_size);
    let n = smoothed.len();
    let max_lag = (n / 4).max(1);

    let mut sum_msd_tau = 0.0;
    let mut sum_tau_sq = 0.0;
    for lag in 1..=max_lag {
        let tau = lag as f64 * dt_yr;
        let mut acc = 0.0;
        let mut count = 0usize;
        for t in 0..(n - lag) {
            let da = smoothed[t + lag] - smoothed[t];
            acc += da * da;
            count += 1;
        }
        let msd = acc / count as f64;
        sum_msd_tau += msd * tau;
        sum_tau_sq += tau * tau;
    }

    if sum_tau_sq <= 0.0 {
        return 0.0;
    }
    sum_msd_tau / sum_tau_sq
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

/// Median of a slice (averaging the two central elements for even n).
fn median(sorted: &[f64]) -> f64 {
    let n = sorted.len();
    if n % 2 == 1 {
        sorted[n / 2]
    } else {
        0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
    }
}

/// Median absolute deviation, scaled by 1.4826 to be a consistent
/// estimator of σ for Gaussian data. (The previous helper computed the
/// *mean* absolute deviation despite its name.)
pub fn median_absolute_deviation(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = median(&sorted);
    let mut devs: Vec<f64> = values.iter().map(|&v| (v - med).abs()).collect();
    devs.sort_by(|a, b| a.partial_cmp(b).unwrap());
    1.4826 * median(&devs)
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

/// Compute mean and std of diffusion measurements across clones, with
/// MAD-based outlier rejection (> 3 robust σ from the median).
pub fn aggregate_clone_diffusion(name: &'static str, d_values: &[f64], q: f64) -> DiffusionResult {
    if d_values.is_empty() {
        return DiffusionResult {
            name,
            d_mean: 0.0,
            d_std: 0.0,
            d_analytical: analytical_diffusion(q),
        };
    }

    let mut sorted = d_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = median(&sorted);
    let mad = median_absolute_deviation(d_values);

    let filtered: Vec<f64> = d_values
        .iter()
        .filter(|&&d| (d - med).abs() <= 3.0 * mad.max(1e-300))
        .copied()
        .collect();
    // Degenerate case (e.g. MAD = 0 with scattered values): keep all.
    let filtered = if filtered.is_empty() {
        d_values.to_vec()
    } else {
        filtered
    };

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
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal};

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
    fn test_analytical_diffusion_consistent_with_p9_core_chirikov() {
        // H3: at fixed a, D_a(q) ∝ K(a, q) — the same exp(−q²/2a_N²).
        use p9_core::analysis::resonance::chirikov_overlap_parameter;
        let a = 300.0;
        let ratio_d = analytical_diffusion(32.0) / analytical_diffusion(42.0);
        let ratio_k = chirikov_overlap_parameter(a, 32.0) / chirikov_overlap_parameter(a, 42.0);
        assert!(
            (ratio_d - ratio_k).abs() < 1e-9 * ratio_k,
            "D ratio {ratio_d} vs K ratio {ratio_k}"
        );
    }

    #[test]
    fn test_measure_diffusion_constant() {
        let series = vec![100.0; 100];
        let d = measure_diffusion(&series, 1e6, 5);
        assert!(d.abs() < 1e-10);
    }

    #[test]
    fn test_measure_diffusion_recovers_random_walk() {
        // Random walk with step σ per sample: true D = σ²/dt
        // (MSD(τ) = σ²·τ/dt). The MSD fit must land within ~2× without
        // the systematic high bias of the old max-over-t estimator.
        let mut rng = rand::rngs::StdRng::seed_from_u64(12345);
        let sigma = 0.05_f64;
        let dt_yr = 1000.0;
        let step = Normal::new(0.0, sigma).unwrap();
        let mut a = 150.0;
        let series: Vec<f64> = (0..4000)
            .map(|_| {
                a += step.sample(&mut rng);
                a
            })
            .collect();
        let d_true = sigma * sigma / dt_yr;
        let d = measure_diffusion(&series, dt_yr, 1);
        assert!(
            d > 0.4 * d_true && d < 2.5 * d_true,
            "D = {d:.3e}, true = {d_true:.3e}"
        );
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
    fn test_median_absolute_deviation() {
        // values: median = 3, |dev| = [2, 1, 0, 1, 2] → median 1 → 1.4826
        let values = [1.0, 2.0, 3.0, 4.0, 5.0];
        let mad = median_absolute_deviation(&values);
        assert!((mad - 1.4826).abs() < 1e-12, "MAD = {mad}");
        // For Gaussian draws, MAD ≈ σ.
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        let normal = Normal::new(10.0, 2.0).unwrap();
        let draws: Vec<f64> = (0..20_000).map(|_| normal.sample(&mut rng)).collect();
        let mad_g = median_absolute_deviation(&draws);
        assert!((mad_g - 2.0).abs() < 0.1, "Gaussian MAD = {mad_g}");
    }

    #[test]
    fn test_aggregate_clone_diffusion() {
        let d_values = vec![1e-4, 1.1e-4, 0.9e-4, 1.05e-4, 0.95e-4];
        let result = aggregate_clone_diffusion("test", &d_values, 35.0);
        assert!(result.d_mean > 0.0);
        assert!(result.d_std < result.d_mean);
    }

    #[test]
    fn test_aggregate_rejects_outlier() {
        let d_values = vec![1e-4, 1.1e-4, 0.9e-4, 1.05e-4, 0.95e-4, 5.0];
        let result = aggregate_clone_diffusion("test", &d_values, 35.0);
        assert!(
            result.d_mean < 2e-4,
            "outlier not rejected: mean = {}",
            result.d_mean
        );
    }
}
