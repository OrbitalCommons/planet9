//! Semi-major axis diffusion estimation for scattered disk objects.
//!
//! `measure_diffusion` is a real diffusion estimator: it computes the
//! mean-squared displacement MSD(τ) averaged over *all time origins* and over
//! an *ensemble* of trajectories, then fits
//!
//!   MSD(τ) = D τ          (the paper's convention: D_a = ⟨Δa²⟩/t)
//!
//! over a range of lags, together with the log-log scaling exponent
//! MSD ∝ τ^γ. γ ≈ 1 indicates diffusive transport (the fit is meaningful);
//! γ ≈ 2 indicates ballistic drift, for which a "diffusion coefficient" is
//! not defined — the previous single-origin estimator could not tell these
//! apart and its own test asserted D > 0 for pure drift.

use serde::{Deserialize, Serialize};

use crate::stability::diffusion_coefficient;

/// Result of a diffusion measurement for a single configuration.
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

/// Fit of the mean-squared displacement curve.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DiffusionFit {
    /// Diffusion coefficient from the MSD(τ) = D·τ fit (AU²/yr)
    pub d: f64,
    /// Log-log scaling exponent γ in MSD ∝ τ^γ (1 ≈ diffusive, 2 ≈ drift)
    pub exponent: f64,
    /// Number of lags used in the fit
    pub n_lags: usize,
}

impl DiffusionFit {
    /// Whether the transport is consistent with diffusion (γ within
    /// `tol` of 1) rather than ballistic drift (γ ≈ 2).
    pub fn is_diffusive(&self, tol: f64) -> bool {
        (self.exponent - 1.0).abs() < tol
    }
}

/// Mean-squared displacement at lag `lag` (in samples), averaged over all
/// time origins of all trajectories in the ensemble.
pub fn msd_at_lag(ensemble: &[Vec<f64>], lag: usize) -> Option<f64> {
    let mut sum = 0.0;
    let mut count = 0usize;
    for series in ensemble {
        if series.len() <= lag {
            continue;
        }
        for t0 in 0..(series.len() - lag) {
            let d = series[t0 + lag] - series[t0];
            sum += d * d;
            count += 1;
        }
    }
    if count == 0 {
        None
    } else {
        Some(sum / count as f64)
    }
}

/// Estimate the diffusion coefficient from an ensemble of a(t) series with
/// uniform sampling interval `dt_yr` (years).
///
/// Lags span 1..=n_max/4 (capped at 64 log-spaced points). Returns `None`
/// when the ensemble is too short to support at least two lags.
pub fn measure_diffusion(ensemble: &[Vec<f64>], dt_yr: f64) -> Option<DiffusionFit> {
    let n_max = ensemble.iter().map(|s| s.len()).max().unwrap_or(0);
    if n_max < 8 || dt_yr <= 0.0 {
        return None;
    }
    let max_lag = (n_max / 4).max(2);

    // Up to 64 distinct lags, log-spaced between 1 and max_lag.
    let mut lags: Vec<usize> = Vec::new();
    let n_pts = 64usize;
    for k in 0..n_pts {
        let f = k as f64 / (n_pts - 1) as f64;
        let lag = ((max_lag as f64).powf(f)).round() as usize;
        let lag = lag.max(1);
        if lags.last() != Some(&lag) {
            lags.push(lag);
        }
    }

    let pairs: Vec<(f64, f64)> = lags
        .iter()
        .filter_map(|&lag| msd_at_lag(ensemble, lag).map(|msd| (lag as f64 * dt_yr, msd)))
        .filter(|&(_, msd)| msd > 0.0)
        .collect();
    if pairs.len() < 2 {
        return None;
    }

    // Least-squares slope through the origin: D = Σ τ·MSD / Σ τ².
    let num: f64 = pairs.iter().map(|&(t, m)| t * m).sum();
    let den: f64 = pairs.iter().map(|&(t, _)| t * t).sum();
    let d = num / den;

    // Log-log regression for the scaling exponent.
    let n = pairs.len() as f64;
    let lx: Vec<f64> = pairs.iter().map(|&(t, _)| t.ln()).collect();
    let ly: Vec<f64> = pairs.iter().map(|&(_, m)| m.ln()).collect();
    let mx = lx.iter().sum::<f64>() / n;
    let my = ly.iter().sum::<f64>() / n;
    let sxy: f64 = lx.iter().zip(&ly).map(|(&x, &y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = lx.iter().map(|&x| (x - mx) * (x - mx)).sum();
    let exponent = if sxx > 0.0 { sxy / sxx } else { f64::NAN };

    Some(DiffusionFit {
        d,
        exponent,
        n_lags: pairs.len(),
    })
}

/// Compare measured vs analytical diffusion coefficients across a sample.
pub fn diffusion_comparison(measurements: &[DiffusionMeasurement]) -> DiffusionComparisonResult {
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

/// Build a `DiffusionMeasurement` from a measured ensemble (e.g. from the
/// N-body validation in `nbody_validation`) against the analytic prediction.
pub fn measurement_from_ensemble(
    a: f64,
    q: f64,
    ensemble: &[Vec<f64>],
    dt_yr: f64,
) -> Option<DiffusionMeasurement> {
    let fit = measure_diffusion(ensemble, dt_yr)?;
    Some(DiffusionMeasurement {
        a_initial: a,
        q_initial: q,
        d_measured: fit.d,
        d_analytical: diffusion_coefficient(q),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{Rng, SeedableRng};

    #[test]
    fn test_measure_diffusion_constant() {
        // Constant semi-major axis => no positive-MSD lags => no fit, and
        // the MSD itself is identically zero.
        let series = vec![vec![100.0; 200]];
        assert_eq!(msd_at_lag(&series, 10), Some(0.0));
        assert!(measure_diffusion(&series, 1.0).is_none());
    }

    #[test]
    fn test_pure_drift_is_not_diffusive() {
        // Linear drift a(t) = 100 + 0.01 t: MSD ∝ τ², so the scaling
        // exponent flags ballistic transport — the old estimator asserted
        // D > 0 here, which is exactly the bug (M5).
        let series: Vec<Vec<f64>> = vec![(0..400).map(|i| 100.0 + 0.01 * i as f64).collect()];
        let fit = measure_diffusion(&series, 1.0).unwrap();
        assert!(
            (fit.exponent - 2.0).abs() < 0.05,
            "drift exponent = {:.3}",
            fit.exponent
        );
        assert!(!fit.is_diffusive(0.3));
    }

    #[test]
    fn test_random_walk_recovers_d() {
        // Unbiased random walk with per-step variance s²: MSD(τ) = s²·τ/dt,
        // so with dt = 1 yr the coefficient in MSD = D·τ is s².
        let mut rng = StdRng::seed_from_u64(7);
        let s = 0.5_f64;
        let ensemble: Vec<Vec<f64>> = (0..20)
            .map(|_| {
                let mut a = 200.0;
                (0..2000)
                    .map(|_| {
                        a += s * if rng.gen::<bool>() { 1.0 } else { -1.0 };
                        a
                    })
                    .collect()
            })
            .collect();
        let fit = measure_diffusion(&ensemble, 1.0).unwrap();
        assert!(
            (fit.exponent - 1.0).abs() < 0.15,
            "random-walk exponent = {:.3}",
            fit.exponent
        );
        assert!(fit.is_diffusive(0.3));
        let expected = s * s;
        assert!(
            (fit.d / expected - 1.0).abs() < 0.3,
            "D = {:.4}, expected {:.4}",
            fit.d,
            expected
        );
    }

    #[test]
    fn test_diffusion_comparison_empty() {
        let result = diffusion_comparison(&[]);
        assert_eq!(result.n_samples, 0);
    }

    #[test]
    fn test_diffusion_comparison_ratios() {
        let ms = vec![
            DiffusionMeasurement {
                a_initial: 200.0,
                q_initial: 35.0,
                d_measured: 2.0e-4,
                d_analytical: 1.0e-4,
            },
            DiffusionMeasurement {
                a_initial: 300.0,
                q_initial: 35.0,
                d_measured: 0.5e-4,
                d_analytical: 1.0e-4,
            },
        ];
        let r = diffusion_comparison(&ms);
        assert_eq!(r.n_samples, 2);
        assert!((r.mean_ratio - 1.25).abs() < 1e-12);
    }
}
