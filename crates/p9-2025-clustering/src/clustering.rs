//! Von Mises fitting and bias-aware significance tests for longitude-of-
//! perihelion clustering.
//!
//! Single mode: f(θ | µ, κ) = exp(κ·cos(θ − µ)) / (2π·I₀(κ))
//! Bimodal: f = w·f_vm(µ₁, κ₁) + (1−w)·f_vm(µ₂, κ₂), fitted by EM.
//!
//! Circular summary statistics (mean, R̄, Rayleigh p) come from
//! `p9_core::analysis::circular`; the observation-bias null uses a
//! documented longitude-dependent selection weighting with a seeded
//! resampling Monte Carlo.

use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::constants::DEG2RAD;
use rand::Rng;
use serde::{Deserialize, Serialize};

/// Parameters of a single von Mises distribution.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VonMisesParams {
    /// Mean direction (rad)
    pub mu: f64,
    /// Concentration parameter (higher = more concentrated)
    pub kappa: f64,
}

/// Parameters of a bimodal von Mises mixture.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BimodalVonMises {
    /// Weight of first component (0 to 1)
    pub w: f64,
    /// First component
    pub comp1: VonMisesParams,
    /// Second component
    pub comp2: VonMisesParams,
}

/// Modified Bessel function I_0(x) approximation.
///
/// Uses polynomial approximation valid for all x >= 0.
pub fn bessel_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let t = (x / 3.75).powi(2);
        1.0 + t
            * (3.5156229
                + t * (3.0899424
                    + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))))
    } else {
        let t = 3.75 / ax;
        (ax.exp() / ax.sqrt())
            * (0.39894228
                + t * (0.01328592
                    + t * (0.00225319
                        + t * (-0.00157565
                            + t * (0.00916281
                                + t * (-0.02057706
                                    + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377))))))))
    }
}

/// Von Mises probability density function.
pub fn von_mises_pdf(theta: f64, mu: f64, kappa: f64) -> f64 {
    (kappa * (theta - mu).cos()).exp() / (std::f64::consts::TAU * bessel_i0(kappa))
}

/// Bimodal von Mises mixture PDF.
pub fn bimodal_pdf(theta: f64, params: &BimodalVonMises) -> f64 {
    params.w * von_mises_pdf(theta, params.comp1.mu, params.comp1.kappa)
        + (1.0 - params.w) * von_mises_pdf(theta, params.comp2.mu, params.comp2.kappa)
}

/// Maximum-likelihood κ from the mean resultant length R̄
/// (Mardia & Jupp 2000 piecewise approximation of A⁻¹).
pub fn kappa_from_r_bar(r_bar: f64) -> f64 {
    let kappa = if r_bar < 0.53 {
        2.0 * r_bar + r_bar.powi(3) + 5.0 * r_bar.powi(5) / 6.0
    } else if r_bar < 0.85 {
        -0.4 + 1.39 * r_bar + 0.43 / (1.0 - r_bar)
    } else {
        let denom = r_bar.powi(3) - 4.0 * r_bar.powi(2) + 3.0 * r_bar;
        if denom.abs() < 1e-10 {
            1000.0 // Very concentrated distribution
        } else {
            1.0 / denom
        }
    };
    kappa.clamp(0.0, 1000.0)
}

/// Fit a single von Mises distribution to a sample of angles via maximum
/// likelihood (circular mean for µ, A⁻¹(R̄) for κ).
pub fn fit_von_mises(angles: &[f64]) -> VonMisesParams {
    if angles.is_empty() {
        return VonMisesParams {
            mu: 0.0,
            kappa: 0.0,
        };
    }
    let mu = circular_mean(angles).unwrap_or(0.0);
    let r_bar = mean_resultant_length(angles);
    VonMisesParams {
        mu,
        kappa: kappa_from_r_bar(r_bar),
    }
}

/// Fit a two-component von Mises mixture by expectation-maximization —
/// an actual fit (the previous "bimodal fit" only evaluated a hard-coded
/// PDF). Initialized from the sample circular mean and its antipode;
/// `n_iter` EM sweeps (50 is plenty for n ≲ 10³).
pub fn fit_bimodal_von_mises(angles: &[f64], n_iter: usize) -> BimodalVonMises {
    assert!(angles.len() >= 4, "need at least 4 angles for a 2-mode fit");

    let mu0 = circular_mean(angles).unwrap_or(0.0);
    let mut fit = BimodalVonMises {
        w: 0.5,
        comp1: VonMisesParams {
            mu: mu0,
            kappa: 1.0,
        },
        comp2: VonMisesParams {
            mu: mu0 + std::f64::consts::PI,
            kappa: 1.0,
        },
    };

    let n = angles.len() as f64;
    for _ in 0..n_iter {
        // E-step: responsibilities of component 1.
        let resp: Vec<f64> = angles
            .iter()
            .map(|&th| {
                let p1 = fit.w * von_mises_pdf(th, fit.comp1.mu, fit.comp1.kappa);
                let p2 = (1.0 - fit.w) * von_mises_pdf(th, fit.comp2.mu, fit.comp2.kappa);
                p1 / (p1 + p2).max(1e-300)
            })
            .collect();

        // M-step: weighted circular means and concentrations.
        let update = |weights: &dyn Fn(usize) -> f64| -> VonMisesParams {
            let mut sum_w = 0.0;
            let (mut sc, mut ss) = (0.0, 0.0);
            for (k, &th) in angles.iter().enumerate() {
                let w = weights(k);
                sum_w += w;
                sc += w * th.cos();
                ss += w * th.sin();
            }
            if sum_w < 1e-12 {
                return VonMisesParams {
                    mu: 0.0,
                    kappa: 0.0,
                };
            }
            let mu = ss.atan2(sc);
            let r_bar = (sc * sc + ss * ss).sqrt() / sum_w;
            VonMisesParams {
                mu,
                kappa: kappa_from_r_bar(r_bar).min(50.0),
            }
        };

        fit.comp1 = update(&|k| resp[k]);
        fit.comp2 = update(&|k| 1.0 - resp[k]);
        fit.w = (resp.iter().sum::<f64>() / n).clamp(0.01, 0.99);
    }
    fit
}

/// Log-likelihood of a bimodal fit (for model comparison).
pub fn bimodal_log_likelihood(angles: &[f64], params: &BimodalVonMises) -> f64 {
    angles
        .iter()
        .map(|&th| bimodal_pdf(th, params).max(1e-300).ln())
        .sum()
}

/// Documented longitude-dependent selection weighting for the
/// bias-resampling null.
///
/// Wide-field TNO surveys preferentially cover fields away from the
/// galactic plane and near the ecliptic opposition seasons, producing a
/// smooth modulation of the discovery probability with longitude of
/// perihelion (high-e objects are discovered near perihelion, so the
/// selection acts on ϖ). We model it as
///
///   w(ϖ) ∝ 1 + A·cos(ϖ − ϖ₀),  A = 0.3, ϖ₀ = 60°,
///
/// a first-harmonic stand-in with the amplitude and phase of the
/// Brown (2017) / Shankman et al. (2017) survey-bias estimates. The exact
/// shape is configurable via `bias_resampled_p_value`'s `weight`.
pub fn default_selection_weight(varpi: f64) -> f64 {
    1.0 + 0.3 * (varpi - 60.0 * DEG2RAD).cos()
}

/// Seeded bias-resampling Monte Carlo p-value for clustering.
///
/// Null hypothesis: the n observed longitudes are drawn independently
/// from the survey selection weighting `weight` (NOT from uniformity —
/// the inadequate null the paper rejects). The statistic is the mean
/// resultant length R̄; p = fraction of synthetic samples with R̄ at
/// least as large as observed.
pub fn bias_resampled_p_value<R: Rng>(
    angles: &[f64],
    weight: &dyn Fn(f64) -> f64,
    n_mc: usize,
    rng: &mut R,
) -> f64 {
    assert!(!angles.is_empty());
    let r_obs = mean_resultant_length(angles);
    let n = angles.len();

    // Rejection-sampling envelope.
    let w_max = (0..720)
        .map(|k| weight(k as f64 * std::f64::consts::TAU / 720.0))
        .fold(f64::MIN, f64::max);

    let mut count = 0usize;
    let mut sample = vec![0.0_f64; n];
    for _ in 0..n_mc {
        for s in sample.iter_mut() {
            *s = loop {
                let th = rng.gen::<f64>() * std::f64::consts::TAU;
                if rng.gen::<f64>() * w_max <= weight(th) {
                    break th;
                }
            };
        }
        if mean_resultant_length(&sample) >= r_obs {
            count += 1;
        }
    }
    // Add-one smoothing keeps p > 0 with finite Monte Carlo samples.
    (count + 1) as f64 / (n_mc + 1) as f64
}

/// Paper-reported clustering values for stable/metastable objects
/// (Pichierri & Batygin 2025) — reference numbers for regression
/// comparison only; the pipeline computes its own fits.
pub fn stable_clustering_paper() -> VonMisesParams {
    VonMisesParams {
        mu: 50.0 * DEG2RAD,
        kappa: 1.0 / (0.5 * 0.5), // sigma ~ 0.5 rad
    }
}

/// Paper-reported bimodal fit for unstable objects — reference numbers
/// for regression comparison only.
pub fn unstable_clustering_paper() -> BimodalVonMises {
    BimodalVonMises {
        w: 0.5,
        comp1: VonMisesParams {
            mu: 25.0 * DEG2RAD,
            kappa: 2.0,
        },
        comp2: VonMisesParams {
            mu: 205.0 * DEG2RAD,
            kappa: 1.0,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal};

    #[test]
    fn test_bessel_i0_at_zero() {
        assert!((bessel_i0(0.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_i0_positive() {
        for x in [0.5, 1.0, 2.0, 5.0, 10.0] {
            assert!(bessel_i0(x) > 0.0, "I_0({}) should be positive", x);
        }
    }

    #[test]
    fn test_von_mises_pdf_normalizes() {
        // Numerical integration should be ~1
        let n = 1000;
        let dt = std::f64::consts::TAU / n as f64;
        let integral: f64 = (0..n)
            .map(|i| von_mises_pdf(i as f64 * dt, 0.0, 2.0) * dt)
            .sum();
        assert!(
            (integral - 1.0).abs() < 0.01,
            "von Mises integral = {}",
            integral
        );
    }

    #[test]
    fn test_fit_von_mises_concentrated() {
        let angles: Vec<f64> = (0..100).map(|_| 1.0).collect();
        let fit = fit_von_mises(&angles);
        assert!((fit.mu - 1.0).abs() < 0.01, "mu = {}", fit.mu);
        assert!(fit.kappa > 10.0, "kappa = {} should be large", fit.kappa);
    }

    #[test]
    fn test_em_recovers_bimodal_mixture() {
        // Seeded sample from two wrapped-normal modes at 0.5 and 3.5 rad
        // (σ = 0.45 ↔ κ ≈ 5), weights 0.5/0.5. EM must find both modes.
        let mut rng = rand::rngs::StdRng::seed_from_u64(2025);
        let noise = Normal::new(0.0, 0.45).unwrap();
        let mut angles = Vec::new();
        for k in 0..400 {
            let mu: f64 = if k % 2 == 0 { 0.5 } else { 3.5 };
            angles.push((mu + noise.sample(&mut rng)).rem_euclid(std::f64::consts::TAU));
        }
        let fit = fit_bimodal_von_mises(&angles, 100);

        // Mode locations recovered (in either component order).
        let mut mus = [
            fit.comp1.mu.rem_euclid(std::f64::consts::TAU),
            fit.comp2.mu.rem_euclid(std::f64::consts::TAU),
        ];
        mus.sort_by(|a, b| a.partial_cmp(b).unwrap());
        assert!((mus[0] - 0.5).abs() < 0.2, "mode 1 = {:.2}", mus[0]);
        assert!((mus[1] - 3.5).abs() < 0.2, "mode 2 = {:.2}", mus[1]);
        assert!((fit.w - 0.5).abs() < 0.15, "w = {:.2}", fit.w);
        assert!(fit.comp1.kappa > 1.5 && fit.comp2.kappa > 1.5);

        // The mixture must beat a single-mode fit in log-likelihood.
        let single = fit_von_mises(&angles);
        let single_as_mixture = BimodalVonMises {
            w: 0.999,
            comp1: single.clone(),
            comp2: VonMisesParams {
                mu: 0.0,
                kappa: 0.0,
            },
        };
        assert!(
            bimodal_log_likelihood(&angles, &fit)
                > bimodal_log_likelihood(&angles, &single_as_mixture)
        );
    }

    #[test]
    fn test_bias_resampled_p_value_detects_clustering() {
        // Tightly clustered sample → small p even against the biased null.
        let mut rng = rand::rngs::StdRng::seed_from_u64(77);
        let noise = Normal::new(0.0, 0.3).unwrap();
        let angles: Vec<f64> = (0..10)
            .map(|_| (0.9_f64 + noise.sample(&mut rng)).rem_euclid(std::f64::consts::TAU))
            .collect();
        let p = bias_resampled_p_value(&angles, &default_selection_weight, 500, &mut rng);
        assert!(p < 0.05, "p = {p}");
    }

    #[test]
    fn test_bias_resampled_p_value_null_calibrated() {
        // A sample drawn from the weighting itself is NOT significant.
        let mut rng = rand::rngs::StdRng::seed_from_u64(78);
        let w_max = 1.3;
        let draw = |rng: &mut rand::rngs::StdRng| loop {
            let th = rng.gen::<f64>() * std::f64::consts::TAU;
            if rng.gen::<f64>() * w_max <= default_selection_weight(th) {
                break th;
            }
        };
        let angles: Vec<f64> = (0..10).map(|_| draw(&mut rng)).collect();
        let p = bias_resampled_p_value(&angles, &default_selection_weight, 500, &mut rng);
        assert!(p > 0.05, "p = {p} — null sample flagged as clustered");
    }

    #[test]
    fn test_kappa_from_r_bar_monotone() {
        assert!(kappa_from_r_bar(0.2) < kappa_from_r_bar(0.6));
        assert!(kappa_from_r_bar(0.6) < kappa_from_r_bar(0.9));
        assert!(kappa_from_r_bar(0.0) < 1e-9);
    }

    #[test]
    fn test_stable_clustering_paper_values() {
        let c = stable_clustering_paper();
        let mu_deg = c.mu / DEG2RAD;
        assert!((mu_deg - 50.0).abs() < 1.0);
    }
}
