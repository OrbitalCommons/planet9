//! Von Mises distribution fitting for longitude of perihelion clustering.
//!
//! Single mode: f(theta | mu, kappa) = exp(kappa * cos(theta - mu)) / (2pi * I_0(kappa))
//! Bimodal: f = w * f_vm(mu1, kappa1) + (1-w) * f_vm(mu2, kappa2)
//!
//! Stable/metastable: varpi* ~ 50 deg, sigma ~ 0.5 rad
//! Unstable bimodal: w ~ 0.5, mu1 ~ 25 deg, mu2 ~ 205 deg

use p9_core::constants::DEG2RAD;
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

/// Fit a single von Mises distribution to a sample of angles via maximum likelihood.
///
/// Uses the circular mean for mu and an approximation for kappa.
pub fn fit_von_mises(angles: &[f64]) -> VonMisesParams {
    if angles.is_empty() {
        return VonMisesParams {
            mu: 0.0,
            kappa: 0.0,
        };
    }

    let n = angles.len() as f64;
    let sum_cos: f64 = angles.iter().map(|&a| a.cos()).sum();
    let sum_sin: f64 = angles.iter().map(|&a| a.sin()).sum();

    let mu = sum_sin.atan2(sum_cos);
    let r_bar = (sum_cos.powi(2) + sum_sin.powi(2)).sqrt() / n;

    // Approximate kappa from R-bar (Mardia & Jupp 2000)
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

    VonMisesParams {
        mu,
        kappa: kappa.max(0.0),
    }
}

/// Circular mean of a set of angles (rad).
pub fn circular_mean(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }
    let sum_cos: f64 = angles.iter().map(|&a| a.cos()).sum();
    let sum_sin: f64 = angles.iter().map(|&a| a.sin()).sum();
    sum_sin.atan2(sum_cos)
}

/// Circular standard deviation (rad).
pub fn circular_std(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }
    let n = angles.len() as f64;
    let sum_cos: f64 = angles.iter().map(|&a| a.cos()).sum();
    let sum_sin: f64 = angles.iter().map(|&a| a.sin()).sum();
    let r_bar = (sum_cos.powi(2) + sum_sin.powi(2)).sqrt() / n;
    (-2.0 * r_bar.ln()).sqrt().max(0.0)
}

/// Paper-reported clustering values for stable/metastable objects.
pub fn stable_clustering() -> VonMisesParams {
    VonMisesParams {
        mu: 50.0 * DEG2RAD,
        kappa: 1.0 / (0.5 * 0.5), // sigma ~ 0.5 rad
    }
}

/// Paper-reported bimodal fit for unstable objects.
pub fn unstable_clustering() -> BimodalVonMises {
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
    fn test_circular_mean_simple() {
        let angles = vec![0.0, 0.0, 0.0];
        let mean = circular_mean(&angles);
        assert!(mean.abs() < 1e-10);
    }

    #[test]
    fn test_circular_std_zero_for_identical() {
        let angles = vec![1.0, 1.0, 1.0];
        let std = circular_std(&angles);
        assert!(std < 0.01, "std = {} should be ~0", std);
    }

    #[test]
    fn test_stable_clustering_values() {
        let c = stable_clustering();
        let mu_deg = c.mu / DEG2RAD;
        assert!((mu_deg - 50.0).abs() < 1.0);
    }
}
