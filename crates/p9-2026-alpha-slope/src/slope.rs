//! The Loeb–Turner α-slope test: fit and classify the power-law index of flux
//! vs. heliocentric distance.
//!
//! Reflected sunlight reaches the observer after two inverse-square legs — the
//! Sun→body leg (incident flux ∝ d⁻²) and the body→observer leg (∝ Δ⁻²). Near
//! opposition Δ ≈ d, so the observed reflected flux falls as
//!
//!   F_refl ∝ d⁻⁴.
//!
//! Self-luminous (thermal or internal) emission travels only the body→observer
//! leg, so for a fixed source brightness it falls as
//!
//!   F_self ∝ d⁻².
//!
//! Fitting α in F ∝ dᵅ by an ordinary least-squares regression of log₁₀F on
//! log₁₀d (slope = α) is therefore a nature diagnostic (Loeb & Turner 2012):
//! α ≈ −4 ⇒ reflected sunlight (natural), α ≈ −2 ⇒ self-luminous (anomalous).
//!
//! Photometry usually arrives as apparent magnitudes; a magnitude maps to a
//! relative flux by f ∝ 10^(−0.4 m) (see [`magnitude_to_relative_flux`]).

use serde::{Deserialize, Serialize};

/// Theoretical α for reflected sunlight (two inverse-square legs): F ∝ d⁻⁴.
pub const ALPHA_REFLECTED: f64 = -4.0;

/// Theoretical α for self-luminous emission (one inverse-square leg): F ∝ d⁻².
pub const ALPHA_SELF_LUMINOUS: f64 = -2.0;

/// Half-width of the acceptance window (in α) around each theoretical value.
///
/// A fitted α within ±[`ALPHA_TOLERANCE`] of −4 is called reflected, within the
/// same window of −2 is called self-luminous, otherwise anomalous. The windows
/// (−4 ± 0.5 and −2 ± 0.5) do not overlap, so the classification is unambiguous.
pub const ALPHA_TOLERANCE: f64 = 0.5;

/// One photometric epoch: a heliocentric distance and the flux observed then.
///
/// Construct from a flux density with [`PhotometryPoint::from_flux`] or from an
/// apparent magnitude with [`PhotometryPoint::from_magnitude`].
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PhotometryPoint {
    /// Heliocentric distance d (AU).
    pub distance_au: f64,
    /// Flux (any consistent linear unit; only ratios across epochs matter).
    pub flux: f64,
}

impl PhotometryPoint {
    /// Epoch from a heliocentric distance (AU) and a flux density.
    pub fn from_flux(distance_au: f64, flux: f64) -> Self {
        Self { distance_au, flux }
    }

    /// Epoch from a heliocentric distance (AU) and an apparent magnitude,
    /// converting the magnitude to a relative flux f ∝ 10^(−0.4 m).
    pub fn from_magnitude(distance_au: f64, magnitude: f64) -> Self {
        Self {
            distance_au,
            flux: magnitude_to_relative_flux(magnitude),
        }
    }
}

/// Relative flux from an apparent magnitude: f = 10^(−0.4 m).
///
/// The zero point is arbitrary (it cancels in the log-flux regression), so the
/// returned value is a *relative* flux only.
pub fn magnitude_to_relative_flux(magnitude: f64) -> f64 {
    10f64.powf(-0.4 * magnitude)
}

/// Result of fitting F ∝ dᵅ to a set of photometry epochs.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SlopeFit {
    /// Fitted power-law index α (slope of log₁₀F vs log₁₀d).
    pub alpha: f64,
    /// 1σ standard error on α from the regression residuals.
    pub alpha_err: f64,
    /// Number of epochs used.
    pub n: usize,
    /// Coefficient of determination R² of the log–log fit.
    pub r_squared: f64,
}

/// Ordinary least-squares fit of α in F ∝ dᵅ.
///
/// Regresses y = log₁₀F on x = log₁₀d; the slope is α. Requires at least two
/// epochs spanning more than one distance (otherwise x has zero variance);
/// returns `None` if that fails or if any flux/distance is non-positive.
pub fn fit_alpha(points: &[PhotometryPoint]) -> Option<SlopeFit> {
    let n = points.len();
    if n < 2 {
        return None;
    }
    if points.iter().any(|p| p.flux <= 0.0 || p.distance_au <= 0.0) {
        return None;
    }

    let xs: Vec<f64> = points.iter().map(|p| p.distance_au.log10()).collect();
    let ys: Vec<f64> = points.iter().map(|p| p.flux.log10()).collect();

    let nf = n as f64;
    let x_mean = xs.iter().sum::<f64>() / nf;
    let y_mean = ys.iter().sum::<f64>() / nf;

    let mut sxx = 0.0;
    let mut sxy = 0.0;
    let mut syy = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let dx = x - x_mean;
        let dy = y - y_mean;
        sxx += dx * dx;
        sxy += dx * dy;
        syy += dy * dy;
    }
    if sxx <= 0.0 {
        return None;
    }

    let alpha = sxy / sxx;
    let intercept = y_mean - alpha * x_mean;

    // Residual sum of squares and the standard error on the slope.
    let mut ss_res = 0.0;
    for (&x, &y) in xs.iter().zip(ys.iter()) {
        let resid = y - (intercept + alpha * x);
        ss_res += resid * resid;
    }
    let alpha_err = if n > 2 {
        (ss_res / ((nf - 2.0) * sxx)).sqrt()
    } else {
        0.0
    };
    let r_squared = if syy > 0.0 { 1.0 - ss_res / syy } else { 1.0 };

    Some(SlopeFit {
        alpha,
        alpha_err,
        n,
        r_squared,
    })
}

/// Loeb–Turner classification of a fitted α.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum SlopeClass {
    /// α ≈ −4: reflected sunlight (natural). F ∝ d⁻⁴.
    Reflected,
    /// α ≈ −2: self-luminous / thermal (anomalous). F ∝ d⁻².
    SelfLuminous,
    /// α outside both acceptance windows.
    Anomalous,
}

/// Classify a fitted α into a [`SlopeClass`] using ±[`ALPHA_TOLERANCE`] windows
/// centered on the theoretical −4 (reflected) and −2 (self-luminous) values.
pub fn classify_alpha(alpha: f64) -> SlopeClass {
    if (alpha - ALPHA_REFLECTED).abs() <= ALPHA_TOLERANCE {
        SlopeClass::Reflected
    } else if (alpha - ALPHA_SELF_LUMINOUS).abs() <= ALPHA_TOLERANCE {
        SlopeClass::SelfLuminous
    } else {
        SlopeClass::Anomalous
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn magnitude_flux_round_trips() {
        // A 1-magnitude difference is a flux ratio of 10^0.4 ≈ 2.512.
        let f0 = magnitude_to_relative_flux(0.0);
        let f1 = magnitude_to_relative_flux(1.0);
        assert!((f0 / f1 - 10f64.powf(0.4)).abs() < 1e-9);
    }

    #[test]
    fn recovers_exact_power_law() {
        // Build flux exactly ∝ d^-4; the fit must recover α = -4 to machine
        // precision and R² = 1.
        let points: Vec<PhotometryPoint> = [30.0, 40.0, 50.0, 60.0, 80.0]
            .iter()
            .map(|&d| PhotometryPoint::from_flux(d, d.powf(-4.0)))
            .collect();
        let fit = fit_alpha(&points).expect("fit");
        assert!((fit.alpha - (-4.0)).abs() < 1e-9, "α = {}", fit.alpha);
        assert!(fit.r_squared > 1.0 - 1e-12);
    }

    #[test]
    fn classifier_buckets_each_regime() {
        assert_eq!(classify_alpha(-4.0), SlopeClass::Reflected);
        assert_eq!(classify_alpha(-3.7), SlopeClass::Reflected);
        assert_eq!(classify_alpha(-2.0), SlopeClass::SelfLuminous);
        assert_eq!(classify_alpha(-2.3), SlopeClass::SelfLuminous);
        assert_eq!(classify_alpha(-3.0), SlopeClass::Anomalous);
        assert_eq!(classify_alpha(-1.0), SlopeClass::Anomalous);
    }

    #[test]
    fn windows_do_not_overlap() {
        // Midpoint between -4 and -2 is -3, outside both windows.
        assert_eq!(classify_alpha(-3.0), SlopeClass::Anomalous);
    }

    #[test]
    fn single_point_fails() {
        let points = [PhotometryPoint::from_flux(40.0, 1.0)];
        assert!(fit_alpha(&points).is_none());
    }

    #[test]
    fn zero_variance_distance_fails() {
        let points = [
            PhotometryPoint::from_flux(40.0, 1.0),
            PhotometryPoint::from_flux(40.0, 2.0),
        ];
        assert!(fit_alpha(&points).is_none());
    }
}
