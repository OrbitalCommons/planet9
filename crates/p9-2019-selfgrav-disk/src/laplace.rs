//! Softened Laplace coefficients b_s^{(j)}(alpha).
//!
//! The classical Laplace coefficient
//!
//!   b_s^{(j)}(alpha) = (1/pi) integral_0^{2pi}
//!                        cos(j psi) / (1 - 2 alpha cos psi + alpha^2)^s  dpsi
//!
//! enters the orbit-averaged (secular) disturbing function of two coplanar
//! orbits with semi-major-axis ratio alpha = a_inner / a_outer < 1
//! (Murray & Dermott 1999, Ch. 7). Sefilian & Touma (2019) model the massive
//! debris disk as a continuum of confocal eccentric rings whose secular
//! coupling to a test particle is built from b_{3/2}^{(1)} and b_{3/2}^{(2)}.
//!
//! The b_{3/2}^{(j)} integrals diverge as alpha -> 1 (rings of nearly equal
//! semi-major axis); a real disk has finite radial width so the kernel is
//! softened. We add a softening to the closest-approach term, replacing
//! (1 - 2 alpha cos psi + alpha^2) by (1 - 2 alpha cos psi + alpha^2 + eps^2),
//! which is the standard Gaussian/Plummer-style regularization used for
//! self-gravitating disks (Hahn 2003; Sefilian & Touma 2019 Sec. 2.2 use a
//! softened-disk kernel for exactly this reason). The integral is evaluated by
//! the (spectrally accurate, periodic) trapezoid rule.

use std::f64::consts::PI;

/// Softened Laplace coefficient b_s^{(j)}(alpha) with a dimensionless
/// softening `eps` added in quadrature to the denominator distance.
///
/// `eps = 0` recovers the classical Laplace coefficient (valid only for
/// alpha bounded away from 1). For a disk of finite radial extent a nonzero
/// softening represents the finite ring thickness / disk scale height.
///
/// `n_quad` trapezoid nodes; 512 is spectrally converged for the smooth
/// softened kernel.
pub fn softened_laplace(s: f64, j: i32, alpha: f64, eps: f64, n_quad: usize) -> f64 {
    let mut sum = 0.0;
    let dpsi = 2.0 * PI / n_quad as f64;
    let eps2 = eps * eps;
    for k in 0..n_quad {
        let psi = k as f64 * dpsi;
        let denom = 1.0 - 2.0 * alpha * psi.cos() + alpha * alpha + eps2;
        sum += (j as f64 * psi).cos() / denom.powf(s);
    }
    // (1/pi) * integral; trapezoid weight dpsi, periodic so endpoints fold in.
    sum * dpsi / PI
}

/// Classical (unsoftened) Laplace coefficient.
pub fn laplace(s: f64, j: i32, alpha: f64, n_quad: usize) -> f64 {
    softened_laplace(s, j, alpha, 0.0, n_quad)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_laplace_small_alpha_series() {
        // Series expansions (Murray & Dermott Eq. 6.69-6.70), small alpha:
        //   b_{1/2}^{(0)} ~ 2 + (1/2) alpha^2 + ...
        //   b_{1/2}^{(1)} ~ alpha + (3/8) alpha^3 + ...
        let alpha = 0.05;
        let b0 = laplace(0.5, 0, alpha, 1024);
        let b1 = laplace(0.5, 1, alpha, 1024);
        assert_relative_eq!(b0, 2.0 + 0.5 * alpha * alpha, epsilon = 1e-4);
        assert_relative_eq!(b1, alpha + 0.375 * alpha.powi(3), epsilon = 1e-4);
    }

    #[test]
    fn test_b32_series_small_alpha() {
        // b_{3/2}^{(1)} ~ 3 alpha + ... ; b_{3/2}^{(2)} ~ (15/4) alpha^2 + ...
        let alpha = 0.04;
        let b1 = laplace(1.5, 1, alpha, 1024);
        let b2 = laplace(1.5, 2, alpha, 1024);
        assert_relative_eq!(b1, 3.0 * alpha, epsilon = 2e-3);
        assert_relative_eq!(b2, 3.75 * alpha * alpha, epsilon = 2e-3);
    }

    #[test]
    fn test_softening_tames_divergence() {
        // Near alpha = 1 the bare b_{3/2} blows up; softening keeps it finite.
        let bare = laplace(1.5, 1, 0.98, 4096);
        let soft = softened_laplace(1.5, 1, 0.98, 0.05, 4096);
        assert!(soft.is_finite());
        assert!(
            soft < bare,
            "softening should reduce the peak: {soft} vs {bare}"
        );
    }

    #[test]
    fn test_quadrature_converged() {
        let a = softened_laplace(1.5, 1, 0.6, 0.02, 256);
        let b = softened_laplace(1.5, 1, 0.6, 0.02, 512);
        assert_relative_eq!(a, b, epsilon = 1e-9);
    }
}
