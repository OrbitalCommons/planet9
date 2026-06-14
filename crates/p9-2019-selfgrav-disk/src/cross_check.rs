//! Cross-validation of the analytic Laplace-Lagrange precession coefficient
//! against `p9_core::analysis::secular` numerical Gauss-ring averaging.
//!
//! The Laplace-Lagrange A coefficient (`secular::solve`) is the second-order
//! (in eccentricity) expansion of the same orbit-averaged disturbing function
//! that p9-core evaluates exactly by double anomaly averaging. For a single
//! perturbing ring at small test-particle eccentricity, the apsidal precession
//! rate dvarpi/dt = A must match the numerical curvature of the averaged
//! Hamiltonian in e, since to second order
//!
//!   <R> = n a^2 (1/2 A e^2 + ...)  and  dvarpi/dt = (1/(n a^2 e)) dR/de.
//!
//! We reuse p9-core's `numerical_secular_hamiltonian` (NOT re-implemented
//! here) and recover dvarpi/dt = (1/(n a^2)) d^2 H_grav/de^2 at e -> 0 with a
//! CIRCULAR ring (so only the A term survives; the forced/B term needs e_ring
//! != 0). This pins that our Laplace-Lagrange operator agrees with the exact
//! ring average to the expected second-order accuracy.

use crate::laplace::softened_laplace;
use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::GM_SUN;

/// Analytic single-ring A coefficient (radians/day): the test particle's
/// free precession rate induced by one circular perturbing ring.
pub fn single_ring_a_coeff(a: f64, a_ring: f64, mass_ring_solar: f64, eps: f64) -> f64 {
    let n = (GM_SUN / a.powi(3)).sqrt();
    let (alpha, abar) = if a < a_ring {
        let al = a / a_ring;
        (al, al)
    } else {
        (a_ring / a, 1.0)
    };
    let b1 = softened_laplace(1.5, 1, alpha, eps, 1024);
    0.25 * n * mass_ring_solar * alpha * abar * b1
}

/// Numerical precession rate from p9-core's exact ring average: the curvature
/// of the gravitational Hamiltonian in e at e -> 0 for a circular perturber.
///
///   dvarpi/dt = (1/(n a^2)) * d^2 <H_grav>/de^2 |_{e=0}
///
/// where <H_grav> = - gm_ring <<1/Delta>> is what p9-core computes (the secular
/// disturbing function is R = -<H_grav> minus constants; the e-curvature is
/// shared). The monopole/constant terms drop under the second derivative.
pub fn numerical_a_coeff(a: f64, a_ring: f64, mass_ring_solar: f64) -> f64 {
    let n = (GM_SUN / a.powi(3)).sqrt();
    let gm_ring = mass_ring_solar * GM_SUN;
    let nq = 200;
    // Softening matched to a small fraction of the larger semi-major axis so
    // the (non-crossing, small-e) average is smooth.
    let soft = 0.01 * a.max(a_ring);
    let h =
        |e: f64| numerical_secular_hamiltonian(a, e, 0.0, 0.0, 0.0, a_ring, 0.0, gm_ring, nq, soft);
    let de = 1e-3;
    // <H_grav> = -gm_ring<<1/Delta>>; the disturbing function R = -<H_grav>
    // (up to the indirect/constant terms). d^2 R/de^2 = -d^2 H/de^2.
    let d2h = (h(de) - 2.0 * h(0.0) + h(-de)) / (de * de);
    let d2r = -d2h;
    // R = n a^2 (1/2 A e^2) => d^2 R/de^2 = n a^2 A.
    d2r / (n * a * a)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_analytic_matches_numerical_ring_average() {
        // Single circular ring well separated from the test particle so the
        // second-order Laplace-Lagrange truncation is accurate. 1 Earth mass.
        let a = 250.0;
        let a_ring = 600.0;
        let m = p9_core::constants::EARTH_MASS_SOLAR;
        let analytic = single_ring_a_coeff(a, a_ring, m, 0.0);
        let numerical = numerical_a_coeff(a, a_ring, m);
        let rel = ((analytic - numerical) / numerical).abs();
        assert!(
            rel < 0.05,
            "Laplace-Lagrange A vs numerical ring average: analytic {analytic:.4e}, \
             numerical {numerical:.4e}, rel {rel:.3e}"
        );
    }

    #[test]
    fn test_analytic_matches_numerical_interior_ring() {
        // Test particle EXTERIOR to the ring (ring at 100 AU, particle 250 AU).
        let a = 250.0;
        let a_ring = 100.0;
        let m = p9_core::constants::EARTH_MASS_SOLAR;
        let analytic = single_ring_a_coeff(a, a_ring, m, 0.0);
        let numerical = numerical_a_coeff(a, a_ring, m);
        let rel = ((analytic - numerical) / numerical).abs();
        assert!(
            rel < 0.05,
            "interior-ring A mismatch: analytic {analytic:.4e}, numerical {numerical:.4e}, \
             rel {rel:.3e}"
        );
    }

    #[test]
    fn test_a_coeff_positive() {
        // The free precession rate from a coplanar ring is prograde (A > 0).
        let a = single_ring_a_coeff(250.0, 600.0, p9_core::constants::EARTH_MASS_SOLAR, 0.0);
        assert!(a > 0.0);
    }
}
