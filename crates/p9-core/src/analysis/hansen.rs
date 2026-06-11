//! Hansen coefficients X_j^{n,m}(e) for the disturbing function.
//!
//! Single workspace implementation (p9-2021-stability and p9-2025-perturbation
//! previously carried verbatim twin copies with no convergence control and a
//! dimensionally wrong asymptotic).
//!
//! The numerical evaluator integrates in *eccentric* anomaly with the
//! Jacobian dM = (1 − e cos E) dE — the integrand is smooth and periodic in E,
//! so the uniform midpoint rule converges spectrally — and doubles the node
//! count until the requested tolerance is met (the required resolution grows
//! like j (1 − e)^{−3/2}).
//!
//! Reference: Hughes (1981); Mardling (2013) Appendix B

use crate::constants::{A_NEPTUNE_AU, TWO_PI};
use crate::types::solve_kepler;

/// Largest node count attempted before giving up on the tolerance.
const MAX_POINTS: usize = 1 << 18;

/// Numerical Hansen coefficient
///
///   X_j^{n,m}(e) = (1/2π) ∮ (r/a)^n cos(m f − j M) dM
///
/// evaluated with automatic convergence control: the quadrature node count
/// doubles from 64 until successive estimates agree to `tol` (relative, with
/// an absolute floor of 1e-13 so that coefficients that are exactly zero
/// converge), panicking in debug builds if 2^18 nodes do not suffice.
pub fn hansen_coefficient(j: i64, n: i32, m: i32, e: f64, tol: f64) -> f64 {
    assert!((0.0..1.0).contains(&e), "eccentricity out of range: {e}");

    let mut n_points = 64usize;
    let mut prev = hansen_fixed(j, n, m, e, n_points);
    while n_points < MAX_POINTS {
        n_points *= 2;
        let next = hansen_fixed(j, n, m, e, n_points);
        if (next - prev).abs() <= tol * next.abs() + 1e-13 {
            return next;
        }
        prev = next;
    }
    debug_assert!(
        false,
        "Hansen quadrature did not converge: j={j}, n={n}, m={m}, e={e}"
    );
    prev
}

/// Fixed-node evaluation on a uniform eccentric-anomaly midpoint grid.
fn hansen_fixed(j: i64, n: i32, m: i32, e: f64, n_points: usize) -> f64 {
    let de = TWO_PI / n_points as f64;
    let mut sum = 0.0;
    for k in 0..n_points {
        let ea = (k as f64 + 0.5) * de;
        let r_over_a = 1.0 - e * ea.cos();
        let mean_anom = ea - e * ea.sin();
        // True anomaly from eccentric anomaly
        let true_anom = 2.0 * (((1.0 + e) / (1.0 - e)).sqrt() * (ea / 2.0).tan()).atan();
        let integrand =
            r_over_a.powi(n) * (m as f64 * true_anom - j as f64 * mean_anom).cos() * r_over_a;
        sum += integrand;
    }
    sum * de / TWO_PI
}

/// Closed form X_0^{-3,0}(e) = (1 − e²)^{−3/2} (the secular orbit average of
/// (a/r)³).
pub fn hansen_x0_neg3_0(e: f64) -> f64 {
    (1.0 - e * e).powf(-1.5)
}

/// Asymptotic X_j^{-3,2} for a particle in the j:2 mean-motion resonance with
/// Neptune:
///
///   X_j^{-3,2} ≈ (2j/5) · exp(−(q/a_N)²)
///
/// where q = a(1 − e) is the particle's perihelion distance in AU and
/// a_N = 30.07 AU. Valid for j ≳ 10 and q ∈ (30, 50) AU. Note the length
/// scale in the exponent is *Neptune's* semi-major axis, not the particle's
/// (using q/a, i.e. 1 − e, is dimensionally inconsistent and off ~4x at
/// e = 0.9 with the opposite e-trend).
pub fn hansen_x_neg3_2_neptune_chain(j: i64, q_au: f64) -> f64 {
    let jf = j.abs() as f64;
    (2.0 * jf / 5.0) * (-(q_au / A_NEPTUNE_AU).powi(2)).exp()
}

/// Mean anomaly → true anomaly via the safeguarded Kepler solver.
pub fn mean_to_true_anomaly(mean_anomaly: f64, e: f64) -> f64 {
    let ea = solve_kepler(e, mean_anomaly);
    2.0 * (((1.0 + e) / (1.0 - e)).sqrt() * (ea / 2.0).tan()).atan()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_hansen_matches_closed_form_x0() {
        for &e in &[0.0, 0.3, 0.7, 0.9] {
            let num = hansen_coefficient(0, -3, 0, e, 1e-10);
            let exact = hansen_x0_neg3_0(e);
            assert_relative_eq!(num, exact, max_relative = 1e-8);
        }
    }

    #[test]
    fn test_hansen_x0_neg2_0_closed_form() {
        // X_0^{-2,0}(e) = (1 - e²)^{-1/2}
        for &e in &[0.0, 0.4, 0.8] {
            let num = hansen_coefficient(0, -2, 0, e, 1e-10);
            let exact = (1.0 - e * e).powf(-0.5);
            assert_relative_eq!(num, exact, max_relative = 1e-8);
        }
    }

    #[test]
    fn test_hansen_x0_1_0_closed_form() {
        // X_0^{1,0}(e) = 1 + e²/2 (orbit average of r/a)
        let num = hansen_coefficient(0, 1, 0, 0.6, 1e-10);
        assert_relative_eq!(num, 1.0 + 0.18, max_relative = 1e-8);
    }

    #[test]
    fn test_hansen_circular_kronecker() {
        // At e = 0, X_j^{n,m} = δ_{jm}.
        assert_relative_eq!(
            hansen_coefficient(2, -3, 2, 0.0, 1e-10),
            1.0,
            epsilon = 1e-10
        );
        assert!(hansen_coefficient(3, -3, 2, 0.0, 1e-10).abs() < 1e-10);
    }

    #[test]
    fn test_hansen_high_e_converges() {
        // High-e, large-j case that aliases on a fixed 256-point grid.
        let v = hansen_coefficient(40, -3, 2, 0.95, 1e-8);
        assert!(v.is_finite());
    }

    #[test]
    fn test_neptune_chain_asymptotic_trends() {
        // Decreasing in q, increasing in j, positive.
        let x = hansen_x_neg3_2_neptune_chain(20, 35.0);
        assert!(x > 0.0);
        assert!(hansen_x_neg3_2_neptune_chain(20, 30.0) > hansen_x_neg3_2_neptune_chain(20, 40.0));
        assert!(hansen_x_neg3_2_neptune_chain(30, 35.0) > x);
    }

    #[test]
    fn test_mean_to_true_anomaly() {
        // Circular: f = M. Perihelion: f = 0 for any e.
        assert_relative_eq!(mean_to_true_anomaly(1.0, 0.0), 1.0, epsilon = 1e-10);
        assert!(mean_to_true_anomaly(0.0, 0.7).abs() < 1e-10);
        // High-e quarter orbit stays finite and ahead of M.
        let f = mean_to_true_anomaly(0.5, 0.9);
        assert!(f > 0.5 && f < std::f64::consts::PI);
    }
}
