//! Hansen coefficients for octupole and higher-order terms.
//!
//! X^{-4,1}_j and X^{-4,3}_j scale as j^{5/3} (steeper than the quadrupole
//! X^{-3,2}_j which scales linearly with j).
//!
//! Numerical evaluation delegates to the single workspace implementation
//! `p9_core::analysis::hansen::hansen_coefficient` (eccentric-anomaly
//! midpoint rule with automatic convergence control — the previous local
//! copy here was a verbatim twin with a fixed node count that aliases at
//! high eccentricity, and a fragile E₀ = M Newton solver).

use p9_core::analysis::hansen::hansen_coefficient;
use serde::{Deserialize, Serialize};

/// Relative tolerance used for the numerical Hansen coefficients.
const HANSEN_TOL: f64 = 1.0e-8;

/// Hansen coefficient X_j^{n,m}(e) via the p9-core convergence-controlled
/// quadrature:
///
/// X_j^{n,m}(e) = (1/2pi) integral_0^{2pi} (r/a)^n cos(m*f - j*M) dM
pub fn hansen_numerical(j: i64, n: i64, m: i64, e: f64) -> f64 {
    hansen_coefficient(j, n as i32, m as i32, e, HANSEN_TOL)
}

/// Approximate Hansen coefficient X_j^{-4,1} for octupole 1:j resonances
/// at perihelion distance `q` (AU):
///
///   X_j^{-4,1} ≈ 0.71 · j^{5/3} · exp(−(q/a_N)²)
///
/// The j^{5/3} scaling follows Belyakov & Batygin (2025); the prefactor is
/// calibrated against the numerical `hansen_numerical(j, -4, 1, e)` with
/// e = 1 − q/(a_N j^{2/3}) over j ∈ [10, 20], q ∈ [30, 40] AU (agreement
/// to ≲ 3%; the previous 1/8 prefactor was ~5.7× too small).
pub fn hansen_neg4_1_approx(j: i64, q: f64, a_neptune: f64) -> f64 {
    let jf = j.abs() as f64;
    0.71 * jf.powf(5.0 / 3.0) * (-(q / a_neptune).powi(2)).exp()
}

/// Approximate Hansen coefficient X_j^{-4,3} for octupole 3:j resonances
/// at perihelion distance `q` (AU):
///
///   X_j^{-4,3} ≈ j^{5/3}/9.5 · exp(−(q/a_N)²)
///
/// Prefactor calibrated against `hansen_numerical(j, -4, 3, e)` with
/// e = 1 − q/(a_N (j/3)^{2/3}) over j ∈ [10, 20], q ≈ 35 AU (agreement to
/// ≲ 15%; the previous 1/12 prefactor was ~25% low).
pub fn hansen_neg4_3_approx(j: i64, q: f64, a_neptune: f64) -> f64 {
    let jf = j.abs() as f64;
    (jf.powf(5.0 / 3.0) / 9.5) * (-(q / a_neptune).powi(2)).exp()
}

/// Compute a table of Hansen coefficients for a range of j values.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HansenTable {
    pub j_values: Vec<i64>,
    pub x_neg3_2: Vec<f64>,
    pub x_neg4_1: Vec<f64>,
    pub x_neg4_3: Vec<f64>,
}

impl HansenTable {
    /// Compute Hansen table for given eccentricity and j range.
    pub fn compute(e: f64, j_min: i64, j_max: i64) -> Self {
        let j_values: Vec<i64> = (j_min..=j_max).collect();
        let x_neg3_2: Vec<f64> = j_values
            .iter()
            .map(|&j| hansen_numerical(j, -3, 2, e))
            .collect();
        let x_neg4_1: Vec<f64> = j_values
            .iter()
            .map(|&j| hansen_numerical(j, -4, 1, e))
            .collect();
        let x_neg4_3: Vec<f64> = j_values
            .iter()
            .map(|&j| hansen_numerical(j, -4, 3, e))
            .collect();

        Self {
            j_values,
            x_neg3_2,
            x_neg4_1,
            x_neg4_3,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::A_NEPTUNE_AU;

    #[test]
    fn test_hansen_numerical_j0_n_neg3_m0() {
        // X_0^{-3,0}(e) = (1 - e²)^{-3/2}
        for &e in &[0.0, 0.5, 0.9] {
            let val = hansen_numerical(0, -3, 0, e);
            let exact = (1.0 - e * e).powf(-1.5);
            assert!(
                (val - exact).abs() < 1e-6 * exact,
                "X_0^(-3,0)({e}) = {val}, expected {exact}"
            );
        }
    }

    #[test]
    fn test_hansen_neg4_1_matches_numerical() {
        // Not a tautology: the approximation is checked against the
        // independent numerical quadrature at the resonance geometry
        // e = 1 − q/(a_N j^{2/3}).
        for &(j, q) in &[(10_i64, 35.0_f64), (20, 35.0), (15, 40.0)] {
            let a = A_NEPTUNE_AU * (j as f64).powf(2.0 / 3.0);
            let e = 1.0 - q / a;
            let num = hansen_numerical(j, -4, 1, e);
            let approx = hansen_neg4_1_approx(j, q, A_NEPTUNE_AU);
            assert!(
                (num - approx).abs() < 0.10 * num.abs(),
                "j={j}, q={q}: numerical {num:.3} vs approx {approx:.3}"
            );
        }
    }

    #[test]
    fn test_hansen_neg4_3_matches_numerical() {
        for &(j, q) in &[(10_i64, 35.0_f64), (20, 35.0)] {
            let a = A_NEPTUNE_AU * (j as f64 / 3.0).powf(2.0 / 3.0);
            let e = 1.0 - q / a;
            assert!(e > 0.0);
            let num = hansen_numerical(j, -4, 3, e);
            let approx = hansen_neg4_3_approx(j, q, A_NEPTUNE_AU);
            assert!(
                (num - approx).abs() < 0.20 * num.abs(),
                "j={j}, q={q}: numerical {num:.3} vs approx {approx:.3}"
            );
        }
    }

    #[test]
    fn test_hansen_neg4_scaling() {
        // X^{-4,1}_j should scale steeper than linearly (j^{5/3})
        let x10 = hansen_neg4_1_approx(10, 35.0, A_NEPTUNE_AU);
        let x20 = hansen_neg4_1_approx(20, 35.0, A_NEPTUNE_AU);
        let ratio = x20 / x10;
        let expected_ratio = 2.0_f64.powf(5.0 / 3.0);
        assert!(
            (ratio - expected_ratio).abs() < 0.1 * expected_ratio,
            "ratio = {ratio}, expected = {expected_ratio}"
        );
    }

    #[test]
    fn test_hansen_high_e_no_aliasing() {
        // The fixed-256-node local copy aliased here; the p9-core
        // quadrature must stay finite and converged.
        let v = hansen_numerical(40, -3, 2, 0.95);
        assert!(v.is_finite());
        assert!(v > 0.0, "X_40^(-3,2)(0.95) = {v}");
    }

    #[test]
    fn test_hansen_table_dimensions() {
        let table = HansenTable::compute(0.7, 5, 15);
        assert_eq!(table.j_values.len(), 11);
        assert_eq!(table.x_neg3_2.len(), 11);
        assert_eq!(table.x_neg4_1.len(), 11);
        assert_eq!(table.x_neg4_3.len(), 11);
    }
}
