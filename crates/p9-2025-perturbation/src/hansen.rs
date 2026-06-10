//! Hansen coefficient computation for octupole and higher-order terms.
//!
//! X^{-4,1}_j and X^{-4,3}_j scale as j^{5/3} (steeper than the quadrupole
//! X^{-3,2}_j which scales linearly with j).
//!
//! Numerical computation via Mardling (2013) Eq. B19.

use serde::{Deserialize, Serialize};

/// Hansen coefficient X_j^{n,m}(e) via numerical quadrature.
///
/// X_j^{n,m}(e) = (1/2pi) integral_0^{2pi} (r/a)^n cos(m*f - j*M) dM
pub fn hansen_numerical(j: i64, n: i64, m: i64, e: f64, n_points: usize) -> f64 {
    let dm = std::f64::consts::TAU / n_points as f64;
    let mut sum = 0.0;

    for k in 0..n_points {
        let mean_anom = k as f64 * dm;
        let ecc_anom = kepler_solve(mean_anom, e);

        let cos_f = (ecc_anom.cos() - e) / (1.0 - e * ecc_anom.cos());
        let sin_f = (1.0 - e * e).sqrt() * ecc_anom.sin() / (1.0 - e * ecc_anom.cos());
        let f = sin_f.atan2(cos_f);

        let r_over_a = (1.0 - e * e) / (1.0 + e * f.cos());
        let integrand = r_over_a.powi(n as i32) * (m as f64 * f - j as f64 * mean_anom).cos();
        sum += integrand;
    }

    sum / n_points as f64
}

/// Approximate Hansen coefficient X_j^{-4,1} for octupole 1:j resonances.
///
/// Scales as ~ j^{5/3} * exp(-(q/a_N)^2) for large j.
pub fn hansen_neg4_1_approx(j: i64, q: f64, a_neptune: f64) -> f64 {
    let jf = j.abs() as f64;
    (jf.powf(5.0 / 3.0) / 8.0) * (-(q / a_neptune).powi(2)).exp()
}

/// Approximate Hansen coefficient X_j^{-4,3} for octupole 3:j resonances.
///
/// Scales as ~ j^{5/3} * exp(-(q/a_N)^2) for large j.
pub fn hansen_neg4_3_approx(j: i64, q: f64, a_neptune: f64) -> f64 {
    let jf = j.abs() as f64;
    (jf.powf(5.0 / 3.0) / 12.0) * (-(q / a_neptune).powi(2)).exp()
}

/// Solve Kepler's equation M = E - e*sin(E) via Newton iteration.
fn kepler_solve(m: f64, e: f64) -> f64 {
    let mut ecc_anom = m;
    for _ in 0..20 {
        let de = (ecc_anom - e * ecc_anom.sin() - m) / (1.0 - e * ecc_anom.cos());
        ecc_anom -= de;
        if de.abs() < 1e-14 {
            break;
        }
    }
    ecc_anom
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
    pub fn compute(e: f64, j_min: i64, j_max: i64, n_quadrature: usize) -> Self {
        let j_values: Vec<i64> = (j_min..=j_max).collect();
        let x_neg3_2: Vec<f64> = j_values
            .iter()
            .map(|&j| hansen_numerical(j, -3, 2, e, n_quadrature))
            .collect();
        let x_neg4_1: Vec<f64> = j_values
            .iter()
            .map(|&j| hansen_numerical(j, -4, 1, e, n_quadrature))
            .collect();
        let x_neg4_3: Vec<f64> = j_values
            .iter()
            .map(|&j| hansen_numerical(j, -4, 3, e, n_quadrature))
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

    #[test]
    fn test_hansen_numerical_j0_n_neg3_m0() {
        // X_0^{-3,0} at e=0 should be 1
        let val = hansen_numerical(0, -3, 0, 0.0, 512);
        assert!((val - 1.0).abs() < 0.01, "X_0^(-3,0) = {}", val);
    }

    #[test]
    fn test_hansen_neg4_1_positive() {
        for j in 10..30 {
            let val = hansen_neg4_1_approx(j, 35.0, 30.07);
            assert!(val > 0.0, "j={}: val = {}", j, val);
        }
    }

    #[test]
    fn test_hansen_neg4_scaling() {
        // X^{-4,1}_j should scale steeper than linearly (j^{5/3})
        let x10 = hansen_neg4_1_approx(10, 35.0, 30.07);
        let x20 = hansen_neg4_1_approx(20, 35.0, 30.07);
        let ratio = x20 / x10;
        let expected_ratio = 2.0_f64.powf(5.0 / 3.0);
        assert!(
            (ratio - expected_ratio).abs() < 0.1 * expected_ratio,
            "ratio = {}, expected = {}",
            ratio,
            expected_ratio
        );
    }

    #[test]
    fn test_kepler_solve_circular() {
        let e_anom = kepler_solve(1.0, 0.0);
        assert!((e_anom - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_hansen_table_dimensions() {
        let table = HansenTable::compute(0.7, 5, 15, 256);
        assert_eq!(table.j_values.len(), 11);
        assert_eq!(table.x_neg3_2.len(), 11);
        assert_eq!(table.x_neg4_1.len(), 11);
        assert_eq!(table.x_neg4_3.len(), 11);
    }
}
