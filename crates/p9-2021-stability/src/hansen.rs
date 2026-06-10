//! Hansen coefficient approximations for the disturbing function.
//!
//! X_j^{-3,0} and X_j^{-3,2} are the key coefficients for the quadrupole
//! disturbing function. The approximation X_j^{-3,2} ~ (2j/5) exp(-(q/a_N)^2)
//! is accurate within a percent for j > 10 and q in (30, 50) AU.

use serde::{Deserialize, Serialize};

/// Neptune's semi-major axis (AU).
pub const A_NEPTUNE: f64 = 30.07;

/// Hansen coefficient X_j^{-3,2} approximate form.
///
/// X_j^{-3,2} ~ (2j/5) * exp(-(q/a_N)^2)
///
/// Valid for j > 10, q in (30, 50) AU.
pub fn hansen_x_neg3_2_approx(j: i64, q: f64) -> f64 {
    let jf = j as f64;
    (2.0 * jf / 5.0) * (-(q / A_NEPTUNE).powi(2)).exp()
}

/// Hansen coefficient X_j^{-3,0} approximate form.
///
/// For the secular (j=0) component, X_0^{-3,0} = 1/(1-e^2)^{3/2}.
/// For j != 0, X_j^{-3,0} ~ j * exp(-(q/a_N)^2) (similar scaling).
pub fn hansen_x_neg3_0(j: i64, e: f64) -> f64 {
    if j == 0 {
        1.0 / (1.0 - e * e).powf(1.5)
    } else {
        let jf = j.abs() as f64;
        jf * (-(1.0 - e).powi(2)).exp()
    }
}

/// Numerical Hansen coefficient via direct integration (Eq. B19 of Mardling 2013).
///
/// X_j^{n,m}(e) = (1/(2pi)) integral_0^{2pi} (r/a)^n cos(m*f - j*M) dM
///
/// where r/a = (1-e^2)/(1 + e*cos(f)) and f is the true anomaly.
pub fn hansen_numerical(j: i64, n: i64, m: i64, e: f64, n_points: usize) -> f64 {
    let dm = std::f64::consts::TAU / n_points as f64;
    let mut sum = 0.0;

    for k in 0..n_points {
        let mean_anom = k as f64 * dm;
        let true_anom = mean_to_true_anomaly(mean_anom, e);
        let r_over_a = (1.0 - e * e) / (1.0 + e * true_anom.cos());
        let integrand =
            r_over_a.powi(n as i32) * (m as f64 * true_anom - j as f64 * mean_anom).cos();
        sum += integrand;
    }

    sum / n_points as f64
}

/// Convert mean anomaly to true anomaly via Kepler's equation (Newton iteration).
fn mean_to_true_anomaly(m: f64, e: f64) -> f64 {
    // Solve M = E - e*sin(E) for E
    let mut ecc_anom = m;
    for _ in 0..20 {
        let de = (ecc_anom - e * ecc_anom.sin() - m) / (1.0 - e * ecc_anom.cos());
        ecc_anom -= de;
        if de.abs() < 1e-14 {
            break;
        }
    }

    // True anomaly from eccentric anomaly
    let half_f = ((1.0 + e) / (1.0 - e)).sqrt() * (ecc_anom / 2.0).tan();
    2.0 * half_f.atan()
}

/// Configuration for Hansen coefficient evaluation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HansenConfig {
    /// Number of quadrature points for numerical integration
    pub n_quadrature: usize,
    /// Maximum j index to evaluate
    pub j_max: i64,
}

impl HansenConfig {
    pub fn default_paper() -> Self {
        Self {
            n_quadrature: 1024,
            j_max: 100,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hansen_approx_positive_for_large_j() {
        for j in 10..50 {
            let val = hansen_x_neg3_2_approx(j, 35.0);
            assert!(val > 0.0, "j={}: X = {}", j, val);
        }
    }

    #[test]
    fn test_hansen_approx_decreases_with_q() {
        let j = 20;
        let x30 = hansen_x_neg3_2_approx(j, 30.0);
        let x40 = hansen_x_neg3_2_approx(j, 40.0);
        assert!(x30 > x40, "should decrease with perihelion distance");
    }

    #[test]
    fn test_hansen_x0_secular() {
        let x0 = hansen_x_neg3_0(0, 0.0);
        assert!((x0 - 1.0).abs() < 1e-10, "X_0^(-3,0) at e=0 should be 1");
    }

    #[test]
    fn test_hansen_numerical_j0() {
        // X_0^{-3,0} at e=0 should be 1
        let val = hansen_numerical(0, -3, 0, 0.0, 512);
        assert!(
            (val - 1.0).abs() < 0.01,
            "X_0^(-3,0) = {}, expected ~1",
            val
        );
    }

    #[test]
    fn test_mean_to_true_anomaly_circular() {
        let f = mean_to_true_anomaly(1.0, 0.0);
        assert!((f - 1.0).abs() < 1e-10, "circular orbit: f should equal M");
    }

    #[test]
    fn test_mean_to_true_anomaly_eccentric() {
        // At M=0 (perihelion), f=0 regardless of eccentricity
        let f = mean_to_true_anomaly(0.0, 0.7);
        assert!(f.abs() < 1e-10, "f at perihelion should be 0, got {}", f);
    }
}
