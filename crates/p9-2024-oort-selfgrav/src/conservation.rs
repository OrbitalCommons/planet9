//! Conservation law tracking for IOC dynamics.
//!
//! Semi-major axis a is conserved (no dissipation).
//! Vertical angular momentum J_z = sqrt(1-e^2)*cos(i) is conserved
//! by axisymmetry of the Miyamoto-Nagai potential.

use serde::{Deserialize, Serialize};

/// Vertical angular momentum (dimensionless).
///
/// J_z = sqrt(1 - e^2) * cos(i)
pub fn vertical_angular_momentum(e: f64, i: f64) -> f64 {
    (1.0 - e * e).sqrt() * i.cos()
}

/// Check conservation of J_z along a trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationCheck {
    pub j_z_initial: f64,
    pub j_z_final: f64,
    pub max_deviation: f64,
    pub relative_error: f64,
}

/// Evaluate conservation quality over a trajectory.
pub fn check_conservation(e_series: &[f64], i_series: &[f64]) -> ConservationCheck {
    if e_series.is_empty() || i_series.is_empty() {
        return ConservationCheck {
            j_z_initial: 0.0,
            j_z_final: 0.0,
            max_deviation: 0.0,
            relative_error: 0.0,
        };
    }

    let j_z_initial = vertical_angular_momentum(e_series[0], i_series[0]);
    let j_z_final = vertical_angular_momentum(*e_series.last().unwrap(), *i_series.last().unwrap());

    let mut max_dev = 0.0_f64;
    for (e, i) in e_series.iter().zip(i_series.iter()) {
        let j_z = vertical_angular_momentum(*e, *i);
        let dev = (j_z - j_z_initial).abs();
        max_dev = max_dev.max(dev);
    }

    let rel_err = if j_z_initial.abs() > 1e-15 {
        max_dev / j_z_initial.abs()
    } else {
        max_dev
    };

    ConservationCheck {
        j_z_initial,
        j_z_final,
        max_deviation: max_dev,
        relative_error: rel_err,
    }
}

/// Kozai-Lidov constant: H_KL = (1-e^2)*(2 - 5*sin^2(i)*sin^2(omega)) / 2.
///
/// This is an additional approximate integral for the quadrupole-order problem.
pub fn kozai_constant(e: f64, i: f64, omega: f64) -> f64 {
    let eta_sq = 1.0 - e * e;
    eta_sq * (2.0 - 5.0 * i.sin().powi(2) * omega.sin().powi(2)) / 2.0
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_jz_circular_orbit() {
        let jz = vertical_angular_momentum(0.0, 30.0 * DEG2RAD);
        assert!((jz - 30.0_f64.to_radians().cos()).abs() < 1e-10);
    }

    #[test]
    fn test_jz_range() {
        let jz = vertical_angular_momentum(0.5, 45.0 * DEG2RAD);
        assert!(jz > 0.0 && jz < 1.0);
    }

    #[test]
    fn test_conservation_perfect() {
        let e_series = vec![0.5, 0.5, 0.5];
        let i_series = vec![0.5, 0.5, 0.5];
        let check = check_conservation(&e_series, &i_series);
        assert!(check.max_deviation < 1e-14);
    }

    #[test]
    fn test_conservation_deviation() {
        let e_series = vec![0.5, 0.6, 0.7];
        let i_series = vec![0.5, 0.4, 0.3];
        let check = check_conservation(&e_series, &i_series);
        assert!(check.max_deviation > 0.0);
    }

    #[test]
    fn test_kozai_constant_finite() {
        let hk = kozai_constant(0.5, 45.0 * DEG2RAD, 90.0 * DEG2RAD);
        assert!(hk.is_finite());
    }
}
