//! Conservation law tracking for IOC dynamics.
//!
//! Semi-major axis a is conserved (orbit-averaged Hamiltonian).
//! Vertical angular momentum J_z = sqrt(1-e^2)*cos(i) is conserved by the
//! axisymmetry of the J2 + Miyamoto-Nagai terms (it is *not* conserved if
//! the optional galactic-tide term is enabled, which breaks axisymmetry).
//!
//! The primary dynamical invariant is the secular Hamiltonian itself —
//! see `check_hamiltonian_conservation` and the trajectory tests in `vzlk`.

use serde::{Deserialize, Serialize};

/// Vertical angular momentum (dimensionless).
///
/// J_z = sqrt(1 - e^2) * cos(i)
pub fn vertical_angular_momentum(e: f64, i: f64) -> f64 {
    (1.0 - e * e).sqrt() * i.cos()
}

/// Check conservation of a scalar series along a trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationCheck {
    pub initial: f64,
    pub final_value: f64,
    pub max_deviation: f64,
    pub relative_error: f64,
}

fn check_series(values: impl Iterator<Item = f64> + Clone) -> ConservationCheck {
    let mut iter = values.clone();
    let initial = match iter.next() {
        Some(v) => v,
        None => {
            return ConservationCheck {
                initial: 0.0,
                final_value: 0.0,
                max_deviation: 0.0,
                relative_error: 0.0,
            }
        }
    };

    let mut final_value = initial;
    let mut max_dev = 0.0_f64;
    for v in values {
        max_dev = max_dev.max((v - initial).abs());
        final_value = v;
    }

    let rel_err = if initial.abs() > 1e-300 {
        max_dev / initial.abs()
    } else {
        max_dev
    };

    ConservationCheck {
        initial,
        final_value,
        max_deviation: max_dev,
        relative_error: rel_err,
    }
}

/// Evaluate J_z conservation along an (e, i) trajectory.
pub fn check_conservation(e_series: &[f64], i_series: &[f64]) -> ConservationCheck {
    let jz = e_series
        .iter()
        .zip(i_series.iter())
        .map(|(&e, &i)| vertical_angular_momentum(e, i))
        .collect::<Vec<_>>();
    check_series(jz.iter().copied())
}

/// Evaluate Hamiltonian conservation along a trajectory's H series. This is
/// the primary invariant check for the secular flow (H generates it).
pub fn check_hamiltonian_conservation(h_series: &[f64]) -> ConservationCheck {
    check_series(h_series.iter().copied())
}

/// Standard Kozai-Lidov integral (quadrupole order, *external* perturber):
///
///   C_KL = e^2 * (1 - 5/2 * sin^2(i) * sin^2(omega))
///
/// Provided for diagnostic comparison only: it is exactly conserved by the
/// external-perturber quadrupole Hamiltonian, and only *approximately*
/// tracked here, where the perturbing potential is the extended (interior +
/// overlapping) Miyamoto-Nagai disk plus the planetary J2 ring. Use
/// `check_hamiltonian_conservation` for the rigorous invariant.
pub fn c_kl(e: f64, i: f64, omega: f64) -> f64 {
    e * e * (1.0 - 2.5 * i.sin().powi(2) * omega.sin().powi(2))
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
    fn test_c_kl_standard_values() {
        // Circular orbit: C_KL = 0 regardless of geometry.
        assert!(c_kl(0.0, 60.0 * DEG2RAD, 90.0 * DEG2RAD).abs() < 1e-15);
        // omega = 0: C_KL = e^2.
        let e = 0.4;
        assert!((c_kl(e, 50.0 * DEG2RAD, 0.0) - e * e).abs() < 1e-15);
        // i = 90 deg, omega = 90 deg: C_KL = e^2 (1 - 5/2) = -1.5 e^2.
        let c = c_kl(e, 90.0 * DEG2RAD, 90.0 * DEG2RAD);
        assert!((c - (-1.5 * e * e)).abs() < 1e-12, "C_KL = {c}");
    }

    #[test]
    fn test_hamiltonian_conservation_along_vzlk_trajectory() {
        use crate::hamiltonian::HamiltonianParams;
        use crate::vzlk::{evolutionary_timescale_gyr, integrate_vzlk};

        let params = HamiltonianParams {
            n_quadrature: 64,
            ..HamiltonianParams::default_paper()
        };
        let (a, j_z, e0, omega0) = (1000.0, 0.4, 0.7, 45.0 * DEG2RAD);
        let eta0 = (1.0_f64 - e0 * e0).sqrt();
        let tau = evolutionary_timescale_gyr(a, e0, (j_z / eta0).acos(), omega0, &params);
        let n_steps = 500;
        let dt = tau * p9_core::constants::GYR_DAYS / n_steps as f64;

        let traj = integrate_vzlk(a, j_z, e0, omega0, &params, dt, n_steps);
        let h: Vec<f64> = traj.iter().map(|p| p.h_value).collect();
        let check = check_hamiltonian_conservation(&h);
        assert!(
            check.relative_error < 1e-6,
            "relative H error = {:.3e}",
            check.relative_error
        );

        // J_z must also hold (axisymmetric default Hamiltonian).
        let e_series: Vec<f64> = traj.iter().map(|p| p.e).collect();
        let i_series: Vec<f64> = traj.iter().map(|p| p.i).collect();
        let jz_check = check_conservation(&e_series, &i_series);
        assert!(
            jz_check.relative_error < 1e-9,
            "relative J_z error = {:.3e}",
            jz_check.relative_error
        );
    }
}
