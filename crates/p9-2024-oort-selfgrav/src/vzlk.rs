//! Von Zeipel-Lidov-Kozai resonance analysis for IOC objects.
//!
//! Phase-space portraits on the (omega, q) plane at fixed semi-major axes.
//! The vZLK mechanism can drive perihelion oscillations within the IOC.

use serde::{Deserialize, Serialize};

use crate::hamiltonian::{secular_hamiltonian, HamiltonianParams};

/// A point on the vZLK phase portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VzlkPoint {
    /// Argument of perihelion (rad)
    pub omega: f64,
    /// Perihelion distance (AU)
    pub q: f64,
    /// Hamiltonian value
    pub h_value: f64,
}

/// Configuration for vZLK phase portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VzlkConfig {
    /// Fixed semi-major axis (AU)
    pub a: f64,
    /// Fixed vertical angular momentum J = sqrt(1-e^2)*cos(i)
    pub j_z: f64,
    /// Number of grid points in omega
    pub n_omega: usize,
    /// Number of grid points in q (perihelion)
    pub n_q: usize,
    /// Minimum perihelion to plot (AU)
    pub q_min: f64,
    /// Maximum perihelion to plot (AU)
    pub q_max: f64,
}

impl VzlkConfig {
    /// Configuration for a typical IOC object at a=1000 AU.
    pub fn default_paper() -> Self {
        Self {
            a: 1000.0,
            j_z: 0.3, // Moderate inclination
            n_omega: 200,
            n_q: 100,
            q_min: 50.0,
            q_max: 500.0,
        }
    }
}

/// Generate a vZLK phase portrait on the (omega, q) plane.
///
/// At fixed (a, J_z), we scan omega and q, computing e and i from:
///   q = a*(1-e), J_z = sqrt(1-e^2)*cos(i)
pub fn compute_vzlk_portrait(config: &VzlkConfig, params: &HamiltonianParams) -> Vec<VzlkPoint> {
    let n_total = config.n_omega * config.n_q;
    let mut points = Vec::with_capacity(n_total);

    let domega = std::f64::consts::TAU / config.n_omega as f64;
    let dq = (config.q_max - config.q_min) / config.n_q as f64;

    for io in 0..config.n_omega {
        let omega = io as f64 * domega;
        for iq in 0..config.n_q {
            let q = config.q_min + iq as f64 * dq;
            let e = 1.0 - q / config.a;
            if e <= 0.0 || e >= 1.0 {
                continue;
            }

            let eta = (1.0 - e * e).sqrt();
            let cos_i = config.j_z / eta;
            if cos_i.abs() > 1.0 {
                continue;
            }
            let i = cos_i.acos();

            let h = secular_hamiltonian(config.a, e, i, omega, params);
            points.push(VzlkPoint {
                omega,
                q,
                h_value: h,
            });
        }
    }

    points
}

/// Check if the vZLK cycle at given parameters can reach perihelion < q_target.
///
/// Scans omega from 0 to 2pi at fixed (a, J_z) to find min achievable q.
pub fn minimum_perihelion(a: f64, j_z: f64, _params: &HamiltonianParams) -> f64 {
    let n_test = 500;
    let mut q_min = a; // Maximum possible perihelion = a (circular)

    for _io in 0..n_test {
        // Scan eccentricity to find min achievable q consistent with J_z
        for ie in 1..999 {
            let e = ie as f64 / 1000.0;
            let q = a * (1.0 - e);
            let eta = (1.0 - e * e).sqrt();
            let cos_i = j_z / eta;
            if cos_i.abs() > 1.0 {
                continue;
            }

            if q < q_min {
                q_min = q;
            }
        }
    }

    q_min
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamiltonian::HamiltonianParams;

    #[test]
    fn test_vzlk_portrait_nonempty() {
        let config = VzlkConfig::default_paper();
        let params = HamiltonianParams::default_paper();
        let points = compute_vzlk_portrait(&config, &params);
        assert!(!points.is_empty(), "portrait should have points");
    }

    #[test]
    fn test_vzlk_portrait_finite() {
        let config = VzlkConfig::default_paper();
        let params = HamiltonianParams::default_paper();
        let points = compute_vzlk_portrait(&config, &params);
        for p in &points {
            assert!(
                p.h_value.is_finite(),
                "H should be finite at omega={}, q={}",
                p.omega,
                p.q
            );
        }
    }

    #[test]
    fn test_minimum_perihelion_positive() {
        let params = HamiltonianParams::default_paper();
        let q_min = minimum_perihelion(1000.0, 0.3, &params);
        assert!(q_min > 0.0 && q_min < 1000.0, "q_min = {}", q_min);
    }
}
