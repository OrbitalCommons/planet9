//! Phase-space portrait generation for e-Delta_varpi dynamics.
//!
//! Generates Hamiltonian level curves showing the three-lobed structure
//! for high-inclination objects and the anti-aligned apsidal confinement
//! for low-inclination scattered disk objects.

use p9_core::constants::DEG2RAD;
use serde::{Deserialize, Serialize};

use crate::hamiltonian::{secular_hamiltonian, SecularHamiltonianParams};

/// A point on a phase-space portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhasePoint {
    /// Delta varpi = varpi - varpi_9 (rad)
    pub delta_varpi: f64,
    /// Eccentricity
    pub e: f64,
    /// Hamiltonian value at this point
    pub h_value: f64,
}

/// Configuration for phase portrait generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhasePortraitConfig {
    /// Semi-major axis (AU) of the test particle
    pub a: f64,
    /// Inclination (rad) of the test particle
    pub i: f64,
    /// Number of grid points in delta_varpi
    pub n_varpi: usize,
    /// Number of grid points in eccentricity
    pub n_e: usize,
    /// Minimum eccentricity
    pub e_min: f64,
    /// Maximum eccentricity
    pub e_max: f64,
}

impl PhasePortraitConfig {
    /// Portrait for low-inclination test particle at a=300 AU.
    pub fn low_inclination() -> Self {
        Self {
            a: 300.0,
            i: 10.0 * DEG2RAD,
            n_varpi: 200,
            n_e: 100,
            e_min: 0.3,
            e_max: 0.95,
        }
    }

    /// Portrait for high-inclination test particle at a=300 AU.
    pub fn high_inclination() -> Self {
        Self {
            a: 300.0,
            i: 60.0 * DEG2RAD,
            n_varpi: 200,
            n_e: 100,
            e_min: 0.3,
            e_max: 0.95,
        }
    }
}

/// Generate a grid of Hamiltonian values over the (delta_varpi, e) plane.
pub fn compute_phase_portrait(
    config: &PhasePortraitConfig,
    params: &SecularHamiltonianParams,
) -> Vec<PhasePoint> {
    let n_total = config.n_varpi * config.n_e;
    let mut points = Vec::with_capacity(n_total);
    let dv_step = std::f64::consts::TAU / config.n_varpi as f64;
    let e_step = (config.e_max - config.e_min) / config.n_e as f64;

    for iv in 0..config.n_varpi {
        let dv = iv as f64 * dv_step;
        for ie in 0..config.n_e {
            let e = config.e_min + ie as f64 * e_step;
            let h = secular_hamiltonian(config.a, e, config.i, 0.0, dv, params);
            points.push(PhasePoint {
                delta_varpi: dv,
                e,
                h_value: h,
            });
        }
    }

    points
}

/// Find the equilibrium points (fixed points) in the phase portrait.
///
/// These are local extrema of the Hamiltonian on the grid.
/// For low-i: single anti-aligned equilibrium near delta_varpi = pi.
/// For high-i: three-lobed structure.
pub fn find_equilibria(
    config: &PhasePortraitConfig,
    params: &SecularHamiltonianParams,
) -> Vec<PhasePoint> {
    let dv_step = std::f64::consts::TAU / config.n_varpi as f64;
    let e_step = (config.e_max - config.e_min) / config.n_e as f64;
    let mut equilibria = Vec::new();

    // Sample at fixed eccentricity e=0.7, scan delta_varpi for extrema
    let e_test = 0.7;
    let mut prev_grad = 0.0_f64;
    for iv in 0..config.n_varpi {
        let dv = iv as f64 * dv_step;
        let h = secular_hamiltonian(config.a, e_test, config.i, 0.0, dv, params);
        let h_next =
            secular_hamiltonian(config.a, e_test, config.i, 0.0, dv + e_step * 0.01, params);
        let grad = h_next - h;

        if iv > 0 && prev_grad * grad < 0.0 {
            equilibria.push(PhasePoint {
                delta_varpi: dv,
                e: e_test,
                h_value: h,
            });
        }
        prev_grad = grad;
    }

    equilibria
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamiltonian::SecularHamiltonianParams;

    #[test]
    fn test_phase_portrait_size() {
        let config = PhasePortraitConfig::low_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_phase_portrait(&config, &params);
        assert_eq!(points.len(), config.n_varpi * config.n_e);
    }

    #[test]
    fn test_phase_portrait_finite_values() {
        let config = PhasePortraitConfig::low_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_phase_portrait(&config, &params);
        for p in &points {
            assert!(
                p.h_value.is_finite(),
                "H should be finite at dv={}, e={}",
                p.delta_varpi,
                p.e
            );
        }
    }

    #[test]
    fn test_high_inclination_portrait() {
        let config = PhasePortraitConfig::high_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_phase_portrait(&config, &params);
        assert!(!points.is_empty());
    }

    #[test]
    fn test_equilibria_found() {
        let config = PhasePortraitConfig::low_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let eq = find_equilibria(&config, &params);
        assert!(!eq.is_empty(), "should find at least one equilibrium");
    }
}
