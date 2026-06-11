//! Phase-space portrait generation for the secular P9–KBO dynamics.
//!
//! Two portraits, following Batygin & Morbidelli (2017):
//!
//! 1. The (Δϖ, e) portrait at fixed semi-major axis and fixed
//!    H_z = √(1−e²)·cos i. H_z is the conserved vertical angular-momentum
//!    action of the Δϖ-averaged problem, so inclination *must* co-vary with
//!    eccentricity along a trajectory: i = acos(H_z/√(1−e²)). (The previous
//!    implementation held i fixed while varying e, which violates the
//!    conservation law the portrait is supposed to respect.)
//!
//! 2. The high-inclination (θ = 2Ω − ϖ − ϖ₉, Θ) portrait, with conjugate
//!    action Θ = √(1−e²)(1 − cos i)/2, evaluated by averaging the
//!    Hamiltonian over the circulating apsidal angle Δϖ at fixed θ
//!    (ω = (Δϖ − θ)/2).

use p9_core::constants::DEG2RAD;
use serde::{Deserialize, Serialize};

use crate::hamiltonian::{secular_hamiltonian, SecularHamiltonianParams};

/// A point on a (Δϖ, e) phase-space portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhasePoint {
    /// Delta varpi = varpi - varpi_9 (rad)
    pub delta_varpi: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination (rad) implied by H_z at this eccentricity
    pub i: f64,
    /// Hamiltonian value at this point
    pub h_value: f64,
}

/// Configuration for phase portrait generation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PhasePortraitConfig {
    /// Semi-major axis (AU) of the test particle
    pub a: f64,
    /// Conserved vertical action H_z = √(1−e²)·cos i labeling the portrait.
    pub h_z: f64,
    /// Number of grid points in delta_varpi
    pub n_varpi: usize,
    /// Number of grid points in eccentricity
    pub n_e: usize,
    /// Minimum eccentricity
    pub e_min: f64,
    /// Maximum eccentricity (clipped to the kinematic bound √(1−H_z²))
    pub e_max: f64,
}

impl PhasePortraitConfig {
    /// H_z for a given (e, i) reference point.
    pub fn h_z_of(e: f64, i: f64) -> f64 {
        (1.0 - e * e).sqrt() * i.cos()
    }

    /// Portrait labeled by the H_z of a low-inclination orbit
    /// (i = 10° at e = 0.7) at a = 300 AU.
    pub fn low_inclination() -> Self {
        Self {
            a: 300.0,
            h_z: Self::h_z_of(0.7, 10.0 * DEG2RAD),
            n_varpi: 200,
            n_e: 100,
            e_min: 0.3,
            e_max: 0.95,
        }
    }

    /// Portrait labeled by the H_z of a high-inclination orbit
    /// (i = 60° at e = 0.7) at a = 300 AU.
    pub fn high_inclination() -> Self {
        Self {
            a: 300.0,
            h_z: Self::h_z_of(0.7, 60.0 * DEG2RAD),
            n_varpi: 200,
            n_e: 100,
            e_min: 0.3,
            e_max: 0.95,
        }
    }

    /// Largest eccentricity compatible with this portrait's H_z
    /// (√(1−e²) ≥ |H_z| requires e ≤ √(1−H_z²)).
    pub fn e_kinematic_max(&self) -> f64 {
        (1.0 - self.h_z * self.h_z).max(0.0).sqrt()
    }

    /// Inclination implied by H_z at eccentricity e.
    pub fn inclination_at(&self, e: f64) -> f64 {
        (self.h_z / (1.0 - e * e).sqrt()).clamp(-1.0, 1.0).acos()
    }
}

/// Generate a grid of Hamiltonian values over the (delta_varpi, e) plane at
/// fixed H_z: each grid point's inclination is i = acos(H_z/√(1−e²)).
/// Eccentricities above the kinematic bound √(1−H_z²) are skipped.
pub fn compute_phase_portrait(
    config: &PhasePortraitConfig,
    params: &SecularHamiltonianParams,
) -> Vec<PhasePoint> {
    let e_cap = config.e_max.min(config.e_kinematic_max() - 1e-9);
    let mut points = Vec::with_capacity(config.n_varpi * config.n_e);
    let dv_step = std::f64::consts::TAU / config.n_varpi as f64;
    let e_step = (e_cap - config.e_min) / config.n_e as f64;

    for iv in 0..config.n_varpi {
        let dv = iv as f64 * dv_step;
        for ie in 0..config.n_e {
            let e = config.e_min + ie as f64 * e_step;
            let i = config.inclination_at(e);
            let h = secular_hamiltonian(config.a, e, i, 0.0, dv, params);
            points.push(PhasePoint {
                delta_varpi: dv,
                e,
                i,
                h_value: h,
            });
        }
    }

    points
}

/// A point on the high-inclination (θ, Θ) portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ThetaPoint {
    /// θ = 2Ω − ϖ − ϖ₉ (rad)
    pub theta: f64,
    /// Conjugate action Θ = √(1−e²)(1 − cos i)/2
    pub capital_theta: f64,
    /// Inclination (rad) implied by Θ at the portrait's eccentricity
    pub i: f64,
    /// Δϖ-averaged Hamiltonian value
    pub h_value: f64,
}

/// Generate the (θ = 2Ω − ϖ − ϖ₉, Θ) portrait at fixed (a, e).
///
/// θ is the resonant combination of the high-inclination secular resonance;
/// its conjugate action Θ = √(1−e²)(1 − cos i)/2 determines i at fixed e.
/// The non-resonant apsidal angle Δϖ circulates, so the portrait Hamiltonian
/// is the Δϖ-average of H at fixed θ, with ω = (Δϖ − θ)/2 (from
/// θ = 2Ω − ϖ − ϖ₉ and Δϖ = ϖ − ϖ₉).
pub fn compute_theta_portrait(
    a: f64,
    e: f64,
    n_theta: usize,
    n_action: usize,
    params: &SecularHamiltonianParams,
) -> Vec<ThetaPoint> {
    let eta = (1.0 - e * e).sqrt();
    let theta_step = std::f64::consts::TAU / n_theta as f64;
    // Θ ranges over (0, η): i from 0 to 180°. Stay off the exact poles.
    let action_step = eta / (n_action as f64 + 1.0);
    let n_dv = 32;
    let dv_step = std::f64::consts::TAU / n_dv as f64;

    let mut points = Vec::with_capacity(n_theta * n_action);
    for it in 0..n_theta {
        let theta = it as f64 * theta_step;
        for ia in 1..=n_action {
            let cap_theta = ia as f64 * action_step;
            let i = (1.0 - 2.0 * cap_theta / eta).clamp(-1.0, 1.0).acos();

            let mut h_sum = 0.0;
            for idv in 0..n_dv {
                let dv = (idv as f64 + 0.5) * dv_step;
                let omega = 0.5 * (dv - theta);
                h_sum += secular_hamiltonian(a, e, i, omega, dv, params);
            }

            points.push(ThetaPoint {
                theta,
                capital_theta: cap_theta,
                i,
                h_value: h_sum / n_dv as f64,
            });
        }
    }

    points
}

/// Find the equilibrium points (fixed points) of the (Δϖ, e) portrait along
/// the e = e_test slice: sign changes of ∂H/∂Δϖ located by central
/// differences with a step proportional to the Δϖ grid spacing (the previous
/// version used the *eccentricity* step as the Δϖ finite-difference step).
///
/// For low inclination there are equilibria at the symmetry points Δϖ = 0
/// and Δϖ = π; the anti-aligned (Δϖ = π) one is the libration center of the
/// confined population.
pub fn find_equilibria(
    config: &PhasePortraitConfig,
    params: &SecularHamiltonianParams,
) -> Vec<PhasePoint> {
    let dv_step = std::f64::consts::TAU / config.n_varpi as f64;
    let fd_step = 0.5 * dv_step;
    let mut equilibria = Vec::new();

    let e_test = 0.7_f64.min(config.e_kinematic_max() - 1e-3);
    let i_test = config.inclination_at(e_test);

    let grad_at = |dv: f64| {
        let h_plus = secular_hamiltonian(config.a, e_test, i_test, 0.0, dv + fd_step, params);
        let h_minus = secular_hamiltonian(config.a, e_test, i_test, 0.0, dv - fd_step, params);
        (h_plus - h_minus) / (2.0 * fd_step)
    };

    let mut prev_grad = grad_at(-0.5 * dv_step);
    for iv in 0..config.n_varpi {
        let dv = (iv as f64 + 0.5) * dv_step;
        let grad = grad_at(dv);
        if prev_grad * grad < 0.0 {
            // Refine by linear interpolation of the gradient zero.
            let frac = prev_grad / (prev_grad - grad);
            let dv_eq = dv - dv_step + frac * dv_step;
            let h = secular_hamiltonian(config.a, e_test, i_test, 0.0, dv_eq, params);
            equilibria.push(PhasePoint {
                delta_varpi: dv_eq.rem_euclid(std::f64::consts::TAU),
                e: e_test,
                i: i_test,
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
    use crate::hamiltonian::{high_inclination_action, SecularHamiltonianParams};
    use std::f64::consts::PI;

    #[test]
    fn test_phase_portrait_size() {
        let config = PhasePortraitConfig::low_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_phase_portrait(&config, &params);
        assert_eq!(points.len(), config.n_varpi * config.n_e);
    }

    #[test]
    fn test_phase_portrait_conserves_h_z() {
        // The defining fix (M6): every grid point of a portrait must share
        // the same H_z = √(1−e²)cos i.
        let config = PhasePortraitConfig::high_inclination();
        let params = SecularHamiltonianParams::default_paper();
        for p in compute_phase_portrait(&config, &params) {
            let h_z = (1.0 - p.e * p.e).sqrt() * p.i.cos();
            assert!(
                (h_z - config.h_z).abs() < 1e-9,
                "H_z = {h_z:.6} != {:.6} at e = {:.3}",
                config.h_z,
                p.e
            );
            assert!(p.h_value.is_finite());
        }
    }

    #[test]
    fn test_inclination_anticorrelated_with_eccentricity() {
        // At fixed H_z > 0, increasing e decreases √(1−e²), so cos i =
        // H_z/√(1−e²) increases and i *decreases* — the Kozai-Lidov e–i
        // trade-off the portrait must respect.
        let config = PhasePortraitConfig::low_inclination();
        let i_low_e = config.inclination_at(0.3);
        let i_high_e = config.inclination_at(0.69);
        assert!(i_high_e < i_low_e);
    }

    #[test]
    fn test_theta_portrait_shape_and_finite() {
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_theta_portrait(300.0, 0.5, 24, 10, &params);
        assert_eq!(points.len(), 24 * 10);
        for p in &points {
            assert!(p.h_value.is_finite());
            // Θ and the implied inclination are consistent.
            let cap = high_inclination_action(0.5, p.i);
            assert!((cap - p.capital_theta).abs() < 1e-9);
        }
    }

    #[test]
    fn test_theta_portrait_depends_on_theta() {
        // The high-inclination resonance term (∝ cos 2ω through the
        // quadrupole) must survive the Δϖ average: H̄ varies with θ.
        let params = SecularHamiltonianParams::default_paper();
        let points = compute_theta_portrait(300.0, 0.5, 16, 6, &params);
        let row: Vec<&ThetaPoint> = points
            .iter()
            .filter(|p| (p.capital_theta - points[3].capital_theta).abs() < 1e-12)
            .collect();
        let h_min = row.iter().map(|p| p.h_value).fold(f64::INFINITY, f64::min);
        let h_max = row
            .iter()
            .map(|p| p.h_value)
            .fold(f64::NEG_INFINITY, f64::max);
        assert!(
            h_max - h_min > 0.0,
            "Δϖ-averaged H should depend on θ (range {:.3e})",
            h_max - h_min
        );
    }

    #[test]
    fn test_equilibria_at_symmetry_points() {
        // Low-inclination portrait: equilibria at the symmetry points of the
        // cos Δϖ octupole term — including the anti-aligned libration center
        // at Δϖ = π.
        let config = PhasePortraitConfig::low_inclination();
        let params = SecularHamiltonianParams::default_paper();
        let eq = find_equilibria(&config, &params);
        assert!(!eq.is_empty(), "should find equilibria");

        let near_pi = eq
            .iter()
            .any(|p| (p.delta_varpi - PI).abs() < 2.0 * std::f64::consts::TAU / 200.0);
        assert!(
            near_pi,
            "anti-aligned equilibrium near Δϖ = π expected, got {:?}",
            eq.iter().map(|p| p.delta_varpi).collect::<Vec<_>>()
        );
    }
}
