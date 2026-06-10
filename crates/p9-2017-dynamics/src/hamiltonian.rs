//! Doubly-averaged secular and resonant Hamiltonians from the paper.
//!
//! Three Hamiltonian frameworks:
//! 1. Doubly-averaged secular (Eq. 1): closed-form phase-averaging with giant-planet
//!    precession, rotating frame, and P9-KBO interaction terms.
//! 2. Resonant (Eq. 10): single averaging over lambda_9 under MMR constraint.
//! 3. Resonant-secular (Eq. 12): combined e-Delta_varpi evolution within MMR.
//!
//! Giant planets modeled as effective J2 moment (Eq. 3).

use p9_core::constants::{GM_SUN, TWO_PI};
use serde::{Deserialize, Serialize};

/// Parameters for the doubly-averaged secular Hamiltonian.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SecularHamiltonianParams {
    /// Planet Nine semi-major axis (AU)
    pub a9: f64,
    /// Planet Nine eccentricity
    pub e9: f64,
    /// Planet Nine mass in solar masses
    pub m9_solar: f64,
    /// Effective J2 coefficient from giant planets (AU^2)
    pub j2_eff: f64,
    /// Planet Nine precession rate (rad/yr)
    pub precession_rate_9: f64,
}

impl SecularHamiltonianParams {
    /// Fiducial parameters from the paper: a_9=700, e_9=0.6, m_9=10 ME.
    pub fn default_paper() -> Self {
        let m9_earth = 10.0;
        let m9_solar = m9_earth * 3.003e-6; // Earth mass in solar masses
        Self {
            a9: 700.0,
            e9: 0.6,
            m9_solar,
            j2_eff: compute_j2_effective(),
            precession_rate_9: 0.0,
        }
    }

    /// Alternative parameter set: a_9=600, e_9=0.5, m_9=10 ME.
    pub fn alternative_600() -> Self {
        let m9_earth = 10.0;
        let m9_solar = m9_earth * 3.003e-6;
        Self {
            a9: 600.0,
            e9: 0.5,
            m9_solar,
            j2_eff: compute_j2_effective(),
            precession_rate_9: 0.0,
        }
    }
}

/// Compute the effective J2 coefficient from the four giant planets (Eq. 3).
///
/// J2_eff = (1/2) sum_j (m_j/M_sun) * a_j^2
/// This captures the quadrupolar potential of the giant planet system.
pub fn compute_j2_effective() -> f64 {
    // Giant planet masses in solar masses and semi-major axes in AU
    let giants: [(f64, f64); 4] = [
        (9.548e-4, 5.203),  // Jupiter
        (2.858e-4, 9.537),  // Saturn
        (4.366e-5, 19.189), // Uranus
        (5.151e-5, 30.070), // Neptune
    ];

    giants.iter().map(|(m, a)| 0.5 * m * a * a).sum()
}

/// Perihelion precession rate from the giant planet J2 potential (Eq. 3).
///
/// dot_varpi_J2 = (3/2) * n * J2_eff / (a^2 * (1-e^2)^2)
pub fn precession_rate_j2(a: f64, e: f64, j2_eff: f64) -> f64 {
    let n = (GM_SUN / (a * a * a)).sqrt(); // mean motion (rad/day)
    let eta_sq = (1.0 - e * e) * (1.0 - e * e);
    1.5 * n * j2_eff / (a * a * eta_sq)
}

/// Secular Hamiltonian value at given test particle orbital elements.
///
/// H = H_J2(precession) + H_frame(rotating) + H_interaction(P9-KBO)
///
/// Inputs are test particle elements: a (AU), e, i (rad), omega (rad),
/// delta_varpi = varpi - varpi_9 (rad).
pub fn secular_hamiltonian(
    a: f64,
    e: f64,
    i: f64,
    _omega: f64,
    delta_varpi: f64,
    params: &SecularHamiltonianParams,
) -> f64 {
    let eta = (1.0 - e * e).sqrt();

    // Term 1: Giant-planet J2 precession
    let h_j2 = -1.5 * GM_SUN * params.j2_eff / (a * a * a * eta * eta * eta);

    // Term 2: Rotating frame contribution (linear in action)
    let h_frame = -params.precession_rate_9 * (GM_SUN * a).sqrt() * eta;

    // Term 3: P9-KBO orbit-averaged interaction (leading order)
    let alpha = a / params.a9;
    let alpha_sq = alpha * alpha;
    let eta9 = (1.0 - params.e9 * params.e9).sqrt();
    let prefactor = -GM_SUN * params.m9_solar * alpha_sq / (4.0 * params.a9 * eta9);

    let cos_i = i.cos();
    let cos_dv = delta_varpi.cos();
    let h_int = prefactor * (1.0 + 1.5 * e * e) * (1.0 + 1.5 * (cos_i * cos_i - 1.0))
        + prefactor * 3.75 * e * e * cos_i.powi(2) * cos_dv;

    h_j2 + h_frame + h_int
}

/// Resonant angle for a p:q mean-motion resonance with Planet Nine.
///
/// phi_res = (p * lambda - q * lambda_9 - (p-q) * varpi) / (p-q)
pub fn resonant_angle(lambda: f64, lambda_9: f64, varpi: f64, p: i64, q: i64) -> f64 {
    let pf = p as f64;
    let qf = q as f64;
    let pmq = pf - qf;
    if pmq.abs() < 1e-10 {
        return 0.0;
    }
    let phi = (pf * lambda - qf * lambda_9 - pmq * varpi) / pmq;
    // Normalize to [0, 2pi)
    ((phi % TWO_PI) + TWO_PI) % TWO_PI
}

/// High-inclination secular resonance angle (theta).
///
/// theta = 2*Omega - varpi - varpi_9
/// Conjugate action: Theta = sqrt(1-e^2)/2 * (1 - cos(i))
pub fn high_inclination_angle(omega_big: f64, varpi: f64, varpi_9: f64) -> f64 {
    let theta = 2.0 * omega_big - varpi - varpi_9;
    ((theta % TWO_PI) + TWO_PI) % TWO_PI
}

/// Conjugate action for the high-inclination secular resonance.
pub fn high_inclination_action(e: f64, i: f64) -> f64 {
    (1.0 - e * e).sqrt() / 2.0 * (1.0 - i.cos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_j2_effective_positive() {
        let j2 = compute_j2_effective();
        assert!(j2 > 0.0, "J2_eff = {}", j2);
        // Should be dominated by Jupiter: ~0.5 * 9.5e-4 * 27 ~ 0.013
        assert!(
            j2 > 0.01 && j2 < 0.1,
            "J2_eff = {} out of expected range",
            j2
        );
    }

    #[test]
    fn test_precession_rate_positive() {
        let j2 = compute_j2_effective();
        let rate = precession_rate_j2(300.0, 0.7, j2);
        assert!(rate > 0.0, "precession rate should be positive");
    }

    #[test]
    fn test_fiducial_params() {
        let params = SecularHamiltonianParams::default_paper();
        assert!((params.a9 - 700.0).abs() < 0.01);
        assert!((params.e9 - 0.6).abs() < 0.01);
    }

    #[test]
    fn test_secular_hamiltonian_finite() {
        let params = SecularHamiltonianParams::default_paper();
        let h = secular_hamiltonian(300.0, 0.7, 20.0 * DEG2RAD, 0.0, 0.5, &params);
        assert!(h.is_finite(), "Hamiltonian should be finite, got {}", h);
    }

    #[test]
    fn test_resonant_angle_normalized() {
        let phi = resonant_angle(1.0, 0.5, 0.3, 3, 2);
        assert!(phi >= 0.0 && phi < TWO_PI, "phi = {}", phi);
    }

    #[test]
    fn test_high_inclination_action_zero_at_zero_i() {
        let theta = high_inclination_action(0.5, 0.0);
        assert!(
            theta.abs() < 1e-10,
            "action at i=0 should be ~0, got {}",
            theta
        );
    }

    #[test]
    fn test_high_inclination_action_positive() {
        let theta = high_inclination_action(0.5, 30.0 * DEG2RAD);
        assert!(theta > 0.0, "action should be positive at i=30deg");
    }
}
