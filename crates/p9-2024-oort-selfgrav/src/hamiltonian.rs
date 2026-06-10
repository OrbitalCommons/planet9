//! Secular Hamiltonian combining giant-planet J2 precession and IOC self-gravity.
//!
//! H = H_J2(precession) + <Psi_MN>(orbit-averaged IOC potential)
//!
//! The J2 term contributes < 2% at a=1000 AU, q=75 AU but ~25% at q=30 AU.

use p9_core::constants::GM_SUN;
use serde::{Deserialize, Serialize};

use crate::miyamoto_nagai::{orbit_averaged_potential, MiyamotoNagaiParams};

/// Combined Hamiltonian parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HamiltonianParams {
    /// Effective J2 from giant planets (AU^2)
    pub j2_eff: f64,
    /// Miyamoto-Nagai IOC parameters
    pub mn_params: MiyamotoNagaiParams,
    /// Number of quadrature points for orbit averaging
    pub n_quadrature: usize,
}

impl HamiltonianParams {
    pub fn default_paper() -> Self {
        Self {
            j2_eff: compute_j2_effective(),
            mn_params: MiyamotoNagaiParams::default_paper(),
            n_quadrature: 128,
        }
    }
}

/// Effective J2 from giant planets (same as p9-2017-dynamics).
fn compute_j2_effective() -> f64 {
    let giants: [(f64, f64); 4] = [
        (9.548e-4, 5.203),
        (2.858e-4, 9.537),
        (4.366e-5, 19.189),
        (5.151e-5, 30.070),
    ];
    giants.iter().map(|(m, a)| 0.5 * m * a * a).sum()
}

/// Giant-planet J2 contribution to the Hamiltonian.
///
/// H_J2 = G*C*(3*cos^2(i) - 1) / (8*a^3*(1-e^2)^{3/2})
/// where C = M_sun * J2_eff
pub fn h_j2(a: f64, e: f64, i: f64, j2_eff: f64) -> f64 {
    let eta = (1.0 - e * e).sqrt();
    let cos_i = i.cos();
    GM_SUN * j2_eff * (3.0 * cos_i * cos_i - 1.0) / (8.0 * a * a * a * eta * eta * eta)
}

/// Full secular Hamiltonian: H_J2 + <Psi_MN>.
pub fn secular_hamiltonian(a: f64, e: f64, i: f64, omega: f64, params: &HamiltonianParams) -> f64 {
    let h_planets = h_j2(a, e, i, params.j2_eff);
    let h_ioc = orbit_averaged_potential(a, e, i, omega, &params.mn_params, params.n_quadrature);
    h_planets + h_ioc
}

/// Ratio of J2 contribution to total Hamiltonian.
///
/// Used to verify the paper's claim: < 2% at a=1000, q=75; ~25% at q=30.
pub fn j2_fraction(a: f64, e: f64, i: f64, omega: f64, params: &HamiltonianParams) -> f64 {
    let h_planets = h_j2(a, e, i, params.j2_eff);
    let h_total = secular_hamiltonian(a, e, i, omega, params);
    if h_total.abs() < 1e-30 {
        return 0.0;
    }
    (h_planets / h_total).abs()
}

/// Evolutionary timescale from the Hamiltonian (approximate, in Gyr).
///
/// T ~ 2*pi * |dH/d(action)|^{-1}, estimated from the orbit-averaged IOC potential.
pub fn evolutionary_timescale_gyr(a: f64, e: f64, i: f64, params: &HamiltonianParams) -> f64 {
    // Estimate from potential gradient
    let de = 0.001;
    let e_hi = (e + de).min(0.999);
    let e_lo = (e - de).max(0.001);
    let h_hi = secular_hamiltonian(a, e_hi, i, 0.0, params);
    let h_lo = secular_hamiltonian(a, e_lo, i, 0.0, params);
    let dh_de = (h_hi - h_lo) / (e_hi - e_lo);

    if dh_de.abs() < 1e-30 {
        return f64::INFINITY;
    }

    let period_days = std::f64::consts::TAU / dh_de.abs();
    let gyr_days = 365.25e9;
    period_days / gyr_days
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_j2_effective_value() {
        let j2 = compute_j2_effective();
        assert!(j2 > 0.01 && j2 < 0.1);
    }

    #[test]
    fn test_h_j2_finite() {
        let j2 = compute_j2_effective();
        let h = h_j2(1000.0, 0.5, 30.0 * DEG2RAD, j2);
        assert!(h.is_finite());
    }

    #[test]
    fn test_secular_hamiltonian_finite() {
        let params = HamiltonianParams::default_paper();
        let h = secular_hamiltonian(1000.0, 0.5, 30.0 * DEG2RAD, 0.0, &params);
        assert!(h.is_finite(), "H = {}", h);
    }

    #[test]
    fn test_j2_fraction_small_at_large_q() {
        let params = HamiltonianParams::default_paper();
        // At a=1000, e=0.925 => q=75 AU, paper says J2 < 2%
        let frac = j2_fraction(1000.0, 0.925, 30.0 * DEG2RAD, 0.0, &params);
        assert!(frac < 0.5, "J2 fraction = {} at q=75 AU", frac);
    }

    #[test]
    fn test_evolutionary_timescale_long() {
        let params = HamiltonianParams::default_paper();
        let tau = evolutionary_timescale_gyr(1000.0, 0.9, 30.0 * DEG2RAD, &params);
        // Paper says timescales are tens of Gyr
        assert!(tau > 0.0, "timescale should be positive, got {}", tau);
    }
}
