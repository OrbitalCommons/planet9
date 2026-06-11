//! Miyamoto-Nagai axisymmetric potential model for the Inner Oort Cloud.
//!
//! Psi_MN = -G * M_IOC / sqrt(r^2 + (a_tilde + sqrt(b_tilde^2 + z^2))^2)
//!
//! Parameters from Batygin & Nesvorny (2024), Celestial Mechanics and
//! Dynamical Astronomy 136, 24 (arXiv:2405.15139), Section 2: the IOC mass
//! distribution of the Nesvorny et al. (2023) `cluster2` simulation snapshot
//! at t = 300 Myr is fit by a Miyamoto-Nagai profile with
//!
//!   M_IOC = 3 M_earth,  a_tilde = 200 AU,  b_tilde / a_tilde = 5.
//!
//! Note that b_tilde >> a_tilde makes the profile quasi-spherical (the MN
//! model limits to a Plummer sphere of scale a_tilde + b_tilde as
//! b_tilde/a_tilde -> inf), with most of the mass concentrated interior to
//! the ~10^3 AU orbits it perturbs. The residual flattening is what breaks
//! spherical symmetry and drives the vZLK-like secular dynamics.

use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};
use serde::{Deserialize, Serialize};

/// Miyamoto-Nagai potential parameters for the IOC disk.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MiyamotoNagaiParams {
    /// IOC mass in solar masses
    pub m_ioc_solar: f64,
    /// Radial scale length (AU)
    pub a_tilde: f64,
    /// Vertical scale height (AU)
    pub b_tilde: f64,
}

impl MiyamotoNagaiParams {
    /// Default parameters from Batygin & Nesvorny (2024): M_IOC = 3 M_earth,
    /// a_tilde = 200 AU, b_tilde/a_tilde = 5 (see module docs for provenance).
    pub fn default_paper() -> Self {
        Self {
            m_ioc_solar: 3.0 * EARTH_MASS_SOLAR,
            a_tilde: 200.0,
            b_tilde: 1000.0, // b_tilde/a_tilde = 5
        }
    }

    /// Low-mass variant: M_IOC = 1 Earth mass.
    pub fn low_mass() -> Self {
        Self {
            m_ioc_solar: EARTH_MASS_SOLAR,
            ..Self::default_paper()
        }
    }
}

/// Miyamoto-Nagai gravitational potential at cylindrical coordinates (r, z).
///
/// Psi = -G * M_IOC / sqrt(r^2 + (a + sqrt(b^2 + z^2))^2)
///
/// Returns potential in AU^2/day^2 (consistent with GM_SUN units).
pub fn potential(r: f64, z: f64, params: &MiyamotoNagaiParams) -> f64 {
    let gm_ioc = GM_SUN * params.m_ioc_solar;
    let zterm = (params.b_tilde * params.b_tilde + z * z).sqrt();
    let denom = ((params.a_tilde + zterm).powi(2) + r * r).sqrt();
    -gm_ioc / denom
}

/// Radial acceleration (inward, toward r=0) from the MN potential.
///
/// a_r = -dPsi/dr = -G*M*r / ((a + sqrt(b^2+z^2))^2 + r^2)^{3/2}
pub fn acceleration_r(r: f64, z: f64, params: &MiyamotoNagaiParams) -> f64 {
    let gm_ioc = GM_SUN * params.m_ioc_solar;
    let zterm = (params.b_tilde * params.b_tilde + z * z).sqrt();
    let denom = ((params.a_tilde + zterm).powi(2) + r * r).powf(1.5);
    -gm_ioc * r / denom
}

/// Vertical acceleration (toward z=0) from the MN potential.
///
/// a_z = -dPsi/dz = -G*M*z*(a + sqrt(b^2+z^2)) / (sqrt(b^2+z^2) * ((a+sqrt(b^2+z^2))^2+r^2)^{3/2})
pub fn acceleration_z(r: f64, z: f64, params: &MiyamotoNagaiParams) -> f64 {
    let gm_ioc = GM_SUN * params.m_ioc_solar;
    let zterm = (params.b_tilde * params.b_tilde + z * z).sqrt();
    let factor = params.a_tilde + zterm;
    let denom = zterm * (factor * factor + r * r).powf(1.5);
    if denom.abs() < 1e-30 {
        return 0.0;
    }
    -gm_ioc * z * factor / denom
}

/// Orbit-averaged potential for a Keplerian orbit in the MN disk.
///
/// <Psi> = (1/2pi) ∮ Psi dM = (1/2pi) ∫_0^{2pi} Psi(nu) * (r^2 / (a^2 eta)) dnu
///
/// where eta = sqrt(1-e^2) and the Jacobian dM/dnu = r^2/(a^2 eta) follows
/// from angular-momentum conservation (r^2 dnu/dt = n a^2 eta). Quadrature in
/// true anomaly avoids solving Kepler's equation and, with the midpoint rule
/// on a periodic integrand, converges spectrally; the Jacobian automatically
/// de-weights the (fast-traversed) perihelion region at high eccentricity.
/// See `test_quadrature_converges_high_e` for the n vs 2n convergence check
/// at e = 0.95.
pub fn orbit_averaged_potential(
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    params: &MiyamotoNagaiParams,
    n_points: usize,
) -> f64 {
    let eta_sq = 1.0 - e * e;
    let eta = eta_sq.sqrt();
    let dnu = std::f64::consts::TAU / n_points as f64;
    let mut sum = 0.0;

    for k in 0..n_points {
        // Midpoint rule on the periodic integrand.
        let nu = (k as f64 + 0.5) * dnu;
        let r_mag = a * eta_sq / (1.0 + e * nu.cos());

        // Cylindrical coordinates from orbital elements (axisymmetric
        // potential: only the argument of latitude u matters, not the node).
        let u = omega + nu;
        let sin_lat = i.sin() * u.sin();
        let z = r_mag * sin_lat;
        let r_cyl = r_mag * (1.0 - sin_lat * sin_lat).sqrt();

        // dM = (r^2 / (a^2 eta)) dnu
        let jacobian = r_mag * r_mag / (a * a * eta);
        sum += potential(r_cyl, z, params) * jacobian;
    }

    sum * dnu / std::f64::consts::TAU
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_potential_negative() {
        let params = MiyamotoNagaiParams::default_paper();
        let psi = potential(500.0, 0.0, &params);
        assert!(psi < 0.0, "gravitational potential should be negative");
    }

    #[test]
    fn test_potential_decreases_with_r() {
        let params = MiyamotoNagaiParams::default_paper();
        let psi_near = potential(100.0, 0.0, &params);
        let psi_far = potential(1000.0, 0.0, &params);
        assert!(
            psi_near < psi_far,
            "potential should be more negative closer"
        );
    }

    #[test]
    fn test_acceleration_r_inward() {
        let params = MiyamotoNagaiParams::default_paper();
        let ar = acceleration_r(500.0, 0.0, &params);
        assert!(ar < 0.0, "radial acceleration should point inward");
    }

    #[test]
    fn test_acceleration_z_restoring() {
        let params = MiyamotoNagaiParams::default_paper();
        let az = acceleration_z(500.0, 100.0, &params);
        assert!(az < 0.0, "z acceleration should point toward midplane");
    }

    #[test]
    fn test_orbit_averaged_finite() {
        let params = MiyamotoNagaiParams::default_paper();
        let avg = orbit_averaged_potential(1000.0, 0.5, 30.0 * DEG2RAD, 0.0, &params, 64);
        assert!(avg.is_finite() && avg < 0.0, "avg potential = {}", avg);
    }

    #[test]
    fn test_quadrature_converges_high_e() {
        // n vs 2n convergence at e = 0.95: the true-anomaly Jacobian keeps
        // the quadrature well-conditioned even for near-radial orbits.
        let params = MiyamotoNagaiParams::default_paper();
        let (a, e, i, omega) = (1000.0, 0.95, 30.0 * DEG2RAD, 70.0 * DEG2RAD);
        let v128 = orbit_averaged_potential(a, e, i, omega, &params, 128);
        let v256 = orbit_averaged_potential(a, e, i, omega, &params, 256);
        let rel = ((v128 - v256) / v256).abs();
        assert!(rel < 1e-8, "n vs 2n relative difference = {rel:.2e}");
    }

    #[test]
    fn test_circular_orbit_average_matches_pointwise() {
        // For e = 0 in the midplane (i = 0), <Psi> is exactly Psi(a, 0).
        let params = MiyamotoNagaiParams::default_paper();
        let avg = orbit_averaged_potential(800.0, 0.0, 0.0, 0.0, &params, 64);
        let point = potential(800.0, 0.0, &params);
        assert!(((avg - point) / point).abs() < 1e-12);
    }

    #[test]
    fn test_default_params() {
        let params = MiyamotoNagaiParams::default_paper();
        assert!((params.a_tilde - 200.0).abs() < 0.01);
        assert!((params.b_tilde / params.a_tilde - 5.0).abs() < 0.01);
        assert!((params.m_ioc_solar - 3.0 * EARTH_MASS_SOLAR).abs() < 1e-12);
    }
}
