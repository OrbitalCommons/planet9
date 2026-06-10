//! Miyamoto-Nagai axisymmetric potential model for the Inner Oort Cloud.
//!
//! Psi_MN = -G * M_IOC / sqrt((a_tilde + sqrt(b_tilde^2 + z^2))^2 + r^2)
//!
//! Parameters: a_tilde = 200 AU (radial scale), b_tilde/a_tilde = 5 (scale height ratio).

use p9_core::constants::GM_SUN;
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
    /// Default parameters from the paper: M_IOC = 3 ME, a_tilde = 200 AU, b_tilde/a_tilde = 5.
    pub fn default_paper() -> Self {
        let m_ioc_earth = 3.0;
        let m_ioc_solar = m_ioc_earth * 3.003e-6;
        Self {
            m_ioc_solar,
            a_tilde: 200.0,
            b_tilde: 1000.0, // b_tilde/a_tilde = 5
        }
    }

    /// Low-mass variant: M_IOC = 1 Earth mass.
    pub fn low_mass() -> Self {
        let mut params = Self::default_paper();
        params.m_ioc_solar = 1.0 * 3.003e-6;
        params
    }
}

/// Miyamoto-Nagai gravitational potential at cylindrical coordinates (r, z).
///
/// Psi = -G * M_IOC / sqrt((a + sqrt(b^2 + z^2))^2 + r^2)
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
/// <Psi> = (1/2pi) integral_0^{2pi} Psi(r(f), z(f)) dM
///
/// Approximated via numerical quadrature over the mean anomaly.
pub fn orbit_averaged_potential(
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    params: &MiyamotoNagaiParams,
    n_points: usize,
) -> f64 {
    let dm = std::f64::consts::TAU / n_points as f64;
    let mut sum = 0.0;

    for k in 0..n_points {
        let m_anom = k as f64 * dm;
        let ecc_anom = kepler_solve(m_anom, e);
        let cos_f = (ecc_anom.cos() - e) / (1.0 - e * ecc_anom.cos());
        let sin_f = (1.0 - e * e).sqrt() * ecc_anom.sin() / (1.0 - e * ecc_anom.cos());
        let f = sin_f.atan2(cos_f);

        let r_mag = a * (1.0 - e * e) / (1.0 + e * f.cos());
        // Cylindrical coordinates from orbital elements
        let u = omega + f; // argument of latitude
        let r_cyl = r_mag * (1.0 - (i.sin() * u.sin()).powi(2)).sqrt();
        let z = r_mag * i.sin() * u.sin();

        sum += potential(r_cyl, z, params);
    }

    sum / n_points as f64
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
    fn test_default_params() {
        let params = MiyamotoNagaiParams::default_paper();
        assert!((params.a_tilde - 200.0).abs() < 0.01);
        assert!((params.b_tilde / params.a_tilde - 5.0).abs() < 0.01);
    }
}
