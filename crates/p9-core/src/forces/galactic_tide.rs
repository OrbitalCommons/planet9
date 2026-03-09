//! Galactic tidal acceleration for outer solar system dynamics.
//!
//! The dominant component is the vertical (z) tide from the Milky Way disk:
//!   a_z = -4π G ρ_MW z
//!
//! where ρ_MW ≈ 0.1 M_sun/pc^3 is the local disk density.
//!
//! This is important for objects with a > 1000 AU (inner Oort Cloud)
//! and becomes the dominant perturbation for a > 10,000 AU.
//!
//! Reference: Heisler & Tremaine (1986), Nesvorny et al. (2017)

use nalgebra::Vector3;

use crate::constants::*;

/// Galactic tidal acceleration on a particle at position `pos` (AU, heliocentric).
/// The z-axis is assumed to be perpendicular to the galactic plane.
///
/// TODO: The ecliptic-to-galactic rotation should be applied for full accuracy,
/// but for the P9 papers the z-tide approximation is used directly.
pub fn galactic_tide_acceleration(pos: &Vector3<f64>) -> Vector3<f64> {
    // 4*pi*G*rho in AU/day^2 per AU
    // G = GM_SUN / M_sun, so G = GM_SUN (in AU^3/(M_sun * day^2))
    // rho = RHO_MW_MSUN_AU3 (in M_sun/AU^3)
    let four_pi_g_rho = 4.0 * std::f64::consts::PI * GM_SUN * RHO_MW_MSUN_AU3;

    // Vertical tide (dominant)
    let az = -four_pi_g_rho * pos.z;

    // Radial tide (subdominant, from differential galactic rotation)
    // a_r ≈ +Ω_0^2 * (x, y, 0) where Ω_0 ~ 26 km/s/kpc is the local angular velocity
    // This is much smaller than the vertical tide for solar system scales
    // TODO: Implement radial tide for completeness

    Vector3::new(0.0, 0.0, az)
}
