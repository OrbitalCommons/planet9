//! J2/J4 secular quadrupole field for orbit-averaged giant planet effects.
//!
//! When a giant planet is not integrated directly (to save CPU), its
//! gravitational effect can be approximated as an enhanced J2 moment
//! on the central body. This is valid for particles with semi-major axes
//! much larger than the planet's orbit.
//!
//! The effective J2 from an orbit-averaged planet is:
//!   J2_eff = (1/2) * (m_planet / M_sun) * (a_planet / R_ref)^2
//!
//! Reference: Murray & Dermott (1999), Section 6.11

use nalgebra::Vector3;

use crate::constants::*;

/// Orbit-averaged J2 coefficient from a planet.
/// `mass_ratio`: m_planet / M_sun
/// `a_planet_au`: semi-major axis of the planet (AU)
/// `r_ref_au`: reference radius for the J2 expansion (typically 1 AU)
pub fn effective_j2(mass_ratio: f64, a_planet_au: f64, r_ref_au: f64) -> f64 {
    0.5 * mass_ratio * (a_planet_au / r_ref_au).powi(2)
}

/// Combined J2 from orbit-averaging Jupiter, Saturn, and Uranus.
/// Uses a reference radius of 1 AU.
pub fn combined_j2_jsu() -> (f64, f64) {
    let r_ref = 1.0; // AU

    let j2_jup = effective_j2(MASS_JUPITER_SOLAR, 5.2026, r_ref);
    let j2_sat = effective_j2(MASS_SATURN_SOLAR, 9.5549, r_ref);
    let j2_ura = effective_j2(MASS_URANUS_SOLAR, 19.2184, r_ref);

    let j2_total = j2_jup + j2_sat + j2_ura;

    // Effective J4 from the same averaging (second order)
    // J4_eff = -(3/8) * (m/M) * (a/R)^4
    let j4_jup = -0.375 * MASS_JUPITER_SOLAR * (5.2026 / r_ref).powi(4);
    let j4_sat = -0.375 * MASS_SATURN_SOLAR * (9.5549 / r_ref).powi(4);
    let j4_ura = -0.375 * MASS_URANUS_SOLAR * (19.2184 / r_ref).powi(4);

    let j4_total = j4_jup + j4_sat + j4_ura;

    (j2_total, j4_total)
}

/// Apply the J2/J4 secular acceleration to a particle.
/// Uses the central body (Sun) as the origin with the effective J2/J4.
pub fn secular_j2_acceleration(
    pos: &Vector3<f64>,
    gm_sun: f64,
    j2: f64,
    j4: f64,
    r_ref: f64,
) -> Vector3<f64> {
    crate::integrator::kick::j2_acceleration(pos, gm_sun, r_ref, j2, Some(j4))
}
