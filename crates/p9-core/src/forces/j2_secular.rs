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
//! Absorbing a planet into the central body also requires the monopole boost
//! ΔGM = GM_planet (ΣGM_JSU ≈ 1.28e-3 GM_sun); without it mean motions and
//! secular frequencies are systematically off by ~6e-4.
//!
//! Validity: the J2+J4 ring expansion is in powers of (a_planet/r)². With
//! a_Uranus = 19.2 AU it converges slowly for perihelia q ≲ 30 AU — at
//! q ≈ 30 AU the J6 term is still at the percent level of J4. Do not use
//! this approximation for Neptune-crossing perihelia; integrate Neptune
//! directly instead.
//!
//! Reference: Murray & Dermott (1999), Section 6.11; Batygin & Brown (2016)

use nalgebra::Vector3;

use crate::constants::*;

/// Orbit-averaged J2 coefficient from a planet.
/// `mass_ratio`: m_planet / M_sun
/// `a_planet_au`: semi-major axis of the planet (AU)
/// `r_ref_au`: reference radius for the J2 expansion (typically 1 AU)
pub fn effective_j2(mass_ratio: f64, a_planet_au: f64, r_ref_au: f64) -> f64 {
    0.5 * mass_ratio * (a_planet_au / r_ref_au).powi(2)
}

/// Combined effective field from orbit-averaging Jupiter, Saturn, and Uranus.
/// Uses a reference radius of 1 AU.
///
/// Returns `(j2_eff, j4_eff, gm_boost)` where `gm_boost = Σ GM_JSU`
/// (AU³/day²) is the monopole that must be added to the central GM.
pub fn combined_j2_jsu() -> (f64, f64, f64) {
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

    // Monopole: the averaged planets' mass joins the central body.
    let gm_boost = GM_JUPITER + GM_SATURN + GM_URANUS;

    (j2_total, j4_total, gm_boost)
}

/// Perturbation acceleration of the J2/J4-averaged planets, relative to bare
/// Keplerian motion about `gm_central`.
///
/// Includes both the ring field (J2/J4 of the boosted central mass) and the
/// extra monopole `-gm_boost * r / r³`. Apply this in the kick stage of an
/// integrator whose Kepler drift uses `gm_central` alone.
pub fn secular_j2_acceleration(
    pos: &Vector3<f64>,
    gm_central: f64,
    j2: f64,
    j4: f64,
    gm_boost: f64,
    r_ref: f64,
) -> Vector3<f64> {
    let r_mag = pos.norm();
    if r_mag < 1e-30 {
        return Vector3::zeros();
    }

    let gm_eff = gm_central + gm_boost;
    let mut accel = crate::integrator::kick::j2_acceleration(pos, gm_eff, r_ref, j2, Some(j4));
    accel -= gm_boost * pos / r_mag.powi(3);
    accel
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gm_boost_magnitude() {
        // ΣGM_JSU ≈ 1.28e-3 GM_sun
        let (_, _, gm_boost) = combined_j2_jsu();
        let ratio = gm_boost / GM_SUN;
        assert!(
            (ratio - 1.28e-3).abs() < 5e-5,
            "GM boost ratio {ratio:.3e} should be ~1.28e-3"
        );
    }

    #[test]
    fn test_secular_acceleration_dominated_by_monopole_at_large_r() {
        // Far from the rings the perturbation tends to the extra monopole.
        let (j2, j4, gm_boost) = combined_j2_jsu();
        let pos = Vector3::new(500.0, 0.0, 0.0);
        let a = secular_j2_acceleration(&pos, GM_SUN, j2, j4, gm_boost, 1.0);
        let monopole = gm_boost / (500.0_f64 * 500.0);
        assert!((a.norm() - monopole).abs() / monopole < 1e-3);
        // Pointing inward (extra attraction).
        assert!(a.x < 0.0);
    }
}
