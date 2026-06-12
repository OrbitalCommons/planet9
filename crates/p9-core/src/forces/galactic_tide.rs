//! Galactic tidal acceleration for outer solar system dynamics.
//!
//! The dominant component is the vertical tide from the Milky Way disk,
//!   a_z̃ = -4π G ρ_MW z̃,
//! where z̃ is the height above the *galactic* mid-plane and
//! ρ_MW ≈ 0.1 M_sun/pc^3 is the local disk density.
//!
//! The galactic plane is inclined ~60.2° to the ecliptic, so the tide must be
//! evaluated in the galactic frame and rotated back; applying it along
//! ecliptic z misdirects the torque by the full obliquity. The galactic
//! orientation comes from starfield's SPICE `GALACTIC` frame via
//! `coords::sky` (previously a hand-rolled λ=180.02°, β=29.81° pole; the
//! starfield-derived pole is λ=180.0232°, β=29.8114°, an ~11 arcsec shift,
//! well inside all test tolerances).
//!
//! The radial tide terms (Heisler & Tremaine's G1, G2, from the galactic
//! rotation curve) are ~10x weaker than the disk term at the Sun's location
//! and are omitted here; only the G3 (vertical) term is applied.
//!
//! This is important for objects with a > 1000 AU (inner Oort Cloud)
//! and becomes the dominant perturbation for a > 10,000 AU.
//!
//! Reference: Heisler & Tremaine (1986), Nesvorny et al. (2017)

use nalgebra::Vector3;

use crate::constants::*;
use crate::coords::sky;

pub use crate::coords::sky::galactic_pole_ecliptic;

/// Galactic tidal acceleration on a particle at position `pos`
/// (AU, heliocentric ecliptic J2000).
///
/// The position is rotated into the galactic frame, the vertical disk tide
/// (0, 0, -4πGρ z̃) is applied along galactic z, and the result is rotated
/// back to the ecliptic frame.
pub fn galactic_tide_acceleration(pos: &Vector3<f64>) -> Vector3<f64> {
    // 4*pi*G*rho in day^-2:
    // G = GM_SUN / M_sun (AU^3/(M_sun*day^2)), rho in M_sun/AU^3.
    let four_pi_g_rho = 4.0 * std::f64::consts::PI * GM_SUN * RHO_MW_MSUN_AU3;

    let ecl_to_gal = sky::ecliptic_to_galactic_matrix();
    let pos_gal = ecl_to_gal * pos;
    let acc_gal = Vector3::new(0.0, 0.0, -four_pi_g_rho * pos_gal.z);

    ecl_to_gal.transpose() * acc_gal
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_pole_is_unit_and_inclined_60deg() {
        let n = galactic_pole_ecliptic();
        assert_relative_eq!(n.norm(), 1.0, epsilon = 1e-12);
        // Obliquity between galactic plane and ecliptic: acos(n.z) ~ 60.2 deg
        let obliquity = n.z.acos() * RAD2DEG;
        assert!((obliquity - 60.19).abs() < 0.1, "obliquity {obliquity}");
    }

    #[test]
    fn test_tide_restoring_toward_galactic_plane() {
        // A point displaced along the galactic pole feels a pure restoring force.
        let n = galactic_pole_ecliptic();
        let pos = 10_000.0 * n;
        let a = galactic_tide_acceleration(&pos);
        // Anti-parallel to displacement
        assert!(a.dot(&pos) < 0.0);
        assert_relative_eq!(a.cross(&pos).norm(), 0.0, epsilon = 1e-20);
    }

    #[test]
    fn test_tide_vanishes_in_galactic_plane() {
        // Any vector orthogonal to the pole sits in the galactic mid-plane.
        let n = galactic_pole_ecliptic();
        let in_plane = n.cross(&Vector3::new(0.0, 0.0, 1.0)).normalize() * 10_000.0;
        let a = galactic_tide_acceleration(&in_plane);
        assert!(a.norm() < 1e-25);
    }

    #[test]
    fn test_tide_matches_pole_projection_form() {
        // The rotate-apply-rotate-back implementation must equal the
        // algebraically collapsed form a = -4πGρ (r·n̂) n̂.
        let four_pi_g_rho = 4.0 * std::f64::consts::PI * GM_SUN * RHO_MW_MSUN_AU3;
        let n = galactic_pole_ecliptic();
        let pos = Vector3::new(3_000.0, -7_000.0, 5_000.0);
        let expected = -four_pi_g_rho * pos.dot(&n) * n;
        let a = galactic_tide_acceleration(&pos);
        assert_relative_eq!((a - expected).norm(), 0.0, epsilon = 1e-30);
    }
}
