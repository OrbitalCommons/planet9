//! Orbit → sky-position mapping for the survey footprint model.
//!
//! Heliocentric ecliptic longitude/latitude come from the Cartesian state of
//! the orbital elements (p9-core conversion); equatorial declination follows
//! from p9-core's `coords::sky` ecliptic → equatorial transform (starfield's
//! ECLIPJ2000 frame matrix). For the distant-object footprint question
//! (hundreds of AU) the parallax between heliocentric and geocentric
//! directions is < 0.2 degrees and is neglected.

use p9_core::constants::{DEG2RAD, GM_SUN};
use p9_core::coords::sky::ecliptic_vec_declination;
use p9_core::types::OrbitalElements;
use p9_core::units::{degrees, Angle};

/// Heliocentric ecliptic longitude, latitude (radians) and distance (AU)
/// for orbital elements at their stored mean anomaly.
pub fn ecliptic_position(elem: &OrbitalElements) -> (f64, f64, f64) {
    let state = elem.to_state_vector(GM_SUN);
    let r = state.pos.norm();
    let lambda = state.pos.y.atan2(state.pos.x);
    let beta = (state.pos.z / r).asin();
    (lambda, beta, r)
}

/// Declination (degrees) of an orbit's current sky position.
pub fn declination_deg(elem: &OrbitalElements) -> f64 {
    let state = elem.to_state_vector(GM_SUN);
    ecliptic_vec_declination(&state.pos) / DEG2RAD
}

/// Declination of an orbit's current sky position as a dimension-checked
/// [`Angle`]. Typed sibling of [`declination_deg`].
pub fn declination(elem: &OrbitalElements) -> Angle {
    degrees(declination_deg(elem))
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::coords::sky::ecliptic_to_equatorial;

    fn elem(i_deg: f64, omega_deg: f64, omega_big_deg: f64, m_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: 400.0,
            e: 0.2,
            i: i_deg * DEG2RAD,
            omega: omega_deg * DEG2RAD,
            omega_big: omega_big_deg * DEG2RAD,
            mean_anomaly: m_deg * DEG2RAD,
        }
    }

    #[test]
    fn typed_declination_matches_f64() {
        use approx::assert_relative_eq;
        use uom::si::angle::degree;
        let e = elem(15.0, 30.0, 45.0, 60.0);
        assert_relative_eq!(
            declination(&e).get::<degree>(),
            declination_deg(&e),
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_ecliptic_orbit_declination_range() {
        // A coplanar (i = 0) orbit's declination must track the ecliptic:
        // |δ| ≤ ε, reaching ±ε at λ = ±90°.
        for m in [0.0, 45.0, 90.0, 180.0, 270.0] {
            let d = declination_deg(&elem(0.0, 0.0, 0.0, m));
            assert!(d.abs() <= 23.44 + 1e-6, "δ = {d:.2} for M = {m}");
        }
        // At λ = 90°, β = 0 the declination equals the obliquity.
        let (_, dec) = ecliptic_to_equatorial(90.0 * DEG2RAD, 0.0);
        assert!((dec / DEG2RAD - 23.439).abs() < 1e-3);
    }

    #[test]
    fn test_inclined_orbit_reaches_higher_declination() {
        // i = 30° orbit can reach |δ| up to ~ i + ε.
        let mut max_d: f64 = 0.0;
        for m in 0..72 {
            let d = declination_deg(&elem(30.0, 90.0, 0.0, m as f64 * 5.0));
            max_d = max_d.max(d.abs());
        }
        assert!(max_d > 30.0, "max |δ| = {max_d:.1}");
        assert!(max_d < 30.0 + 23.5, "max |δ| = {max_d:.1}");
    }

    #[test]
    fn test_distance_at_perihelion() {
        let e = elem(10.0, 0.0, 0.0, 0.0);
        let (_, _, r) = ecliptic_position(&e);
        assert!((r - 320.0).abs() < 0.5, "r = {r:.1}");
    }
}
