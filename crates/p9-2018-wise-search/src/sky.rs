//! Orbit → sky-position mapping for the WISE footprint model.
//!
//! WISE's footprint cut is in galactic latitude (the masked galactic-plane
//! confusion band), so we map orbital elements to a heliocentric ecliptic
//! Cartesian state (p9-core) and then to galactic coordinates via
//! `p9_core::coords::sky` (starfield ECLIPJ2000 → GALACTIC frame matrices).
//! At hundreds of AU the heliocentric/geocentric parallax (< 0.2°) is
//! negligible for the footprint question and is neglected.

use p9_core::constants::{GM_SUN, RAD2DEG};
use p9_core::coords::sky::{ecliptic_to_galactic_matrix, ecliptic_vec_to_equatorial_deg};
use p9_core::types::OrbitalElements;

/// Equatorial declination (degrees) of an orbit's current sky position.
pub fn declination_deg(elem: &OrbitalElements) -> f64 {
    let state = elem.to_state_vector(GM_SUN);
    let (_, dec) = ecliptic_vec_to_equatorial_deg(&state.pos);
    dec
}

/// Galactic latitude b (degrees) of an orbit's current sky position.
pub fn galactic_latitude_deg(elem: &OrbitalElements) -> f64 {
    let state = elem.to_state_vector(GM_SUN);
    let gal = ecliptic_to_galactic_matrix() * state.pos;
    (gal.z / gal.norm()).asin() * RAD2DEG
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    fn elem(i_deg: f64, omega_deg: f64, omega_big_deg: f64, m_deg: f64) -> OrbitalElements {
        OrbitalElements {
            a: 500.0,
            e: 0.2,
            i: i_deg * DEG2RAD,
            omega: omega_deg * DEG2RAD,
            omega_big: omega_big_deg * DEG2RAD,
            mean_anomaly: m_deg * DEG2RAD,
        }
    }

    #[test]
    fn galactic_latitude_in_range() {
        for m in 0..360 {
            let b = galactic_latitude_deg(&elem(20.0, 60.0, 120.0, m as f64));
            assert!(b.abs() <= 90.0 + 1e-6, "|b| = {b}");
        }
    }

    #[test]
    fn a_low_inclination_orbit_crosses_the_galactic_plane() {
        // Over a full orbit a near-ecliptic object sweeps a range of galactic
        // latitudes (the ecliptic is inclined ~60° to the galactic plane), so
        // some mean anomalies fall in the masked band and some do not.
        let mut min_b: f64 = 90.0;
        let mut max_b: f64 = -90.0;
        for m in 0..360 {
            let b = galactic_latitude_deg(&elem(5.0, 0.0, 0.0, m as f64));
            min_b = min_b.min(b);
            max_b = max_b.max(b);
        }
        assert!(
            min_b < 10.0 && max_b > 10.0,
            "b range [{min_b:.1}, {max_b:.1}]"
        );
    }
}
