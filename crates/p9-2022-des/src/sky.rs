//! Orbit -> apparent sky-position mapping for the DES footprint model.
//!
//! The footprint of `survey_model` is defined in *equatorial* coordinates
//! (RA/dec), so orbital sky positions must be rotated from the heliocentric
//! ecliptic frame to the equatorial frame before any footprint test (the
//! analogous fix to p9-2021-ztf's `sky.rs`; substituting ecliptic latitude
//! for declination misclassifies exactly the survey-boundary regions that
//! drive unique-coverage estimates). The rotation itself comes from
//! p9-core's `coords::sky` (starfield's ECLIPJ2000 frame matrix).
//!
//! Apparent (geocentric) positions are computed per observing night: the
//! geocentric vector is the heliocentric vector minus Earth's (circular,
//! coplanar) orbital position, which captures both the annual parallactic
//! oscillation (~1/r radians, i.e. +/-0.14 deg at 400 AU) and the slow
//! orbital drift of the object across the 6-year survey baseline.

use nalgebra::Vector3;
use p9_core::constants::{GM_SUN, YEAR_DAYS};
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::types::OrbitalElements;

/// Apparent geocentric equatorial (RA, dec) in degrees of an orbit
/// `t_days` after its stored epoch, plus the heliocentric distance (AU).
///
/// The object's mean anomaly is advanced by its mean motion; Earth is
/// modeled on a circular coplanar 1-AU orbit (phase zero at t = 0), which is
/// accurate to ~0.03 deg in parallax — ample for footprint membership.
pub fn apparent_position_deg(elem: &OrbitalElements, t_days: f64) -> (f64, f64, f64) {
    let n = (GM_SUN / (elem.a * elem.a * elem.a)).sqrt(); // rad/day
    let mut advanced = *elem;
    advanced.mean_anomaly = (elem.mean_anomaly + n * t_days).rem_euclid(std::f64::consts::TAU);
    let state = advanced.to_state_vector(GM_SUN);
    let r_helio = state.pos.norm();

    let theta = std::f64::consts::TAU * t_days / YEAR_DAYS;
    let (sin_t, cos_t) = theta.sin_cos();
    let geo = state.pos - Vector3::new(cos_t, sin_t, 0.0);

    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&geo);
    (ra, dec, r_helio)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_ecliptic_pole_maps_to_obliquity_complement() {
        // The north ecliptic pole sits at dec = 90 - 23.4393 = 66.56 deg.
        let (_, dec) = ecliptic_vec_to_equatorial_deg(&Vector3::new(0.0, 0.0, 1.0));
        assert!((dec - (90.0 - 23.439_291)).abs() < 1e-6, "dec = {dec}");
    }

    #[test]
    fn test_equinox_direction_unchanged() {
        // The vernal-equinox direction is shared by both frames.
        let (ra, dec) = ecliptic_vec_to_equatorial_deg(&Vector3::new(1.0, 0.0, 0.0));
        assert!(ra.abs() < 1e-9 && dec.abs() < 1e-9);
    }

    #[test]
    fn test_ecliptic_plane_declination_bounded_by_obliquity() {
        // Positions in the ecliptic plane reach at most |dec| = obliquity.
        for k in 0..36 {
            let lambda = k as f64 * 10.0 * DEG2RAD;
            let (_, dec) =
                ecliptic_vec_to_equatorial_deg(&Vector3::new(lambda.cos(), lambda.sin(), 0.0));
            assert!(dec.abs() <= 23.44 + 1e-9, "dec = {dec}");
        }
    }

    #[test]
    fn test_parallax_amplitude_scales_inversely_with_distance() {
        // Over one year the apparent RA of a distant object oscillates with
        // amplitude ~ (1 AU / r) radians.
        let elem = OrbitalElements {
            a: 400.0,
            e: 0.0,
            i: 0.0,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 90.0 * DEG2RAD,
        };
        let mut ra_min = f64::INFINITY;
        let mut ra_max = f64::NEG_INFINITY;
        for k in 0..=24 {
            let (ra, _, _) = apparent_position_deg(&elem, k as f64 * YEAR_DAYS / 24.0);
            ra_min = ra_min.min(ra);
            ra_max = ra_max.max(ra);
        }
        let half_amplitude = (ra_max - ra_min) / 2.0;
        let expected = (1.0_f64 / 400.0).asin() / DEG2RAD; // ~0.143 deg
        assert!(
            (half_amplitude - expected).abs() < 0.05,
            "parallax half-amplitude = {half_amplitude:.3} deg, expected ~{expected:.3}"
        );
    }
}
