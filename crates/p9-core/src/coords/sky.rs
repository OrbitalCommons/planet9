//! Sky-frame transforms (ecliptic ↔ equatorial ↔ galactic, J2000), backed by
//! starfield's `framelib` rotation matrices.
//!
//! Planet9's native representation is heliocentric ecliptic-J2000 Cartesian
//! (nalgebra `Vector3<f64>`); these helpers convert to/from equatorial and
//! galactic frames at the boundary. All matrices come from starfield's SPICE
//! frame table (`ECLIPJ2000`, `GALACTIC`), replacing the obliquity and
//! galactic-pole literals that were previously hand-rolled in the survey
//! crates and the galactic-tide force.
//!
//! Only starfield's offline machinery is used here: the frame matrices and
//! typed coordinates are static data and require no network access.

use std::sync::LazyLock;

use nalgebra::{Matrix3, Vector3};
use starfield::coordinates::cartesian::Cartesian3;
use starfield::framelib::inertial::{Ecliptic, Equatorial, Galactic};
use starfield::framelib::{Frame, ECLIPJ2000_MATRIX, ECLIPJ2000_MATRIX_INV, GALACTIC};
use starfield::time::Timescale;

use crate::constants::{GM_SUN, J2000, RAD2DEG, YEAR_DAYS};
use crate::coords::observer::{apparent_direction_from, EarthState};
use crate::types::OrbitalElements;

/// Equatorial (ICRF/J2000) → galactic rotation, from starfield's SPICE
/// `GALACTIC` frame (static; the epoch argument is ignored by the frame).
static EQ_TO_GAL: LazyLock<Matrix3<f64>> =
    LazyLock::new(|| GALACTIC.rotation_at(&Timescale::default().tt_jd(J2000, None)));

/// Ecliptic J2000 → galactic rotation.
static ECL_TO_GAL: LazyLock<Matrix3<f64>> = LazyLock::new(|| *EQ_TO_GAL * *ECLIPJ2000_MATRIX_INV);

/// Rotation matrix: equatorial J2000 → ecliptic J2000.
pub fn equatorial_to_ecliptic_matrix() -> Matrix3<f64> {
    *ECLIPJ2000_MATRIX
}

/// Rotation matrix: ecliptic J2000 → equatorial J2000.
pub fn ecliptic_to_equatorial_matrix() -> Matrix3<f64> {
    *ECLIPJ2000_MATRIX_INV
}

/// Rotation matrix: equatorial J2000 → galactic.
pub fn equatorial_to_galactic_matrix() -> Matrix3<f64> {
    *EQ_TO_GAL
}

/// Rotation matrix: ecliptic J2000 → galactic.
pub fn ecliptic_to_galactic_matrix() -> Matrix3<f64> {
    *ECL_TO_GAL
}

/// Rotation matrix: galactic → ecliptic J2000.
pub fn galactic_to_ecliptic_matrix() -> Matrix3<f64> {
    ECL_TO_GAL.transpose()
}

/// Unit vector toward the north galactic pole in ecliptic J2000 coordinates
/// (the galactic z-axis row of the ecliptic→galactic rotation).
///
/// λ ≈ 180.023°, β ≈ +29.812°; the complement of β gives the ≈60.19°
/// ecliptic–galactic obliquity.
pub fn galactic_pole_ecliptic() -> Vector3<f64> {
    ECL_TO_GAL.row(2).transpose()
}

/// Equatorial (RA, dec) in radians from ecliptic (λ, β) in radians.
pub fn ecliptic_to_equatorial(lon: f64, lat: f64) -> (f64, f64) {
    let eq: Equatorial = Ecliptic { lon, lat }.into();
    (eq.ra, eq.dec)
}

/// Equatorial (RA, dec) in degrees from ecliptic (λ, β) in degrees.
/// RA is normalized to [0, 360).
pub fn ecliptic_to_equatorial_deg(lon_deg: f64, lat_deg: f64) -> (f64, f64) {
    let (ra, dec) = ecliptic_to_equatorial(lon_deg.to_radians(), lat_deg.to_radians());
    (ra * RAD2DEG, dec * RAD2DEG)
}

/// Ecliptic (λ, β) in radians from equatorial (RA, dec) in radians.
pub fn equatorial_to_ecliptic(ra: f64, dec: f64) -> (f64, f64) {
    let ecl: Ecliptic = Equatorial::new(ra, dec).into();
    (ecl.lon, ecl.lat)
}

/// Ecliptic (λ, β) in degrees from equatorial (RA, dec) in degrees.
/// λ is normalized to [0, 360).
pub fn equatorial_to_ecliptic_deg(ra_deg: f64, dec_deg: f64) -> (f64, f64) {
    let (lon, lat) = equatorial_to_ecliptic(ra_deg.to_radians(), dec_deg.to_radians());
    (
        lon.rem_euclid(std::f64::consts::TAU) * RAD2DEG,
        lat * RAD2DEG,
    )
}

/// Galactic (l, b) in radians from equatorial (RA, dec) in radians.
pub fn equatorial_to_galactic(ra: f64, dec: f64) -> (f64, f64) {
    let gal: Galactic = Equatorial::new(ra, dec).into();
    (gal.lon, gal.lat)
}

/// Equatorial (RA, dec) in degrees of a Cartesian vector in heliocentric
/// (or geocentric) ecliptic-J2000 coordinates. RA is normalized to [0, 360).
pub fn ecliptic_vec_to_equatorial_deg(v: &Vector3<f64>) -> (f64, f64) {
    let eq = *ECLIPJ2000_MATRIX_INV * v;
    let (ra, dec, _) = Cartesian3::from_vector3(eq).to_spherical();
    (ra * RAD2DEG, dec * RAD2DEG)
}

/// Equatorial declination (radians) of a Cartesian vector in ecliptic-J2000
/// coordinates.
pub fn ecliptic_vec_declination(v: &Vector3<f64>) -> f64 {
    let eq = *ECLIPJ2000_MATRIX_INV * v;
    (eq.z / eq.norm()).asin()
}

/// True spherical angular distance (radians) between two direction vectors
/// (any common frame), via starfield's `Cartesian3::angular_distance`.
pub fn angular_distance(a: &Vector3<f64>, b: &Vector3<f64>) -> f64 {
    Cartesian3::from_vector3(*a).angular_distance(&Cartesian3::from_vector3(*b))
}

/// Heliocentric ecliptic position (AU) of an orbit `t_days` after its
/// stored epoch (mean anomaly advanced at the mean motion).
fn heliocentric_position(elem: &OrbitalElements, t_days: f64) -> Vector3<f64> {
    let n = (GM_SUN / (elem.a * elem.a * elem.a)).sqrt(); // rad/day
    let mut advanced = *elem;
    advanced.mean_anomaly = (elem.mean_anomaly + n * t_days).rem_euclid(std::f64::consts::TAU);
    advanced.to_state_vector(GM_SUN).pos
}

/// Apparent geocentric equatorial (RA, dec) in degrees of an orbit
/// `t_days` after its stored epoch, plus the heliocentric distance (AU).
///
/// The object's mean anomaly is advanced by its mean motion; Earth is
/// modeled on a circular coplanar 1-AU orbit (phase zero at t = 0), which is
/// accurate to ~0.03 deg in parallax — ample for footprint membership.
pub fn apparent_position_deg(elem: &OrbitalElements, t_days: f64) -> (f64, f64, f64) {
    let pos = heliocentric_position(elem, t_days);
    let r_helio = pos.norm();

    let theta = std::f64::consts::TAU * t_days / YEAR_DAYS;
    let (sin_t, cos_t) = theta.sin_cos();
    let geo = pos - Vector3::new(cos_t, sin_t, 0.0);

    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&geo);
    (ra, dec, r_helio)
}

/// Ephemeris-grade apparent geocentric equatorial (RA, dec) in degrees of
/// an orbit `t_days` after its stored epoch, observed from the given Earth
/// state, plus the heliocentric distance (AU).
///
/// Unlike the analytic `apparent_position_deg`, this uses a real Earth
/// position and the full apparent chain: light-time retardation (~2-4 days
/// at P9 distances) and stellar aberration.
pub fn apparent_position_with_earth_deg(
    elem: &OrbitalElements,
    earth_state: &EarthState,
    t_days: f64,
) -> (f64, f64, f64) {
    let r_helio = heliocentric_position(elem, t_days).norm();
    let geo = apparent_direction_from(earth_state, |dt| heliocentric_position(elem, t_days + dt));
    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&geo);
    (ra, dec, r_helio)
}

/// Sun–object–Earth phase angle (radians) of an orbit `t_days` after its
/// stored epoch, observed from the given Earth state (instantaneous
/// geometry — the few-day light-time displacement changes α by far less
/// than the H-G law resolves at these distances). Bounded by
/// asin(1 AU / r): below 0.35° everywhere in the BB21 reference
/// population (minimum r ≈ 176 AU), ~0.11° at the fiducial 500 AU.
pub fn phase_angle_with_earth(
    elem: &OrbitalElements,
    earth_state: &EarthState,
    t_days: f64,
) -> f64 {
    let helio = heliocentric_position(elem, t_days);
    let delta = (helio - earth_state.0).norm();
    crate::analysis::photometry::phase_angle(helio.norm(), delta, earth_state.0.norm())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::DEG2RAD;
    use approx::assert_relative_eq;

    /// Mean obliquity of the ecliptic at J2000 implied by starfield's
    /// ECLIPJ2000 matrix (≈ 23.4392911°).
    fn obliquity_deg() -> f64 {
        equatorial_to_ecliptic_matrix()[(1, 1)].acos() * RAD2DEG
    }

    #[test]
    fn test_matrices_are_rotations_and_mutually_inverse() {
        for m in [
            equatorial_to_ecliptic_matrix(),
            equatorial_to_galactic_matrix(),
            ecliptic_to_galactic_matrix(),
        ] {
            assert_relative_eq!(m.determinant(), 1.0, epsilon = 1e-12);
            assert_relative_eq!(
                (m * m.transpose() - Matrix3::identity()).norm(),
                0.0,
                epsilon = 1e-12
            );
        }
        assert_relative_eq!(
            (equatorial_to_ecliptic_matrix() * ecliptic_to_equatorial_matrix()
                - Matrix3::identity())
            .norm(),
            0.0,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            (ecliptic_to_galactic_matrix() * galactic_to_ecliptic_matrix() - Matrix3::identity())
                .norm(),
            0.0,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_equinox_direction_shared_by_ecliptic_and_equatorial() {
        let (ra, dec) = ecliptic_to_equatorial_deg(0.0, 0.0);
        assert!(ra.abs() < 1e-9 || (ra - 360.0).abs() < 1e-9);
        assert!(dec.abs() < 1e-9);
    }

    #[test]
    fn test_north_ecliptic_pole_declination() {
        // The north ecliptic pole sits at dec = 90° − ε, RA = 270°.
        let (ra, dec) = ecliptic_vec_to_equatorial_deg(&Vector3::new(0.0, 0.0, 1.0));
        assert_relative_eq!(dec, 90.0 - obliquity_deg(), epsilon = 1e-9);
        assert_relative_eq!(ra, 270.0, epsilon = 1e-9);
    }

    #[test]
    fn test_ecliptic_equatorial_round_trip_quadrants() {
        for &lon in &[10.0, 100.0, 190.0, 280.0, 359.0] {
            for &lat in &[-80.0, -30.0, 0.0, 45.0, 85.0] {
                let (ra, dec) = ecliptic_to_equatorial_deg(lon, lat);
                let (lon2, lat2) = equatorial_to_ecliptic_deg(ra, dec);
                assert_relative_eq!(lon2, lon, epsilon = 1e-9);
                assert_relative_eq!(lat2, lat, epsilon = 1e-9);
            }
        }
    }

    #[test]
    fn test_galactic_round_trip_through_typed_coordinates() {
        for &(ra_deg, dec_deg) in &[(0.0, 0.0), (192.86, 27.13), (266.4, -28.9), (45.0, -60.0)] {
            let (l, b) = equatorial_to_galactic(ra_deg * DEG2RAD, dec_deg * DEG2RAD);
            let eq: Equatorial = Galactic { lon: l, lat: b }.into();
            // Compare as directions (RA wraps 360 → 0 at the equinox).
            let sep = angular_distance(
                &Cartesian3::from_spherical(eq.ra, eq.dec, 1.0).to_vector3(),
                &Cartesian3::from_spherical(ra_deg * DEG2RAD, dec_deg * DEG2RAD, 1.0).to_vector3(),
            );
            assert!(
                sep < 1e-9 * DEG2RAD,
                "sep = {sep} rad at ({ra_deg}, {dec_deg})"
            );
        }
    }

    #[test]
    fn test_north_galactic_pole_equatorial() {
        // IAU NGP: RA 192.85948°, dec +27.12825° (J2000).
        let ngp = equatorial_to_galactic_matrix().row(2).transpose();
        let (ra, dec, _) = Cartesian3::from_vector3(ngp).to_spherical();
        assert_relative_eq!(ra * RAD2DEG, 192.85948, epsilon = 2e-4);
        assert_relative_eq!(dec * RAD2DEG, 27.12825, epsilon = 2e-4);
    }

    #[test]
    fn test_galactic_pole_ecliptic_orientation() {
        // λ ≈ 180.02°, β ≈ +29.81°; ecliptic-galactic obliquity ≈ 60.19°.
        let n = galactic_pole_ecliptic();
        assert_relative_eq!(n.norm(), 1.0, epsilon = 1e-12);
        let (lon, lat, _) = Cartesian3::from_vector3(n).to_spherical();
        assert_relative_eq!(lon * RAD2DEG, 180.02, epsilon = 0.01);
        assert_relative_eq!(lat * RAD2DEG, 29.81, epsilon = 0.01);
        assert_relative_eq!(n.z.acos() * RAD2DEG, 60.19, epsilon = 0.01);
    }

    #[test]
    fn test_galactic_center_near_galactic_origin() {
        // Sgr A* (RA 266.41683°, dec −29.00781°) lies within ~0.1° of l = b = 0.
        let (l, b) = equatorial_to_galactic(266.41683 * DEG2RAD, -29.00781 * DEG2RAD);
        let sep = angular_distance(
            &Vector3::new(l.cos() * b.cos(), l.sin() * b.cos(), b.sin()),
            &Vector3::new(1.0, 0.0, 0.0),
        );
        assert!(sep * RAD2DEG < 0.1, "sep = {} deg", sep * RAD2DEG);
    }

    #[test]
    fn test_declination_helper_matches_full_conversion() {
        let v = Vector3::new(0.3, -1.2, 0.7);
        let (_, dec_deg) = ecliptic_vec_to_equatorial_deg(&v);
        assert_relative_eq!(
            ecliptic_vec_declination(&v) * RAD2DEG,
            dec_deg,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_angular_distance_orthogonal_and_parallel() {
        let x = Vector3::new(2.0, 0.0, 0.0);
        let z = Vector3::new(0.0, 0.0, 0.5);
        assert_relative_eq!(
            angular_distance(&x, &z),
            std::f64::consts::FRAC_PI_2,
            epsilon = 1e-12
        );
        assert_relative_eq!(angular_distance(&x, &x), 0.0, epsilon = 1e-15);
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
