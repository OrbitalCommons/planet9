//! The MOND External Field Effect (EFE) quadrupole tidal field.
//!
//! Brown & Mathur (2023) argue that, in Modified Newtonian Dynamics, the Sun's
//! internal gravity in the outer solar system is modified by the *external*
//! galactic field `g_ext` in which the solar system is embedded (the EFE). In
//! the regime where the internal Newtonian field `g_N` is small compared with
//! `g_ext` (which itself is of order the MOND acceleration scale `a0`), the
//! leading anomalous force is a *quadrupole tide* whose symmetry axis points
//! along the external-field direction — i.e. toward the galactic center.
//!
//! Concretely, the anomalous acceleration on a test body at heliocentric
//! position `r` is, to leading order, an axisymmetric tidal field about the
//! Sun→galactic-center unit vector `n̂`, derived as `a = -∇Φ` from the
//! perturbing potential per unit mass
//!
//!   Φ_efe(r) = -(A_efe/2) [ (r·n̂)² - r²/3 ]             (potential form)
//!
//! giving
//!
//!   a_efe(r) = A_efe [ (r·n̂) n̂ - r/3 ].                 (force form)
//!
//! The amplitude `A_efe` (units day⁻²) is the EFE tidal coefficient; its
//! magnitude is fixed by `a0` and the galactic field strength. We carry it as
//! a *labelled* free parameter — the physics of the forced alignment (the
//! headline of the paper) depends only on the *direction* `n̂` and the *sign*
//! of `A_efe`, not its numerical value; the precession *rate* scales linearly
//! with it (pinned in the tests).
//!
//! The galactic-center direction in ecliptic-J2000 coordinates is obtained
//! from `p9_core::coords::sky` (starfield's SPICE `GALACTIC` frame), never
//! hand-rolled.

use nalgebra::Vector3;

use p9_core::coords::sky;
use p9_core::units::{degrees, Angle};

/// MOND acceleration scale `a0`, in m/s² (Milgrom 1983; Brown & Mathur 2023
/// adopt the canonical value). Labelled reference constant; not used to set the
/// tidal amplitude in this reduced model (that is carried as `A_efe`), but
/// pinned so the published number lives in the crate.
pub const MOND_A0_M_S2: f64 = 1.2e-10;

/// Galactic-center alignment claim (Brown & Mathur 2023, abstract): the EFE
/// forces ETNO perihelia toward the galactic-center direction. The galactic
/// center sits at galactic longitude `l = 0`, ecliptic longitude ≈ 267°
/// (the Sgr A* direction). Labelled reference for the headline claim; the
/// test computes the actual ecliptic longitude from the coordinate transforms
/// and compares.
pub const GALACTIC_CENTER_ECLIPTIC_LON_DEG_REF: f64 = 266.8;

/// Unit vector toward the galactic center in heliocentric ecliptic-J2000
/// Cartesian coordinates.
///
/// The galactic center is the galactic-frame +x axis (`l = b = 0`). The
/// ecliptic→galactic rotation's *first row* is the galactic x-axis expressed
/// in ecliptic components, so its transpose is `n̂` in the ecliptic frame.
pub fn galactic_center_ecliptic() -> Vector3<f64> {
    sky::ecliptic_to_galactic_matrix().row(0).transpose()
}

/// Ecliptic longitude (in [0, 360)) of the galactic-center direction as a typed
/// [`Angle`].
pub fn galactic_center_ecliptic_lon() -> Angle {
    let n = galactic_center_ecliptic();
    degrees(
        n.y.atan2(n.x)
            .rem_euclid(std::f64::consts::TAU)
            .to_degrees(),
    )
}

/// Ecliptic latitude of the galactic-center direction as a typed [`Angle`].
pub fn galactic_center_ecliptic_lat() -> Angle {
    let n = galactic_center_ecliptic();
    degrees((n.z / n.norm()).asin().to_degrees())
}

/// EFE quadrupole perturbing potential per unit mass (AU²/day²) at a
/// heliocentric ecliptic position `r` (AU):
///
///   Φ_efe = -(A/2) [ (r·n̂)² - r²/3 ].
///
/// `a_efe` is the EFE tidal amplitude (day⁻²); `n` the (unit) symmetry axis.
pub fn efe_potential(a_efe: f64, n: &Vector3<f64>, r: &Vector3<f64>) -> f64 {
    let rn = r.dot(n);
    -0.5 * a_efe * (rn * rn - r.dot(r) / 3.0)
}

/// EFE anomalous tidal acceleration (AU/day²) at heliocentric position `r`:
///
///   a_efe(r) = -∇Φ = A [ (r·n̂) n̂ - r/3 ].
///
/// This is a stretching (prolate) tide along `n̂` for `A > 0`: it pulls orbits
/// toward the galactic-center axis, the configuration Brown & Mathur identify
/// with the observed apsidal clustering.
pub fn efe_acceleration(a_efe: f64, n: &Vector3<f64>, r: &Vector3<f64>) -> Vector3<f64> {
    a_efe * (r.dot(n) * n - r / 3.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::coords::sky;

    #[test]
    fn typed_galactic_center_angles_match_formula() {
        use uom::si::angle::degree;
        let n = galactic_center_ecliptic();
        let lon_deg =
            n.y.atan2(n.x)
                .rem_euclid(std::f64::consts::TAU)
                .to_degrees();
        let lat_deg = (n.z / n.norm()).asin().to_degrees();
        assert_relative_eq!(
            galactic_center_ecliptic_lon().get::<degree>(),
            lon_deg,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            galactic_center_ecliptic_lat().get::<degree>(),
            lat_deg,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_galactic_center_is_unit_and_round_trips() {
        let n = galactic_center_ecliptic();
        assert_relative_eq!(n.norm(), 1.0, epsilon = 1e-12);
        // Round-trip ecliptic → galactic: n̂ must map to galactic +x (l=b=0).
        let g = sky::ecliptic_to_galactic_matrix() * n;
        assert_relative_eq!(g.x, 1.0, epsilon = 1e-12);
        assert_relative_eq!(g.y, 0.0, epsilon = 1e-12);
        assert_relative_eq!(g.z, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_galactic_center_ecliptic_longitude_matches_reference() {
        use uom::si::angle::degree;
        // Sgr A* lies near ecliptic longitude ~267°; pin the computed value.
        let lon = galactic_center_ecliptic_lon().get::<degree>();
        assert_relative_eq!(lon, GALACTIC_CENTER_ECLIPTIC_LON_DEG_REF, epsilon = 1.0);
        // Galactic center is well south of the ecliptic (β ≈ −5.5°).
        assert!(galactic_center_ecliptic_lat().get::<degree>() < 0.0);
    }

    #[test]
    fn test_acceleration_is_negative_gradient_of_potential() {
        let n = galactic_center_ecliptic();
        let a_efe = 3.7e-13;
        let r = Vector3::new(120.0, -340.0, 80.0);
        let h = 1e-3;
        // Finite-difference -∇Φ vs analytic acceleration.
        let mut grad = Vector3::zeros();
        for k in 0..3 {
            let mut rp = r;
            let mut rm = r;
            rp[k] += h;
            rm[k] -= h;
            grad[k] = (efe_potential(a_efe, &n, &rp) - efe_potential(a_efe, &n, &rm)) / (2.0 * h);
        }
        let a = efe_acceleration(a_efe, &n, &r);
        assert_relative_eq!((a + grad).norm() / a.norm(), 0.0, epsilon = 1e-7);
    }

    #[test]
    fn test_tide_stretches_along_axis_compresses_transverse() {
        // For A>0 the acceleration is outward along n̂ and inward transverse.
        let n = galactic_center_ecliptic();
        let a_efe = 1.0e-12;
        let along = 100.0 * n;
        let a_along = efe_acceleration(a_efe, &n, &along);
        assert!(a_along.dot(&n) > 0.0); // outward along axis
                                        // A transverse displacement (orthogonal to n̂) is pulled back inward.
        let t = n.cross(&Vector3::new(0.0, 0.0, 1.0)).normalize() * 100.0;
        let a_t = efe_acceleration(a_efe, &n, &t);
        assert!(a_t.dot(&t) < 0.0);
    }
}
