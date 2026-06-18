//! Orbital elements of 2015 BP519, the high-inclination ETNO of Becker et
//! al. (2018), arXiv:1805.05355.
//!
//! Two labelled element sets are pinned:
//!
//! * [`discovery_2018`] — the paper's discovery-epoch solution (abstract /
//!   Table values): a ≈ 449 AU, e ≈ 0.92, i ≈ 54°. These are the numbers the
//!   paper's dynamical analysis used and the ones the headline tests pin.
//! * [`jpl_current`] — the current JPL Small-Body Database osculating
//!   solution (looked up 2026-06-14; designation 768325): a ≈ 486 AU,
//!   e ≈ 0.927, i ≈ 54.1°, Ω ≈ 135°, ω ≈ 348°, H = 4.32. As with the rest of
//!   the workspace (`p9_core::data::etno`), JPL keeps refitting as the arc
//!   lengthens, so the live a/e have drifted up the fit degeneracy from the
//!   discovery solution while q, i, Ω, H are stable.
//!
//! The argument of perihelion and node are not in the 2018 abstract; the
//! discovery set therefore borrows ω/Ω from the JPL solution (i, q and the
//! angles are the stable elements) and is flagged as such. What matters for
//! the paper's headline is (a, e, i): a detached (q ≈ 35 AU), extremely
//! eccentric, highly inclined large-a orbit.

use p9_core::constants::{DEG2RAD, TWO_PI};
use p9_core::types::OrbitalElements;
use p9_core::units::{radians, Angle, Length};

/// A labelled orbit solution for 2015 BP519.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bp519 {
    /// Provenance label.
    pub source: &'static str,
    /// Semi-major axis (AU).
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination (degrees).
    pub i_deg: f64,
    /// Argument of perihelion (degrees).
    pub omega_deg: f64,
    /// Longitude of ascending node (degrees).
    pub omega_big_deg: f64,
    /// Absolute magnitude H_r.
    pub h_mag: f64,
}

impl Bp519 {
    /// Keplerian elements (radians/AU), mean anomaly set to 0.
    pub fn elements(&self) -> OrbitalElements {
        OrbitalElements {
            a: self.a,
            e: self.e,
            i: self.i_deg * DEG2RAD,
            omega: self.omega_deg * DEG2RAD,
            omega_big: self.omega_big_deg * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    /// Semi-major axis as a typed [`Length`] (via the core
    /// [`OrbitalElements`] accessor).
    pub fn semi_major_axis(&self) -> Length {
        self.elements().semi_major_axis()
    }

    /// Inclination as a typed [`Angle`] (via the core [`OrbitalElements`]
    /// accessor).
    pub fn inclination(&self) -> Angle {
        self.elements().inclination()
    }

    /// Perihelion distance q = a(1 − e) as a typed [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        self.elements().perihelion_typed()
    }

    /// Aphelion distance Q = a(1 + e) as a typed [`Length`].
    pub fn aphelion_typed(&self) -> Length {
        self.elements().aphelion_typed()
    }

    /// Longitude of perihelion ϖ = ω + Ω as a typed [`Angle`], wrapped to
    /// [0, 2π).
    pub fn longitude_of_perihelion_typed(&self) -> Angle {
        radians(((self.omega_deg + self.omega_big_deg) * DEG2RAD).rem_euclid(TWO_PI))
    }
}

/// Becker et al. (2018) discovery-epoch solution. a ≈ 449 AU, e ≈ 0.92,
/// i ≈ 54°. ω/Ω are taken from the JPL fit (not quoted in the abstract).
pub const fn discovery_2018() -> Bp519 {
    Bp519 {
        source: "Becker et al. 2018 (arXiv:1805.05355)",
        a: 449.0,
        e: 0.92,
        i_deg: 54.0,
        omega_deg: 348.0,
        omega_big_deg: 135.0,
        h_mag: 4.3,
    }
}

/// Current JPL Small-Body Database osculating solution (768325, 2026-06-14).
pub const fn jpl_current() -> Bp519 {
    Bp519 {
        source: "JPL SBDB 768325 (lookup 2026-06-14)",
        a: 486.0,
        e: 0.927,
        i_deg: 54.1,
        omega_deg: 348.0,
        omega_big_deg: 135.0,
        h_mag: 4.32,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn discovery_elements_pinned() {
        use uom::si::length::astronomical_unit;
        let bp = discovery_2018();
        // Headline: highly inclined, eccentric, large-a detached ETNO.
        assert_relative_eq!(bp.a, 449.0, epsilon = 1.0);
        assert_relative_eq!(bp.e, 0.92, epsilon = 0.01);
        assert_relative_eq!(bp.i_deg, 54.0, epsilon = 1.0);
        // Detached: q well outside Neptune's reach (~35 AU).
        let q = bp.perihelion_typed().get::<astronomical_unit>();
        assert!((33.0..40.0).contains(&q), "q = {q:.1} AU");
        // Aphelion deep in the P9-clustered distance range.
        let aph = bp.aphelion_typed().get::<astronomical_unit>();
        assert!(aph > 800.0, "Q = {aph:.0} AU");
    }

    #[test]
    fn jpl_solution_drift_is_documented() {
        use uom::si::length::astronomical_unit;
        let disc = discovery_2018();
        let jpl = jpl_current();
        // i, q stable across solutions; a, e drift up the fit degeneracy.
        assert_relative_eq!(disc.i_deg, jpl.i_deg, epsilon = 0.5);
        assert_relative_eq!(
            disc.perihelion_typed().get::<astronomical_unit>(),
            jpl.perihelion_typed().get::<astronomical_unit>(),
            epsilon = 1.5
        );
        assert!(jpl.a > disc.a, "JPL a should drift upward");
    }

    #[test]
    fn elements_roundtrip_radians() {
        let bp = discovery_2018();
        let el = bp.elements();
        assert_relative_eq!(el.i, 54.0 * DEG2RAD, epsilon = 1e-12);
        assert_relative_eq!(el.a, 449.0, epsilon = 1e-12);
    }

    #[test]
    fn typed_accessors_match_f64_sources() {
        use uom::si::angle::radian;
        use uom::si::length::astronomical_unit;
        let bp = discovery_2018();
        assert_relative_eq!(
            bp.semi_major_axis().get::<astronomical_unit>(),
            bp.a,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            bp.inclination().get::<radian>(),
            bp.i_deg * DEG2RAD,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            bp.perihelion_typed().get::<astronomical_unit>(),
            bp.a * (1.0 - bp.e),
            epsilon = 1e-9
        );
        assert_relative_eq!(
            bp.aphelion_typed().get::<astronomical_unit>(),
            bp.a * (1.0 + bp.e),
            epsilon = 1e-9
        );
        assert_relative_eq!(
            bp.longitude_of_perihelion_typed().get::<radian>(),
            ((bp.omega_deg + bp.omega_big_deg) * DEG2RAD).rem_euclid(TWO_PI),
            epsilon = 1e-9
        );
    }
}
