//! Typed units of measure (`uom`) scaffold for the workspace.
//!
//! Astronomy-friendly quantity aliases, unit constructors, and typed physical
//! constants built on `uom`'s SI system. This is the foundation for migrating
//! the workspace off bare `f64` toward dimensionally-checked quantities.
//!
//! During the migration the existing `f64` constants and APIs in
//! [`crate::constants`] / [`crate::types`] are kept ALONGSIDE these typed
//! versions; nothing here changes existing behavior.
//!
//! Conventions:
//! - Lengths in astronomical units, times in days, masses in solar / Earth
//!   masses, angles in radians (with degree / arcsecond constructors). `uom`
//!   stores everything in SI base internally, so values are exact regardless of
//!   the unit you read them out in.
//! - `GM` (the standard gravitational parameter) has dimension L³·T⁻²; it has
//!   no named SI quantity, so [`GravitationalParameter`] is a typed alias.

use core::marker::PhantomData;

pub use uom::si::f64::{Angle, AngularVelocity, Length, Mass, Time, Velocity};
use uom::si::{Quantity, ISQ, SI};
use uom::typenum::{N2, P3, Z0};

use uom::si::angle::{degree, radian, second as arcsecond};
use uom::si::length::{astronomical_unit, kilometer, meter};
use uom::si::mass::kilogram;
use uom::si::time::day;

/// Standard gravitational parameter `GM`, dimension L³·T⁻² (no named SI
/// quantity exists, so it is spelled out via the ISQ dimension exponents:
/// length³, time⁻²).
pub type GravitationalParameter = Quantity<ISQ<P3, Z0, N2, Z0, Z0, Z0, Z0>, SI<f64>, f64>;

/// One solar mass in kilograms (IAU nominal `GM_sun`/`G`).
pub const SOLAR_MASS_KG: f64 = 1.988_409_87e30;
/// One Earth mass in kilograms.
pub const EARTH_MASS_KG: f64 = 5.972_19e24;
/// Days per Julian century.
pub const DAYS_PER_CENTURY: f64 = 36_525.0;

// ---- unit constructors (read with `.get::<unit>()`) ----

/// Length from astronomical units.
pub fn au(x: f64) -> Length {
    Length::new::<astronomical_unit>(x)
}
/// Length from kilometers.
pub fn km(x: f64) -> Length {
    Length::new::<kilometer>(x)
}
/// Length from meters.
pub fn meters(x: f64) -> Length {
    Length::new::<meter>(x)
}
/// Mass from solar masses.
pub fn solar_masses(x: f64) -> Mass {
    Mass::new::<kilogram>(x * SOLAR_MASS_KG)
}
/// Mass from Earth masses.
pub fn earth_masses(x: f64) -> Mass {
    Mass::new::<kilogram>(x * EARTH_MASS_KG)
}
/// Time from days.
pub fn days(x: f64) -> Time {
    Time::new::<day>(x)
}
/// Time from Julian centuries.
pub fn julian_centuries(x: f64) -> Time {
    Time::new::<day>(x * DAYS_PER_CENTURY)
}
/// Angle from radians.
pub fn radians(x: f64) -> Angle {
    Angle::new::<radian>(x)
}
/// Angle from degrees.
pub fn degrees(x: f64) -> Angle {
    Angle::new::<degree>(x)
}
/// Angle from arcseconds.
pub fn arcseconds(x: f64) -> Angle {
    Angle::new::<arcsecond>(x)
}

/// Construct a `GM` value from its SI base value (m³·s⁻²).
pub fn gm_si(value_m3_s2: f64) -> GravitationalParameter {
    GravitationalParameter {
        dimension: PhantomData,
        units: PhantomData,
        value: value_m3_s2,
    }
}

// ---- typed physical constants (mirror of crate::constants) ----

/// Heliocentric gravitational parameter `GM_sun` (IAU 2015 nominal,
/// 1.327124400e20 m³·s⁻²).
pub fn gm_sun() -> GravitationalParameter {
    gm_si(1.327_124_400_18e20)
}

/// One astronomical unit as a [`Length`].
pub fn astronomical_unit_length() -> Length {
    au(1.0)
}

// ---- typed planetary constants ----
//
// These are built from the existing `f64` constants in [`crate::constants`] so
// the numeric source of truth stays in one place; the test module pins the
// round-trip back to those values.

use crate::constants;

/// Convert a `GM` expressed in AU³·day⁻² (the workspace's native gravitational
/// unit) into a typed [`GravitationalParameter`] (SI base, m³·s⁻²).
pub fn gm_from_au3_day2(gm: f64) -> GravitationalParameter {
    let au_m = constants::AU_M;
    let day_s = constants::DAY_S;
    gm_si(gm * au_m * au_m * au_m / (day_s * day_s))
}

/// Heliocentric `GM_sun` as a typed quantity (from [`constants::GM_SUN`]).
pub fn gm_sun_typed() -> GravitationalParameter {
    gm_from_au3_day2(constants::GM_SUN)
}
/// Jovian `GM` (from [`constants::GM_JUPITER`]).
pub fn gm_jupiter() -> GravitationalParameter {
    gm_from_au3_day2(constants::GM_JUPITER)
}
/// Saturnian `GM` (from [`constants::GM_SATURN`]).
pub fn gm_saturn() -> GravitationalParameter {
    gm_from_au3_day2(constants::GM_SATURN)
}
/// Uranian `GM` (from [`constants::GM_URANUS`]).
pub fn gm_uranus() -> GravitationalParameter {
    gm_from_au3_day2(constants::GM_URANUS)
}
/// Neptunian `GM` (from [`constants::GM_NEPTUNE`]).
pub fn gm_neptune() -> GravitationalParameter {
    gm_from_au3_day2(constants::GM_NEPTUNE)
}

/// Neptune's semi-major axis as a [`Length`] (from [`constants::A_NEPTUNE_AU`]).
pub fn a_neptune() -> Length {
    au(constants::A_NEPTUNE_AU)
}
/// Solar radius as a [`Length`] (from [`constants::RADIUS_SUN_AU`]).
pub fn radius_sun() -> Length {
    au(constants::RADIUS_SUN_AU)
}
/// Earth equatorial radius as a [`Length`] (from [`constants::EARTH_RADIUS_KM`]).
pub fn radius_earth() -> Length {
    km(constants::EARTH_RADIUS_KM)
}

/// Neptune's orbital period as a [`Time`] (from [`constants::NEPTUNE_PERIOD_DAYS`]).
pub fn neptune_period() -> Time {
    days(constants::NEPTUNE_PERIOD_DAYS)
}
/// One Julian year as a [`Time`] (from [`constants::YEAR_DAYS`]).
pub fn julian_year() -> Time {
    days(constants::YEAR_DAYS)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn au_round_trips_to_meters() {
        assert_relative_eq!(au(1.0).get::<meter>(), 1.495_978_707e11, epsilon = 1.0);
        assert_relative_eq!(au(30.07).get::<astronomical_unit>(), 30.07, epsilon = 1e-12);
    }

    #[test]
    fn masses_construct_in_kg() {
        assert_relative_eq!(
            solar_masses(1.0).get::<kilogram>(),
            SOLAR_MASS_KG,
            epsilon = 1.0
        );
        assert_relative_eq!(
            earth_masses(6.2).get::<kilogram>(),
            6.2 * EARTH_MASS_KG,
            epsilon = 1.0
        );
        // Earth/Sun mass ratio ~ 3.0e-6.
        let ratio = earth_masses(1.0).get::<kilogram>() / solar_masses(1.0).get::<kilogram>();
        assert_relative_eq!(ratio, 3.003e-6, epsilon = 1e-8);
    }

    #[test]
    fn angles_and_times_convert() {
        assert_relative_eq!(
            degrees(180.0).get::<radian>(),
            std::f64::consts::PI,
            epsilon = 1e-12
        );
        // 1 arcsecond = 1/206264.806 rad.
        assert_relative_eq!(arcseconds(206_264.806).get::<radian>(), 1.0, epsilon = 1e-6);
        assert_relative_eq!(julian_centuries(1.0).get::<day>(), 36_525.0, epsilon = 1e-9);
    }

    #[test]
    fn gm_has_l3_t2_dimension_and_value() {
        assert_relative_eq!(gm_sun().value, 1.327_124_400_18e20, epsilon = 1e10);
        // Kepler: n = sqrt(GM / a^3) is an angular velocity (1/time). That this
        // type-checks proves the GM dimension (L³·T⁻²) is correct — a wrong
        // exponent would fail to unify with `AngularVelocity`.
        let a = au(1.0);
        let n: uom::si::f64::Frequency = (gm_sun() / (a * a * a)).sqrt();
        // Earth's mean motion ≈ 2π / 365.25 d ≈ 1.99e-7 s⁻¹.
        assert_relative_eq!(
            n.get::<uom::si::frequency::hertz>(),
            1.99e-7,
            epsilon = 1e-9
        );
    }

    #[test]
    fn typed_planetary_constants_round_trip() {
        // GM_sun typed (from AU³/day²) matches the literal SI nominal.
        assert_relative_eq!(gm_sun_typed().value, gm_sun().value, max_relative = 1e-3);
        // GM ordering Jupiter > Saturn > Neptune > Uranus is preserved.
        assert!(gm_jupiter().value > gm_saturn().value);
        assert!(gm_saturn().value > gm_neptune().value);
        assert!(gm_neptune().value > gm_uranus().value);
        // Neptune sits at ~30 AU and the Sun's radius is well under 0.01 AU.
        assert_relative_eq!(
            a_neptune().get::<astronomical_unit>(),
            constants::A_NEPTUNE_AU,
            epsilon = 1e-12
        );
        assert!(radius_sun().get::<astronomical_unit>() < 0.01);
        // Neptune's period is ~165 yr.
        assert_relative_eq!(
            (neptune_period() / julian_year()).value,
            164.8,
            epsilon = 0.2
        );
    }
}
