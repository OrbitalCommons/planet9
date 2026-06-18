//! Das & Batygin (2023) suppression criterion (arXiv:2307.00378).
//!
//! Das & Batygin show that the Madigan & McCourt inclination instability is
//! *suppressed* unless the disc self-gravity is strong enough to overcome the
//! competing differential apsidal/nodal precession imposed by the giant
//! planets. A coherent collective mode can only grow if all disc orbits
//! precess together; an external precession source (the giant planets, acting
//! as an orbit-averaged quadrupole / J2 ring) makes orbits at different
//! (a, e, i) precess at different rates, detuning the mode. The instability
//! survives only when the disc's own secular frequency dominates:
//!
//!   omega_disk(M_d, a)  >  g_planet(a, e, i).
//!
//! Setting them equal defines the critical disc mass M_crit: below it the mode
//! is stabilised, above it unstable. The giant-planet precession rate is the
//! standard secular J2-ring rate (Murray & Dermott 1999, Ch. 7), reusing
//! `p9_core::forces::j2_secular::combined_j2_jsu` (Jupiter+Saturn+Uranus) plus
//! Neptune for J2_eff:
//!
//!   g_planet = (3/2) n J2_eff / [ a^2 (1 - e^2)^2 ],
//!
//! where J2_eff = (1/2) sum_p (m_p/M_sun) a_p^2 has units of AU^2 and
//! 2 J2_eff/a^2 = sum_p (m_p/M_sun)(a_p/a)^2 is the usual dimensionless ring
//! strength. (This is the same J2_eff used by `p9-2024-oort-selfgrav`.)
//!
//! Equating omega_disk = (M_d/M_sun) n with g_planet and solving for the mass:
//!
//!   M_crit / M_sun = (3/2) J2_eff / [ a^2 (1 - e^2)^2 ].
//!
//! Note the disc mean motion n cancels: the criterion is a ratio of secular
//! rates, both proportional to n. M_crit is therefore independent of a's
//! orbital-time scaling and set by the ring geometry and the disc eccentricity.

use p9_core::constants::{A_NEPTUNE_AU, EARTH_MASS_SOLAR, MASS_NEPTUNE_SOLAR};
use p9_core::forces::j2_secular::{combined_j2_jsu, effective_j2};
use p9_core::units::{days, radians, solar_masses, AngularVelocity, Mass};

use crate::growth::{mean_motion, secular_frequency};

/// Effective J2 (AU^2) of the four giant planets: (1/2) sum_p (m_p/M_sun) a_p^2,
/// built from the shared `p9_core` J2 helpers (J+S+U) plus Neptune. Identical
/// construction to `p9-2024-oort-selfgrav`'s `compute_j2_effective`.
pub fn giant_planet_j2_eff() -> f64 {
    let (j2_jsu, _j4, _gm_boost) = combined_j2_jsu();
    j2_jsu + effective_j2(MASS_NEPTUNE_SOLAR, A_NEPTUNE_AU, 1.0)
}

/// Giant-planet secular apsidal precession rate for a disc orbit at (a, e)
/// [rad/day]:  g = (3/2) n J2_eff / [a^2 (1-e^2)^2].
pub fn giant_planet_precession_rate(a: f64, e: f64) -> f64 {
    let j2_eff = giant_planet_j2_eff();
    let eta2 = 1.0 - e * e;
    1.5 * mean_motion(a) * j2_eff / (a * a * eta2 * eta2)
}

/// Critical disc mass (solar masses) above which the inclination instability
/// is unstable per Das & Batygin (2023): set omega_disk = g_planet.
///
///   M_crit/M_sun = (3/2) J2_eff / [a^2 (1-e^2)^2].
pub fn critical_mass_solar(a: f64, e: f64) -> f64 {
    let j2_eff = giant_planet_j2_eff();
    let eta2 = 1.0 - e * e;
    1.5 * j2_eff / (a * a * eta2 * eta2)
}

/// Critical disc mass in Earth masses.
pub fn critical_mass_earth(a: f64, e: f64) -> f64 {
    critical_mass_solar(a, e) / EARTH_MASS_SOLAR
}

/// Giant-planet precession rate [`g`](giant_planet_precession_rate) as a typed
/// [`AngularVelocity`].
pub fn giant_planet_precession_rate_typed(a: f64, e: f64) -> AngularVelocity {
    (radians(giant_planet_precession_rate(a, e)) / days(1.0)).into()
}

/// Critical disc mass [`M_crit`](critical_mass_solar) as a typed [`Mass`].
pub fn critical_mass(a: f64, e: f64) -> Mass {
    solar_masses(critical_mass_solar(a, e))
}

/// Stability verdict for a disc of mass `m_disk_solar` at (a, e).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stability {
    /// Self-gravity dominates giant-planet precession: instability grows.
    Unstable,
    /// Giant-planet differential precession dominates: mode suppressed.
    Stable,
}

/// Apply the Das & Batygin suppression criterion: compare the disc self-gravity
/// secular frequency against the giant-planet precession rate at (a, e).
pub fn classify(m_disk_solar: f64, a: f64, e: f64) -> Stability {
    let omega_disk = secular_frequency(m_disk_solar, a);
    let g_planet = giant_planet_precession_rate(a, e);
    if omega_disk > g_planet {
        Stability::Unstable
    } else {
        Stability::Stable
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::units::SOLAR_MASS_KG;
    use uom::si::angular_velocity::radian_per_second;
    use uom::si::mass::kilogram;

    const A_AU: f64 = 250.0;
    const E_DISK: f64 = 0.6;

    #[test]
    fn typed_accessors_match_f64_sources() {
        assert_relative_eq!(
            giant_planet_precession_rate_typed(A_AU, E_DISK).get::<radian_per_second>(),
            giant_planet_precession_rate(A_AU, E_DISK) / 86_400.0,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            critical_mass(A_AU, E_DISK).get::<kilogram>(),
            critical_mass_solar(A_AU, E_DISK) * SOLAR_MASS_KG,
            epsilon = 1e10
        );
    }

    #[test]
    fn test_critical_mass_is_positive() {
        let m_crit = critical_mass_earth(A_AU, E_DISK);
        assert!(m_crit > 0.0, "critical mass must be positive, got {m_crit}");
        assert!(m_crit.is_finite());
    }

    #[test]
    fn test_critical_mass_published_order() {
        // Das & Batygin (2023) find the instability requires discs of several
        // to tens of Earth masses (consistent with Madigan & McCourt's massive
        // discs). The critical mass for a fiducial eccentric disc at a few
        // hundred AU should land in roughly the 0.1-100 M_earth band.
        let m_crit = critical_mass_earth(A_AU, E_DISK);
        assert!(
            (0.1..=100.0).contains(&m_crit),
            "critical mass {m_crit} M_earth outside the expected order"
        );
    }

    #[test]
    fn test_threshold_separates_stable_unstable() {
        let m_crit = critical_mass_solar(A_AU, E_DISK);
        // Just below: stable; just above: unstable.
        assert_eq!(classify(0.99 * m_crit, A_AU, E_DISK), Stability::Stable);
        assert_eq!(classify(1.01 * m_crit, A_AU, E_DISK), Stability::Unstable);
    }

    #[test]
    fn test_a_massive_disk_is_unstable_a_light_one_stable() {
        use p9_core::constants::EARTH_MASS_SOLAR;
        let m_crit_e = critical_mass_earth(A_AU, E_DISK);
        // A disc 10x above critical is unstable; one 10x below is stable.
        let heavy = 10.0 * m_crit_e * EARTH_MASS_SOLAR;
        let light = 0.1 * m_crit_e * EARTH_MASS_SOLAR;
        assert_eq!(classify(heavy, A_AU, E_DISK), Stability::Unstable);
        assert_eq!(classify(light, A_AU, E_DISK), Stability::Stable);
    }

    #[test]
    fn test_more_eccentric_disk_has_higher_critical_mass() {
        // (1-e^2)^-2 grows with e, so the giant-planet precession is faster
        // for eccentric orbits and a heavier disc is needed to overcome it.
        let m_low_e = critical_mass_earth(A_AU, 0.2);
        let m_high_e = critical_mass_earth(A_AU, 0.7);
        assert!(
            m_high_e > m_low_e,
            "eccentric disc should need more mass: {m_high_e} vs {m_low_e}"
        );
    }

    #[test]
    fn test_critical_mass_independent_of_orbital_time() {
        // M_crit ~ J2_eff / a^2: the n factors cancel, so M_crit depends on a
        // only through the 1/a^2 ring geometry, not through the period.
        let m1 = critical_mass_solar(A_AU, E_DISK);
        let m2 = critical_mass_solar(2.0 * A_AU, E_DISK);
        // 1/a^2 scaling => factor 4.
        assert!((m1 / m2 - 4.0).abs() < 1e-9, "ratio {}", m1 / m2);
    }
}
