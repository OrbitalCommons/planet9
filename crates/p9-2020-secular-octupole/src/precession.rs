//! Apsidal precession rate and quadrupole-vs-octupole comparison.
//!
//! In coplanar secular theory the test particle's motion is governed by the
//! one-degree-of-freedom Hamiltonian H(Gamma, dvarpi), where the canonical
//! momentum conjugate to dvarpi is
//!
//!   Gamma = sqrt(GM_sun a) * (1 - sqrt(1 - e^2))
//!
//! (the Delaunay / AMD-like action; only the e-dependent factor matters for
//! the precession rate). The apsidal precession rate is
//!
//!   d(dvarpi)/dt = + dH / dGamma = (dH/de) (de/dGamma).
//!
//! At QUADRUPOLE order H depends on e only, so d(dvarpi)/dt is a fixed
//! (e-dependent) number independent of dvarpi: uniform circulation. The
//! OCTUPOLE term adds a dvarpi-dependence, so the precession rate is modulated
//! around the circle and, when the octupole is strong enough, the rate changes
//! sign within a dvarpi range — the signature of a libration island.
//!
//! Reference: Kohne & Batygin (2020); Murray & Dermott Ch. 7 for the canonical
//! secular variables.

use crate::octupole::{hamiltonian_octupole, hamiltonian_quadrupole, octupole_coupling};
use p9_core::constants::GM_SUN;
use p9_core::units::{days, radians, AngularVelocity};

/// Canonical apsidal momentum Gamma = sqrt(GM_sun a)(1 - sqrt(1-e^2)).
pub fn gamma_of_e(a: f64, e: f64) -> f64 {
    (GM_SUN * a).sqrt() * (1.0 - (1.0 - e * e).sqrt())
}

/// d e / d Gamma at fixed a. From Gamma above, dGamma/de = sqrt(GM a) * e /
/// sqrt(1-e^2), so de/dGamma = sqrt(1-e^2) / (sqrt(GM a) e).
fn de_dgamma(a: f64, e: f64) -> f64 {
    (1.0 - e * e).sqrt() / ((GM_SUN * a).sqrt() * e)
}

/// QUADRUPOLE apsidal precession rate d(dvarpi)/dt (radians/day). A function of
/// e only (dvarpi-independent): uniform circulation, no libration. Computed
/// from the e-derivative of the reused `p9_core` coplanar quadrupole.
pub fn precession_rate_quadrupole(a: f64, e: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let de = 1e-5;
    let dh_de = (hamiltonian_quadrupole(a, e + de, a_p, e_p, gm_p)
        - hamiltonian_quadrupole(a, e - de, a_p, e_p, gm_p))
        / (2.0 * de);
    dh_de * de_dgamma(a, e)
}

/// QUADRUPOLE apsidal precession rate as a typed [`AngularVelocity`] (the
/// dimension-checked view of [`precession_rate_quadrupole`], rad/day).
pub fn precession_rate_quadrupole_typed(
    a: f64,
    e: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> AngularVelocity {
    (radians(precession_rate_quadrupole(a, e, a_p, e_p, gm_p)) / days(1.0)).into()
}

/// OCTUPOLE apsidal precession rate d(dvarpi)/dt (radians/day) at apsidal angle
/// `dvarpi`. Carries the dvarpi-dependence from the octupole term; equals the
/// quadrupole rate plus an octupole modulation proportional to cos(dvarpi).
pub fn precession_rate_octupole(a: f64, e: f64, dvarpi: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let de = 1e-5;
    let dh_de = (hamiltonian_octupole(a, e + de, dvarpi, a_p, e_p, gm_p)
        - hamiltonian_octupole(a, e - de, dvarpi, a_p, e_p, gm_p))
        / (2.0 * de);
    dh_de * de_dgamma(a, e)
}

/// OCTUPOLE apsidal precession rate as a typed [`AngularVelocity`] (the
/// dimension-checked view of [`precession_rate_octupole`], rad/day).
pub fn precession_rate_octupole_typed(
    a: f64,
    e: f64,
    dvarpi: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> AngularVelocity {
    (radians(precession_rate_octupole(a, e, dvarpi, a_p, e_p, gm_p)) / days(1.0)).into()
}

/// Peak octupole modulation of the precession rate (radians/day): the
/// amplitude of the dvarpi-dependent part. Equals
/// C_oct e_p * de/dGamma (the e-derivative of +C_oct e e_p cos at the cos
/// extremum has magnitude C_oct e_p), a clean closed form we can pin.
pub fn octupole_modulation_amplitude(a: f64, e: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    let c_oct = octupole_coupling(a, a_p, e_p, gm_p);
    // d/de [ +C_oct e e_p cos(dvarpi) ] = C_oct e_p cos(dvarpi); peak |.| at
    // cos = +/-1 is C_oct e_p.
    c_oct * e_p * de_dgamma(a, e)
}

/// Peak octupole modulation of the precession rate as a typed
/// [`AngularVelocity`] (the dimension-checked view of
/// [`octupole_modulation_amplitude`], rad/day).
pub fn octupole_modulation_amplitude_typed(
    a: f64,
    e: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> AngularVelocity {
    (radians(octupole_modulation_amplitude(a, e, a_p, e_p, gm_p)) / days(1.0)).into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};
    use std::f64::consts::PI;
    use uom::si::angular_velocity::radian_per_second;

    fn gm_p9(mass_earth: f64) -> f64 {
        mass_earth * EARTH_MASS_SOLAR * GM_SUN
    }

    #[test]
    fn typed_precession_rates_match_f64() {
        use approx::assert_relative_eq;
        // 1 day = 86400 s, so a rate in rad/day reads back as rate/86400 in
        // rad/s through the typed AngularVelocity.
        const SEC_PER_DAY: f64 = 86_400.0;
        let gm = gm_p9(10.0);
        let (a, e, a_p, e_p, dv) = (250.0, 0.3, 700.0, 0.6, 0.4);

        assert_relative_eq!(
            precession_rate_quadrupole_typed(a, e, a_p, e_p, gm).get::<radian_per_second>(),
            precession_rate_quadrupole(a, e, a_p, e_p, gm) / SEC_PER_DAY,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            precession_rate_octupole_typed(a, e, dv, a_p, e_p, gm).get::<radian_per_second>(),
            precession_rate_octupole(a, e, dv, a_p, e_p, gm) / SEC_PER_DAY,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            octupole_modulation_amplitude_typed(a, e, a_p, e_p, gm).get::<radian_per_second>(),
            octupole_modulation_amplitude(a, e, a_p, e_p, gm) / SEC_PER_DAY,
            max_relative = 1e-12
        );
    }

    #[test]
    fn quadrupole_precession_is_apsidal_independent() {
        // By construction the quadrupole rate has no dvarpi argument; confirm
        // the octupole rate reduces to it when e_p = 0 (octupole off) at any
        // dvarpi.
        let gm = gm_p9(10.0);
        let q = precession_rate_quadrupole(250.0, 0.4, 700.0, 0.0, gm);
        let o0 = precession_rate_octupole(250.0, 0.4, 0.0, 700.0, 0.0, gm);
        let o1 = precession_rate_octupole(250.0, 0.4, PI, 700.0, 0.0, gm);
        assert!((q - o0).abs() < 1e-18);
        assert!((o0 - o1).abs() < 1e-18);
    }

    #[test]
    fn octupole_modulates_precession_with_apsidal_angle() {
        // With an eccentric perturber the octupole makes the precession rate
        // depend on dvarpi: aligned (0) and anti-aligned (pi) rates differ.
        let gm = gm_p9(10.0);
        let (a, a_p, e_p, e) = (250.0, 700.0, 0.6, 0.3);
        let r_aligned = precession_rate_octupole(a, e, 0.0, a_p, e_p, gm);
        let r_anti = precession_rate_octupole(a, e, PI, a_p, e_p, gm);
        assert!(
            (r_aligned - r_anti).abs() > 1e-12,
            "expected apsidal modulation, aligned {r_aligned:e} anti {r_anti:e}"
        );
    }

    #[test]
    fn modulation_amplitude_matches_aligned_minus_quadrupole() {
        // The octupole modulation amplitude equals half the aligned-vs-anti
        // precession-rate spread (the cos goes +1 -> -1).
        let gm = gm_p9(10.0);
        let (a, a_p, e_p, e) = (250.0, 700.0, 0.6, 0.3);
        let amp = octupole_modulation_amplitude(a, e, a_p, e_p, gm);
        let r_aligned = precession_rate_octupole(a, e, 0.0, a_p, e_p, gm);
        let r_anti = precession_rate_octupole(a, e, PI, a_p, e_p, gm);
        let half_spread = 0.5 * (r_aligned - r_anti).abs();
        let rel = ((amp - half_spread) / amp).abs();
        assert!(
            rel < 1e-3,
            "modulation amp {amp:e} vs half-spread {half_spread:e}"
        );
    }

    #[test]
    fn octupole_modulation_grows_with_strength() {
        // Modulation amplitude scales with e_p (octupole strength).
        let gm = gm_p9(10.0);
        let (a, a_p, e) = (250.0, 700.0, 0.3);
        let lo = octupole_modulation_amplitude(a, e, a_p, 0.3, gm);
        let hi = octupole_modulation_amplitude(a, e, a_p, 0.6, gm);
        assert!(hi > lo, "modulation hi {hi:e} not > lo {lo:e}");
    }
}
