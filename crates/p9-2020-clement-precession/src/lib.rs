//! Reproduction of Clement & Kaib (2020), "Orbital Precession in the
//! Distant Solar System" (arXiv:2005.05326), and the companion 2021 study of
//! Neptune's distant mean-motion resonances under a Planet Nine
//! (arXiv:2105.01065).
//!
//! Two headline results are reproduced from real computations (no hard-coded
//! answers), reusing `p9-core`'s secular, resonance, constants, and ETNO data:
//!
//! 1. [`precession`] — the orbit-averaged quadrupole apsidal precession rate
//!    dϖ/dt imposed on a distant ETNO by the known giant planets, the extra
//!    rate from a nominal Planet Nine, and the combined rate. The
//!    giant-planet-only precession period is of order ~1 Gyr for the clustered
//!    ETNOs (slow, coherent circulation); an exterior P9 adds a prograde term
//!    that shortens it. The leading-order analytic rates are validated against
//!    a finite difference of `p9_core`'s exact Gauss-ring secular Hamiltonian.
//!
//! 2. [`resonance`] — the semi-major-axis locations of Neptune's distant
//!    exterior n:1 (n = 2..12) and n:2 resonances, computed from the shared
//!    `resonance_semi_major_axis` (a = a_N (n/m)^{2/3}), and the nearest
//!    distant resonance for each clustered ETNO in `p9_core::data::etno`.

pub mod precession;
pub mod resonance;

use p9_core::constants::GYR_DAYS;
use p9_core::units::{au, days, radians, AngularVelocity, Length, Time};

/// Convenience summary of the two headline results for a test ETNO and a
/// nominal Planet Nine. Returned as plain numbers so callers (and the binary
/// report below) can serialize or print them.
#[derive(Debug, Clone)]
pub struct Summary {
    /// Test ETNO semi-major axis (AU).
    pub a_au: f64,
    /// Test ETNO eccentricity.
    pub e: f64,
    /// Giant-planet-induced precession rate (rad/day).
    pub giant_rate: f64,
    /// Giant-planet precession period (Gyr).
    pub giant_period_gyr: f64,
    /// Combined giant-planet + P9 precession rate (rad/day).
    pub combined_rate: f64,
    /// Combined precession period (Gyr).
    pub combined_period_gyr: f64,
    /// Nearest-resonance matches for the clustered sample.
    pub matches: Vec<resonance::EtnoResonance>,
}

impl Summary {
    /// Test ETNO semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }

    /// Giant-planet-induced apsidal precession rate as a typed
    /// [`AngularVelocity`] (from rad/day).
    pub fn giant_precession_rate(&self) -> AngularVelocity {
        (radians(self.giant_rate) / days(1.0)).into()
    }

    /// Combined giant-planet + P9 apsidal precession rate as a typed
    /// [`AngularVelocity`] (from rad/day).
    pub fn combined_precession_rate(&self) -> AngularVelocity {
        (radians(self.combined_rate) / days(1.0)).into()
    }

    /// Giant-planet precession period as a typed [`Time`] (from Gyr).
    pub fn giant_period(&self) -> Time {
        days(self.giant_period_gyr * GYR_DAYS)
    }

    /// Combined precession period as a typed [`Time`] (from Gyr).
    pub fn combined_period(&self) -> Time {
        days(self.combined_period_gyr * GYR_DAYS)
    }
}

/// Build the [`Summary`] for an ETNO at (a, e) and a Planet Nine.
pub fn summarize(a_au: f64, e: f64, p9: &p9_core::types::P9Params) -> Summary {
    let giant_rate = precession::giant_planet_precession_rate(a_au, e);
    let combined_rate = precession::combined_precession_rate(a_au, e, p9);
    Summary {
        a_au,
        e,
        giant_rate,
        giant_period_gyr: precession::precession_period_gyr(giant_rate),
        combined_rate,
        combined_period_gyr: precession::precession_period_gyr(combined_rate),
        matches: resonance::match_sample_to_resonances(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::types::P9Params;

    #[test]
    fn test_summary_is_self_consistent() {
        // a = 250 AU, q = 40 AU ⇒ e = 0.84.
        let a = 250.0;
        let e = 1.0 - 40.0 / a;
        let s = summarize(a, e, &P9Params::nominal_2016());

        assert!(s.giant_rate > 0.0);
        assert!(s.combined_rate > s.giant_rate);
        assert!(s.combined_period_gyr < s.giant_period_gyr);
        assert_eq!(
            s.matches.len(),
            p9_core::data::etno::BROWN_2017_SAMPLE.len()
        );
    }

    #[test]
    fn test_typed_accessors_match_f64() {
        use approx::assert_relative_eq;
        use p9_core::constants::{DAY_S, GYR_DAYS};
        use uom::si::angular_velocity::radian_per_second;
        use uom::si::length::astronomical_unit;
        use uom::si::time::day;

        let a = 250.0;
        let e = 1.0 - 40.0 / a;
        let s = summarize(a, e, &P9Params::nominal_2016());
        assert_relative_eq!(
            s.semi_major_axis().get::<astronomical_unit>(),
            s.a_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            s.giant_precession_rate().get::<radian_per_second>(),
            s.giant_rate / DAY_S,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            s.combined_precession_rate().get::<radian_per_second>(),
            s.combined_rate / DAY_S,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            s.giant_period().get::<day>(),
            s.giant_period_gyr * GYR_DAYS,
            epsilon = 1e-3
        );
        assert_relative_eq!(
            s.combined_period().get::<day>(),
            s.combined_period_gyr * GYR_DAYS,
            epsilon = 1e-3
        );
    }
}
