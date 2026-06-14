//! Reproduction of Clement & Sheppard (2020), "Orbital Precession in the
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
}
