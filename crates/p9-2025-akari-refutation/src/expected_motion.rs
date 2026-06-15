//! Bound-object half-year parallax and proper motion (Chen et al. 2025,
//! §3.2; their equations 4 and 7 of Cowan, Holder & Kaib 2016).
//!
//! This is the "consistent with a bound Planet Nine" envelope on the
//! six-month timescale that motivates the Chen et al. AKARI-only search.
//! Reproduced here from first principles so the published 10–25′ parallax /
//! 0.6–3.2′ proper-motion ranges can be pinned against real geometry rather
//! than copied.
//!
//! ## Parallax
//!
//! A bound object at heliocentric distance `d` (AU) has heliocentric
//! parallax `π = 1/d` rad (half-amplitude). Over six months Earth moves to
//! the opposite side of its orbit, so the apparent position swings by the
//! full *diameter* `2/d` rad. The proper motion is the heliocentric
//! transverse orbital rate `µ = √(GM·a(1−e²))/r²` integrated over the
//! interval — reused from `p9_core::coords::candidate_pair`.

use p9_core::constants::RAD2DEG;
use p9_core::coords::candidate_pair::transverse_rate_arcmin_yr;

/// Arcminutes per radian.
const ARCMIN_PER_RAD: f64 = RAD2DEG * 60.0;

/// Half a Julian year in years (the parallax/proper-motion interval used by
/// Chen et al. for the AKARI six-month confirmation cadence).
pub const HALF_YEAR: f64 = 0.5;

/// Heliocentric parallax *half-amplitude* `1/d` (arcmin) for a bound object
/// at distance `d_au`. At d = 600 AU this is `1/600 rad ≈ 5.73′`.
pub fn parallax_half_amplitude_arcmin(d_au: f64) -> f64 {
    (1.0 / d_au) * ARCMIN_PER_RAD
}

/// Six-month parallactic swing (arcmin): Earth crosses the full `2 AU`
/// diameter of its orbit, so the apparent shift is `2/d` rad. This is the
/// quantity Chen et al. quote as "10 to 25 arcminutes".
pub fn parallax_six_month_arcmin(d_au: f64) -> f64 {
    2.0 * parallax_half_amplitude_arcmin(d_au)
}

/// Six-month proper motion (arcmin) for an object at heliocentric distance
/// `r_au` on an orbit (`a_au`, `e`), from the heliocentric transverse rate.
/// For a circular orbit pass `a_au == r_au`, `e == 0`.
pub fn proper_motion_six_month_arcmin(r_au: f64, a_au: f64, e: f64) -> f64 {
    transverse_rate_arcmin_yr(r_au, a_au, e) * HALF_YEAR
}

/// Circular-orbit six-month proper motion (arcmin) at distance `d_au`.
pub fn proper_motion_six_month_circular_arcmin(d_au: f64) -> f64 {
    proper_motion_six_month_arcmin(d_au, d_au, 0.0)
}

/// The six-month parallax range (arcmin) over a distance interval
/// `[d_near, d_far]`. Parallax is monotonically decreasing in distance, so
/// the maximum is at `d_near` and the minimum at `d_far`.
///
/// Returns `(min, max)`.
pub fn parallax_range_arcmin(d_near: f64, d_far: f64) -> (f64, f64) {
    (
        parallax_six_month_arcmin(d_far),
        parallax_six_month_arcmin(d_near),
    )
}

/// The six-month proper-motion range (arcmin) over a distance interval,
/// scanning circular orbits at the far edge (slowest) up to the perihelion
/// of an `e = e_max` orbit at the near edge (fastest). At perihelion
/// `r = a(1−e)`, so `a = d_near/(1−e_max)`.
///
/// Returns `(min, max)`.
pub fn proper_motion_range_arcmin(d_near: f64, d_far: f64, e_max: f64) -> (f64, f64) {
    let slow = proper_motion_six_month_circular_arcmin(d_far);
    let a_peri = d_near / (1.0 - e_max);
    let fast = proper_motion_six_month_arcmin(d_near, a_peri, e_max);
    (slow, fast)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use approx::assert_relative_eq;

    #[test]
    fn parallax_half_amplitude_at_600au_is_5_73_arcmin() {
        // Headline arithmetic: 1 AU / 600 AU = 1/600 rad ≈ 5.73'.
        let amp = parallax_half_amplitude_arcmin(600.0);
        assert_relative_eq!(amp, 5.7296, epsilon = 1e-3);
    }

    #[test]
    fn parallax_six_month_swing_at_600au_is_about_11_5_arcmin() {
        let swing = parallax_six_month_arcmin(600.0);
        assert_relative_eq!(swing, 11.459, epsilon = 1e-2);
    }

    #[test]
    fn six_month_parallax_spans_chen_range() {
        // Chen et al. quote 10-25' for the parallax over six months. The
        // 25' end maps to ~300 AU (near edge), the 10' end to ~700 AU.
        let near = parallax_six_month_arcmin(300.0);
        let far = parallax_six_month_arcmin(700.0);
        assert!(
            near >= published::PARALLAX_6MO_MAX_ARCMIN - 3.0,
            "near-edge parallax = {near:.1}', expected ~25'"
        );
        assert!(
            far <= published::PARALLAX_6MO_MAX_ARCMIN
                && far >= published::PARALLAX_6MO_MIN_ARCMIN - 1.0,
            "far-edge parallax = {far:.1}', expected ~10'"
        );
    }

    #[test]
    fn parallax_range_brackets_chen_10_to_25() {
        // The 10-25' band is reproduced by the ~280-700 AU interval.
        let (lo, hi) = parallax_range_arcmin(286.0, 700.0);
        assert!(
            (lo - published::PARALLAX_6MO_MIN_ARCMIN).abs() < 1.5,
            "low end = {lo:.1}', expected ~10'"
        );
        assert!(
            (hi - published::PARALLAX_6MO_MAX_ARCMIN).abs() < 1.5,
            "high end = {hi:.1}', expected ~25'"
        );
    }

    #[test]
    fn six_month_proper_motion_spans_chen_range() {
        // Chen et al. quote 0.6-3.2' for the proper motion over six months.
        // Circular at 700 AU gives the slow (~0.6') end; perihelion of an
        // eccentric orbit near the inner edge of the search (~270 AU) gives
        // the fast (~3.2') end.
        let slow = proper_motion_six_month_circular_arcmin(700.0);
        assert!(
            (slow - published::PROPER_MOTION_6MO_MIN_ARCMIN).abs() < 0.15,
            "circular pm at 700 AU = {slow:.2}', expected ~0.6'"
        );
        // Perihelion of an e = 0.6 orbit at r = 270 AU (a = 675 AU): this is
        // the inner, fastest configuration of the search band.
        let a = 270.0 / (1.0 - 0.6);
        let fast = proper_motion_six_month_arcmin(270.0, a, 0.6);
        assert!(
            (fast - published::PROPER_MOTION_6MO_MAX_ARCMIN).abs() < 0.2,
            "perihelion pm at 270 AU (e=0.6) = {fast:.2}', expected ~3.2'"
        );
    }

    #[test]
    fn parallax_dominates_proper_motion_over_six_months() {
        // The central Chen et al. claim: on the half-year timescale parallax
        // dominates, so a bound P9 is found by its parallax, not its drift.
        for &d in &[400.0, 500.0, 600.0, 700.0] {
            let par = parallax_six_month_arcmin(d);
            let pm = proper_motion_six_month_circular_arcmin(d);
            assert!(
                par > 4.0 * pm,
                "at {d} AU parallax {par:.2}' should dominate pm {pm:.2}' by >4x"
            );
        }
    }

    #[test]
    fn proper_motion_range_over_search_band_is_subdominant() {
        let (pm_lo, pm_hi) =
            proper_motion_range_arcmin(published::DISTANCE_MIN_AU, published::DISTANCE_MAX_AU, 0.6);
        let (par_lo, _) =
            parallax_range_arcmin(published::DISTANCE_MIN_AU, published::DISTANCE_MAX_AU);
        // Even the fastest proper motion in the 500-700 AU band stays below
        // the smallest parallax there.
        assert!(
            pm_hi < par_lo,
            "max pm {pm_hi:.2}' should be < min parallax {par_lo:.2}'"
        );
        assert!(pm_lo > 0.0 && pm_lo < pm_hi);
    }
}
