//! Reflex parallax of a distant body over a chosen epoch separation.
//!
//! All of the geometry is the baseline/distance relation already implemented in
//! `p9_core::coords::parallax`; this module only *reuses* it and re-expresses it in
//! the units the parallax search cares about (arcmin, and amplitude vs distance).
//!
//! Two quantities matter for the search:
//!
//! * the **full annual peak-to-peak** swing (Earth on opposite sides of the Sun,
//!   a 2 AU baseline) — `p9_core::coords::parallax::annual_parallax_arcsec`, and
//! * the **partial-baseline** shift over an arbitrary epoch separation. For two
//!   epochs separated by a fraction of Earth's orbit the projected baseline is
//!   `2 AU * sin(Δθ/2)` where `Δθ = 2π * Δt / 1 yr`; a ~6-month (quadrature)
//!   separation gives the full 2 AU, two consecutive nights give a tiny sliver.

use p9_core::constants::{AU_KM, PC_AU};
use p9_core::coords::parallax::{annual_parallax_arcsec, parallax_arcsec};

/// Arcseconds per arcminute.
pub const ARCSEC_PER_ARCMIN: f64 = 60.0;
/// Days in a (Julian) year, the Earth orbital period used for the baseline.
pub const YEAR_DAYS: f64 = p9_core::constants::YEAR_DAYS;

/// Reflex-parallax **half-angle** (the classic "parallax", arcsec) of a body at
/// `distance_au`: the angle Earth's 1 AU radius subtends from the body.
///
/// This is `(1 AU / d)` in radians times arcsec-per-radian, i.e. half the annual
/// peak-to-peak swing. It is what most people mean by "the parallax" of an
/// object and falls as `1/d`.
pub fn half_parallax_arcsec(distance_au: f64) -> f64 {
    // PC_AU == arcsec per radian (the parsec definition).
    PC_AU / distance_au
}

/// Reflex-parallax half-angle in **arcminutes** — the paper's order-of-arcmin
/// headline for a 500–700 AU body.
pub fn half_parallax_arcmin(distance_au: f64) -> f64 {
    half_parallax_arcsec(distance_au) / ARCSEC_PER_ARCMIN
}

/// Full annual peak-to-peak parallax in arcminutes (2 AU baseline), reusing
/// `p9_core::coords::parallax::annual_parallax_arcsec`.
pub fn annual_parallax_arcmin(distance_au: f64) -> f64 {
    annual_parallax_arcsec(distance_au) / ARCSEC_PER_ARCMIN
}

/// Projected Earth-orbit baseline (AU) between two epochs separated by
/// `epoch_separation_days`, for a circular 1 AU orbit. Two positions separated
/// by orbital phase `Δθ` are `2 * sin(Δθ/2)` AU apart (chord length). A
/// half-year separation (Δθ = π) gives the full 2 AU; a quarter year gives
/// `√2 ≈ 1.414` AU.
pub fn projected_baseline_au(epoch_separation_days: f64) -> f64 {
    let dtheta = std::f64::consts::TAU * epoch_separation_days / YEAR_DAYS;
    2.0 * (dtheta / 2.0).sin().abs()
}

/// Parallactic displacement (arcsec) of a body at `distance_au` observed at two
/// epochs separated by `epoch_separation_days`, using the chord baseline above.
/// This is the quantity the search actually measures between two visits.
pub fn epoch_parallax_arcsec(distance_au: f64, epoch_separation_days: f64) -> f64 {
    let baseline_km = projected_baseline_au(epoch_separation_days) * AU_KM;
    parallax_arcsec(baseline_km, distance_au)
}

/// Same, in arcminutes.
pub fn epoch_parallax_arcmin(distance_au: f64, epoch_separation_days: f64) -> f64 {
    epoch_parallax_arcsec(distance_au, epoch_separation_days) / ARCSEC_PER_ARCMIN
}

/// The amplitude-vs-distance relation as a sampled curve: for each distance in
/// `distances_au`, the annual half-parallax in arcmin. Useful for showing the
/// `1/d` falloff explicitly.
pub fn half_parallax_curve_arcmin(distances_au: &[f64]) -> Vec<(f64, f64)> {
    distances_au
        .iter()
        .map(|&d| (d, half_parallax_arcmin(d)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use approx::assert_relative_eq;

    #[test]
    fn half_parallax_is_arcminutes_at_p9_distance() {
        // ~7.5 arcmin half-angle at 460 AU (the paper's worked example).
        let p = half_parallax_arcmin(published::REFERENCE_DISTANCE_AU);
        assert_relative_eq!(p, published::HALF_PARALLAX_AT_REF_ARCMIN, epsilon = 0.05);
        // Order of arcminutes across the whole 500-700 AU prior.
        assert!(half_parallax_arcmin(500.0) > 5.0);
        assert!(half_parallax_arcmin(700.0) > 4.0);
        assert!(half_parallax_arcmin(700.0) < 6.0);
    }

    #[test]
    fn half_is_exactly_half_the_annual_swing() {
        let d = 600.0;
        assert_relative_eq!(
            annual_parallax_arcmin(d),
            2.0 * half_parallax_arcmin(d),
            epsilon = 1e-9
        );
    }

    #[test]
    fn parallax_falls_as_one_over_distance() {
        // Doubling the distance halves the parallax (inverse-linear).
        let p1 = half_parallax_arcmin(500.0);
        let p2 = half_parallax_arcmin(1000.0);
        assert_relative_eq!(p1 / p2, 2.0, epsilon = 1e-9);
        // d * parallax is constant (= 1 AU in arcsec * 1/60).
        let c500 = 500.0 * half_parallax_arcmin(500.0);
        let c700 = 700.0 * half_parallax_arcmin(700.0);
        assert_relative_eq!(c500, c700, epsilon = 1e-9);
        // That constant is PC_AU / 60 arcmin.
        assert_relative_eq!(c500, PC_AU / ARCSEC_PER_ARCMIN, epsilon = 1e-6);
    }

    #[test]
    fn quadrature_baseline_is_root_two_and_half_year_is_full() {
        // Quarter-year separation -> chord 2 sin(45 deg) = sqrt(2) AU.
        assert_relative_eq!(
            projected_baseline_au(YEAR_DAYS / 4.0),
            std::f64::consts::SQRT_2,
            epsilon = 1e-9
        );
        // Half-year separation -> full 2 AU diameter.
        assert_relative_eq!(projected_baseline_au(YEAR_DAYS / 2.0), 2.0, epsilon = 1e-9);
    }

    #[test]
    fn six_month_epoch_parallax_matches_annual_swing() {
        // A ~6-month baseline reaches the full peak-to-peak annual parallax.
        let d = published::REFERENCE_DISTANCE_AU;
        assert_relative_eq!(
            epoch_parallax_arcmin(d, YEAR_DAYS / 2.0),
            annual_parallax_arcmin(d),
            epsilon = 1e-9
        );
        // ~15 arcmin peak-to-peak at 460 AU.
        assert_relative_eq!(
            epoch_parallax_arcmin(d, YEAR_DAYS / 2.0),
            14.95,
            epsilon = 0.1
        );
    }

    #[test]
    fn consecutive_night_shift_is_still_arcseconds_for_p9() {
        // One-day baseline: chord 2 sin(pi/365.25) ~ 0.0172 AU -> at 460 AU the
        // shift is still ~7.7", a robustly measurable per-night motion that a
        // ~kpc field star (mas/yr proper motion) cannot mimic in a single night.
        let s = epoch_parallax_arcsec(published::REFERENCE_DISTANCE_AU, 1.0);
        assert!(s > 5.0 && s < 12.0, "one-night P9 shift = {s:.2} arcsec");
    }

    #[test]
    fn curve_is_monotone_decreasing() {
        let ds = [400.0, 500.0, 600.0, 700.0, 800.0];
        let curve = half_parallax_curve_arcmin(&ds);
        for w in curve.windows(2) {
            assert!(w[1].1 < w[0].1, "parallax must fall with distance");
        }
    }
}
