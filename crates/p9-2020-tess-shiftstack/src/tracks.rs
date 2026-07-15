//! Number of trial orbital tracks for the TESS shift-stack search.
//!
//! Shift-and-stacking needs one co-add per hypothesised on-sky track. The grid
//! spacing is set by the orbit-parameter metric of `p9-2025-stacking` (a trial
//! that misses the true rate loses SNR quadratically), so the trial-track count
//! is *its* `n_trial_orbits` evaluated with TESS's PSF, baseline and the
//! distant-mover rate range. This module only supplies the TESS geometry; the
//! metric itself is not re-derived.
//!
//! TESS's coarse 21-arcsec pixels and short single-sector (27.4-day) baseline
//! make the rate cells *huge*, so the trial-track count stays modest — of
//! order 10³ over a single sector with the real ±50″/day parallax rate box
//! (Rice & Laughlin's grid used 748 shift vectors per sector), growing as T²
//! with more sectors — nothing like the billions a multi-year arcsec-PSF
//! survey (ZTF) demands.

use p9_core::analysis::stacking::orbit_metric::{
    apparent_sky_rate_at_opposition, n_trial_orbits, rate_resolution, snr_retained_exact,
};
use p9_core::types::P9Params;
use p9_core::units::{arcseconds, days, Angle, AngularVelocity};

use crate::tess;

/// On-sky rate range (arcsec/day) the search must span. The dominant apparent
/// motion is Earth's parallactic reflex, ~51·(70 AU/d) arcsec/day near
/// opposition: ≈45 arcsec/day for the closest (70 AU) targets Rice & Laughlin
/// consider, ≈8 arcsec/day at Planet Nine distances. Their actual grid spans
/// 4–58 px per ~27-day sector (≈44–52 arcsec/day peak rate, 748 shift
/// vectors). A ±50 arcsec/day box per axis covers it. (A previous version
/// used 5 arcsec/day — sized from the heliocentric mean-motion rate, 10×
/// under the real parallax scale — a box the true apparent motion of every
/// target, P9 included, falls OUTSIDE of.)
pub const RATE_RANGE_ARCSEC_PER_DAY: f64 = 50.0;

/// SNR tolerance per grid cell (retain ≥90% of the matched-filter SNR).
pub const SNR_TOLERANCE: f64 = 0.9;

/// Effective TESS PSF 1σ for track-grid purposes (≈0.5 px ≈ 10.5 arcsec).
pub fn track_psf_arcsec() -> f64 {
    tess::psf_sigma_arcsec(0.5)
}

/// Effective TESS track-grid PSF 1σ as a typed [`Angle`] (dimension-checked
/// view of [`track_psf_arcsec`]).
pub fn track_psf() -> Angle {
    arcseconds(track_psf_arcsec())
}

/// Number of trial tracks for a TESS shift-stack over `n_sectors` sectors,
/// from the shared orbit metric (TESS PSF, baseline, distant-mover rate box).
pub fn n_trial_tracks(n_sectors: f64) -> f64 {
    let baseline_days = tess::SECTOR_DAYS * n_sectors;
    n_trial_orbits(
        RATE_RANGE_ARCSEC_PER_DAY,
        track_psf_arcsec(),
        baseline_days,
        SNR_TOLERANCE,
    )
}

/// Rate-cell half-width (arcsec/day) for a TESS stack over `n_sectors`.
pub fn rate_cell_arcsec_per_day(n_sectors: f64) -> f64 {
    rate_resolution(
        track_psf_arcsec(),
        tess::SECTOR_DAYS * n_sectors,
        SNR_TOLERANCE,
    )
}

/// Rate-cell half-width as a typed [`AngularVelocity`] (dimension-checked view
/// of [`rate_cell_arcsec_per_day`], arcsec/day).
pub fn rate_cell(n_sectors: f64) -> AngularVelocity {
    (arcseconds(rate_cell_arcsec_per_day(n_sectors)) / days(1.0)).into()
}

/// Planet Nine's peak APPARENT on-sky rate (arcsec/day) for a model — the
/// parallax-dominated opposition rate from the shared metric, evaluated at
/// the orbit's perihelion (the fastest, closest configuration the box must
/// contain). This is the quantity to compare against the rate box; the
/// heliocentric mean-motion rate is ~20× smaller and sizing the box from it
/// is exactly the trap the shared metric now warns about.
pub fn p9_sky_rate(p9: &P9Params) -> f64 {
    apparent_sky_rate_at_opposition(p9.a * (1.0 - p9.e))
}

/// Planet Nine's peak apparent on-sky rate as a typed
/// [`AngularVelocity`] (dimension-checked view of [`p9_sky_rate`], arcsec/day).
pub fn p9_sky_rate_typed(p9: &P9Params) -> AngularVelocity {
    (arcseconds(p9_sky_rate(p9)) / days(1.0)).into()
}

/// Retained matched-filter SNR for a track whose rate is offset by `rate_error`
/// over an `n_sectors` baseline (exact Gaussian-overlap form from the metric).
pub fn snr_retained(n_sectors: f64, rate_error: f64) -> f64 {
    snr_retained_exact(
        track_psf_arcsec(),
        tess::SECTOR_DAYS * n_sectors,
        rate_error,
        2001,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn single_sector_track_count_is_finite_and_tiny() {
        // Coarse TESS pixels + 27-day baseline keep the count modest even
        // with the real parallax rate box, far from the billions a fine-PSF
        // multi-year survey needs.
        let n = n_trial_tracks(1.0);
        assert!(n.is_finite() && n > 0.0);
        // With the real ±50 arcsec/day parallax box a single sector needs
        // hundreds-to-thousands of tracks — the same order as Rice &
        // Laughlin's 748 shift vectors per sector.
        assert!((1.0e2..1.0e4).contains(&n), "single-sector tracks = {n:.0}");
        // A year of stacking grows as T² but stays far below the billions a
        // fine-PSF multi-year survey needs.
        let year = n_trial_tracks(tess::SECTORS_PER_YEAR);
        assert!((1.0e5..1.0e6).contains(&year), "year tracks = {year:.0}");
    }

    #[test]
    fn track_count_grows_as_baseline_squared() {
        // Two sectors quadruple the trial-track count (metric ∝ T²).
        let n1 = n_trial_tracks(1.0);
        let n2 = n_trial_tracks(2.0);
        assert_relative_eq!(n2 / n1, 4.0, epsilon = 1e-6);
        // … and a year grows it further still, monotonically.
        assert!(n_trial_tracks(tess::SECTORS_PER_YEAR) > n2);
    }

    #[test]
    fn rate_cells_sharpen_with_longer_baseline() {
        let c1 = rate_cell_arcsec_per_day(1.0);
        let c2 = rate_cell_arcsec_per_day(2.0);
        assert_relative_eq!(c1 / c2, 2.0, epsilon = 1e-9);
    }

    #[test]
    fn typed_track_quantities_match_f64() {
        use uom::si::angle::{radian, second as arcsecond};
        use uom::si::angular_velocity::radian_per_second;
        use uom::si::time::second;

        // PSF angle round-trips to its arcsec source.
        assert_relative_eq!(
            track_psf().get::<arcsecond>(),
            track_psf_arcsec(),
            epsilon = 1e-12
        );

        // A rate in arcsec/day expressed through the typed AngularVelocity
        // (read in rad/s) equals the source value converted by hand.
        let sec_per_day = days(1.0).get::<second>();
        let to_rad_per_s =
            |arcsec_per_day: f64| arcseconds(arcsec_per_day).get::<radian>() / sec_per_day;

        assert_relative_eq!(
            rate_cell(1.0).get::<radian_per_second>(),
            to_rad_per_s(rate_cell_arcsec_per_day(1.0)),
            max_relative = 1e-12
        );

        let p9 = P9Params::mcmc_2021();
        assert_relative_eq!(
            p9_sky_rate_typed(&p9).get::<radian_per_second>(),
            to_rad_per_s(p9_sky_rate(&p9)),
            max_relative = 1e-12
        );
    }

    #[test]
    fn planet_nine_sits_inside_the_tess_rate_box() {
        // The APPARENT (parallax-dominated) rate at the BB21 median orbit's
        // perihelion is ~12 arcsec/day — inside the ±50 arcsec/day box but
        // far outside the old 5 arcsec/day one.
        let rate = p9_sky_rate(&P9Params::mcmc_2021());
        assert!(
            rate > 5.0 && rate < RATE_RANGE_ARCSEC_PER_DAY,
            "P9 apparent rate = {rate:.1} arcsec/day"
        );
    }

    #[test]
    fn closest_targets_also_inside_the_box() {
        // Rice & Laughlin's nearest recovered targets sit at ~70 AU:
        // apparent rate ≈45 arcsec/day must be inside the box.
        let near = apparent_sky_rate_at_opposition(70.0);
        assert!(
            near > 40.0 && near < RATE_RANGE_ARCSEC_PER_DAY,
            "70 AU rate = {near:.1}"
        );
    }

    #[test]
    fn on_grid_track_keeps_most_snr() {
        // A track offset by one rate cell retains ≥ the chosen tolerance.
        let cell = rate_cell_arcsec_per_day(1.0);
        let retained = snr_retained(1.0, cell);
        assert!(retained >= SNR_TOLERANCE - 0.03, "retained = {retained:.3}");
        // A wildly wrong rate loses essentially all the signal.
        assert!(snr_retained(1.0, 50.0 * cell) < 0.2);
    }
}
