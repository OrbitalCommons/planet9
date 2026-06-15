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
//! make the rate cells *huge*, so the trial-track count is tiny — of order ten
//! over a single sector, ~1500 over a year — nothing like the billions a
//! multi-year arcsec-PSF survey (ZTF) demands. Stacking across more sectors
//! lengthens the baseline and grows the count as T².

use p9_core::analysis::stacking::orbit_metric::{
    n_trial_orbits, p9_orbital_sky_rate, rate_resolution, snr_retained_exact,
};
use p9_core::types::P9Params;

use crate::tess;

/// On-sky rate range (arcsec/day) the search must span. Distant TNOs / Planet
/// Nine move at well under an arcsec/day; the dominant apparent motion is
/// Earth's parallax, giving up to a few arcsec/day for the closest (70 au)
/// targets Rice & Laughlin consider. A 5 arcsec/day box per axis covers it.
pub const RATE_RANGE_ARCSEC_PER_DAY: f64 = 5.0;

/// SNR tolerance per grid cell (retain ≥90% of the matched-filter SNR).
pub const SNR_TOLERANCE: f64 = 0.9;

/// Effective TESS PSF 1σ for track-grid purposes (≈0.5 px ≈ 10.5 arcsec).
pub fn track_psf_arcsec() -> f64 {
    tess::psf_sigma_arcsec(0.5)
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

/// Planet Nine's heliocentric on-sky angular rate (arcsec/day) for a model —
/// reused from the shared metric to confirm P9 sits inside the search box.
pub fn p9_sky_rate(p9: &P9Params) -> f64 {
    p9_orbital_sky_rate(p9)
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
        // Coarse TESS pixels + 27-day baseline → only of order ten trial
        // tracks, far from the billions a fine-PSF multi-year survey needs.
        let n = n_trial_tracks(1.0);
        assert!(n.is_finite() && n > 0.0);
        assert!((1.0..1.0e2).contains(&n), "single-sector tracks = {n:.1}");
        // A year of stacking still stays in the modest thousands.
        let year = n_trial_tracks(tess::SECTORS_PER_YEAR);
        assert!((1.0e3..1.0e4).contains(&year), "year tracks = {year:.0}");
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
    fn planet_nine_sits_inside_the_tess_rate_box() {
        let rate = p9_sky_rate(&P9Params::mcmc_2021());
        assert!(
            rate > 0.0 && rate < RATE_RANGE_ARCSEC_PER_DAY,
            "P9 rate = {rate}"
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
