//! TESS full-frame-image (FFI) cadence and survey-geometry bookkeeping.
//!
//! Rice & Laughlin (2020) co-add TESS FFIs along trial orbital tracks. The
//! depth gain of that co-add is set entirely by the *number of frames* that
//! land on the moving target over an observing baseline, so the first thing to
//! get right is how many FFIs TESS records per sector and per year.
//!
//! TESS observes the sky in 27.4-day "sectors". During the prime mission
//! (Sectors 1–26, which include the Sectors 18–19 that Rice & Laughlin search)
//! the full frame images were read out on a **30-minute cadence**. Each sector
//! spans two 13.7-day spacecraft orbits with a brief perigee data downlink gap,
//! and TESS completes one ecliptic-pole-to-equator strip of 13 sectors per
//! year. The large 21 arcsec pixels and ~1-pixel PSF make TESS shallow in any
//! single frame (T ≈ 15) but extremely well sampled in time — the regime where
//! shift-and-stacking pays off.

use p9_core::units::{arcseconds, Angle};

/// TESS full-frame-image cadence during the prime mission, in minutes
/// (Sectors 1–26 read the FFIs out every 30 minutes).
pub const FFI_CADENCE_MIN: f64 = 30.0;

/// Length of one TESS observing sector, in days (two ~13.7-day orbits).
pub const SECTOR_DAYS: f64 = 27.4;

/// Number of sectors TESS tiles per year (one hemisphere strip of 13).
pub const SECTORS_PER_YEAR: f64 = 13.0;

/// TESS detector pixel scale, in arcsec per pixel.
pub const PIXEL_SCALE_ARCSEC: f64 = 21.0;

/// Fraction of a sector lost to perigee downlinks / instrument overheads.
/// Rice & Laughlin retain the vast majority of the cadence; we adopt a modest
/// 10% loss so the frame count is a realistic *usable* count, not the naive
/// duty-cycle-1 maximum.
pub const SECTOR_DUTY_CYCLE: f64 = 0.90;

/// Number of usable FFI frames recorded over `n_sectors` sectors at the prime-
/// mission 30-minute cadence (with the [`SECTOR_DUTY_CYCLE`] gap allowance).
pub fn frames_in_sectors(n_sectors: f64) -> f64 {
    let minutes_per_sector = SECTOR_DAYS * 24.0 * 60.0;
    (minutes_per_sector / FFI_CADENCE_MIN) * SECTOR_DUTY_CYCLE * n_sectors
}

/// FFI frames recorded over a single sector — the natural shift-stack baseline
/// for a slow distant mover that stays in one sector's field.
pub fn frames_per_sector() -> f64 {
    frames_in_sectors(1.0)
}

/// FFI frames recorded over a full year of TESS observing (13 sectors).
pub fn frames_per_year() -> f64 {
    frames_in_sectors(SECTORS_PER_YEAR)
}

/// Effective Gaussian PSF 1σ width on the sky, in arcsec. The TESS PSF is
/// roughly one 21-arcsec pixel across (FWHM ≈ 1–2 px); a 1σ ≈ 0.5–1 px is a
/// fair matched-template width. Returned in arcsec from a pixel count.
pub fn psf_sigma_arcsec(sigma_pixels: f64) -> f64 {
    sigma_pixels * PIXEL_SCALE_ARCSEC
}

/// Effective Gaussian PSF 1σ width as a typed [`Angle`] (the dimension-checked
/// view of [`psf_sigma_arcsec`]).
pub fn psf_sigma(sigma_pixels: f64) -> Angle {
    arcseconds(psf_sigma_arcsec(sigma_pixels))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn one_sector_holds_roughly_a_thousand_frames() {
        // 27.4 d × 48 frames/d × 0.9 ≈ 1.18e3 usable FFIs per sector.
        let n = frames_per_sector();
        assert!((1.0e3..1.4e3).contains(&n), "frames/sector = {n:.0}");
    }

    #[test]
    fn two_sectors_double_the_frame_count() {
        // Rice & Laughlin stack the Sectors 18+19 pair.
        assert_relative_eq!(
            frames_in_sectors(2.0),
            2.0 * frames_per_sector(),
            epsilon = 1e-9
        );
    }

    #[test]
    fn a_year_holds_tens_of_thousands_of_frames() {
        let n = frames_per_year();
        assert_relative_eq!(n, 13.0 * frames_per_sector(), epsilon = 1e-6);
        assert!((1.0e4..2.0e4).contains(&n), "frames/year = {n:.0}");
    }

    #[test]
    fn psf_sigma_scales_with_pixel_scale() {
        assert_relative_eq!(psf_sigma_arcsec(1.0), PIXEL_SCALE_ARCSEC, epsilon = 1e-12);
        assert_relative_eq!(psf_sigma_arcsec(0.5), 10.5, epsilon = 1e-12);
    }

    #[test]
    fn typed_psf_sigma_matches_f64() {
        use uom::si::angle::second as arcsecond;
        assert_relative_eq!(
            psf_sigma(0.5).get::<arcsecond>(),
            psf_sigma_arcsec(0.5),
            epsilon = 1e-12
        );
    }
}
