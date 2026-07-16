//! IRAS and AKARI far-infrared survey models.
//!
//! IRAS (Infrared Astronomical Satellite, 1983):
//! - Observed at 12, 25, 60, and 100 µm
//! - The Faint Source Catalogue (FSC) contains ~173,000 sources at 60 µm
//! - Point source sensitivity ~0.2 Jy at 60 µm
//! - All-sky coverage ~96%
//! - Positional accuracy ~20" at 60 µm
//!
//! AKARI (2006–2007):
//! - Far-Infrared Surveyor (FIS) at 65, 90, 140, 160 µm
//! - The Monthly Unconfirmed Source List includes sources detected on
//!   hourly timescales but not confirmed after months — exactly the
//!   signature expected from a slow-moving solar system object
//! - AKARI-MUSL detection limit ~0.21 Jy at 90 µm (Phan et al. 2025
//!   Table 1 — the MUSL's "twice-deeper flux detection limit" is the
//!   paper's stated reason for using it; the Bright Source Catalogue
//!   limits are 0.55/0.44 Jy for v1/v2 and are kept as labelled
//!   constants below)
//! - All-sky coverage ~94%
//! - Positional accuracy ~30" at 90 µm

use p9_core::constants::YEAR_DAYS;
use p9_core::coords::observer::{Time, Timescale};
use p9_core::units::{self, arcseconds, julian_year, Angle};
use serde::{Deserialize, Serialize};

/// A far-infrared point source detection from either survey.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FirSource {
    /// Right ascension in degrees
    pub ra_deg: f64,
    /// Declination in degrees
    pub dec_deg: f64,
    /// Flux density in Janskys at the primary detection band
    /// (IRAS: 60 µm; AKARI: 90 µm)
    pub flux_jy: f64,
    /// Flux density in Janskys at the survey's secondary band, when
    /// catalogued (IRAS: 100 µm; AKARI: 65 µm). Drives the paper's
    /// within-survey colour cuts; `None` skips those cuts for this source.
    pub flux_secondary_jy: Option<f64>,
    /// Catalogue flux-quality acceptable (the paper rejects sources whose
    /// flux quality equals 1 in any of the 60/65/90/100 µm bands).
    pub flux_quality_ok: bool,
    /// Positional uncertainty in arcseconds
    pub pos_err_arcsec: f64,
    /// Survey epoch (fractional year, e.g., 1983.5 for IRAS)
    pub epoch_year: f64,
    /// Which survey this source comes from
    pub survey: Survey,
}

/// Survey identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Survey {
    Iras,
    Akari,
}

/// IRAS survey characteristics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IrasSurvey {
    /// Detection band wavelength in µm
    pub wavelength_um: f64,
    /// Point source sensitivity limit in Jy (Faint Source Catalogue)
    pub sensitivity_jy: f64,
    /// Sky coverage fraction
    pub sky_coverage: f64,
    /// Positional uncertainty in arcseconds
    pub position_error_arcsec: f64,
}

impl Default for IrasSurvey {
    fn default() -> Self {
        Self {
            wavelength_um: 60.0,
            sensitivity_jy: 0.2,
            sky_coverage: 0.96,
            position_error_arcsec: 20.0,
        }
    }
}

/// AKARI FIS survey characteristics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AkariFisSurvey {
    /// Detection band wavelength in µm
    pub wavelength_um: f64,
    /// Point source sensitivity limit in Jy (Monthly Unconfirmed Source List)
    pub sensitivity_jy: f64,
    /// Sky coverage fraction
    pub sky_coverage: f64,
    /// Positional uncertainty in arcseconds
    pub position_error_arcsec: f64,
}

/// AKARI-FIS Bright Source Catalogue v1 90 µm detection limit (Jy). NOT the
/// MUSL limit the paper's search uses — see `AkariFisSurvey::default`.
pub const AKARI_BSC_V1_LIMIT_JY: f64 = 0.55;
/// AKARI-FIS Bright Source Catalogue v2 90 µm detection limit (Jy).
pub const AKARI_BSC_V2_LIMIT_JY: f64 = 0.44;

impl Default for AkariFisSurvey {
    fn default() -> Self {
        Self {
            wavelength_um: 90.0,
            // AKARI-MUSL: ~0.21 Jy (Phan et al. 2025 Table 1). A previous
            // version used the BSC limit 0.55 Jy, making the 90 µm
            // detectability model ~2.6x too shallow and misrepresenting the
            // paper's premise (their Fig. 3 has even 7 M⊕ crossing both
            // limits near 500-520 AU).
            sensitivity_jy: 0.21,
            sky_coverage: 0.94,
            position_error_arcsec: 30.0,
        }
    }
}

impl IrasSurvey {
    /// Returns true if a source with the given flux would be catalogued.
    pub fn is_detectable(&self, flux_jy: f64) -> bool {
        flux_jy >= self.sensitivity_jy
    }

    /// Observing window of the IRAS all-sky survey: 1983-01-25 (launch)
    /// to 1983-11-21 (cryogen exhaustion). Catalogued sources were
    /// detected at unknown epochs inside this window.
    pub fn window(&self, ts: &Timescale) -> (Time, Time) {
        (ts.utc((1983, 1, 25)), ts.utc((1983, 11, 21)))
    }

    /// Survey mid-epoch, 1983-06-24 ≈ 1983.48 — the real midpoint of the
    /// observing window, refining the nominal "1983.5" quoted by
    /// Phan et al. (2025).
    pub fn mid_epoch(&self, ts: &Timescale) -> Time {
        ts.utc((1983, 6, 24))
    }
}

impl AkariFisSurvey {
    /// Returns true if a source with the given flux would be catalogued.
    pub fn is_detectable(&self, flux_jy: f64) -> bool {
        flux_jy >= self.sensitivity_jy
    }

    /// Observing window of the AKARI FIS cryogenic all-sky survey:
    /// 2006-05-08 to 2007-08-28 (liquid-helium boil-off). The Monthly
    /// Unconfirmed Source List samples this window.
    pub fn window(&self, ts: &Timescale) -> (Time, Time) {
        (ts.utc((2006, 5, 8)), ts.utc((2007, 8, 28)))
    }

    /// Survey mid-epoch, 2006-12-31 ≈ 2007.0 — the real midpoint of the
    /// cryogenic survey, refining the nominal "2006.5" quoted by
    /// Phan et al. (2025).
    pub fn mid_epoch(&self, ts: &Timescale) -> Time {
        ts.utc((2006, 12, 31))
    }
}

/// The epoch baseline between the IRAS and AKARI survey mid-epochs in
/// Julian years (≈ 23.5; the paper rounds to "23 years").
pub fn epoch_baseline(iras: &IrasSurvey, akari: &AkariFisSurvey) -> f64 {
    let ts = Timescale::default();
    (akari.mid_epoch(&ts).tdb() - iras.mid_epoch(&ts).tdb()) / YEAR_DAYS
}

/// The IRAS–AKARI epoch baseline as a dimension-checked [`units::Time`].
pub fn epoch_baseline_time(iras: &IrasSurvey, akari: &AkariFisSurvey) -> units::Time {
    julian_year() * epoch_baseline(iras, akari)
}

/// Angular separation between two sky positions in arcminutes.
///
/// Uses the Vincenty formula for great-circle distance, accurate at all
/// angular separations.
pub fn angular_separation_arcmin(ra1_deg: f64, dec1_deg: f64, ra2_deg: f64, dec2_deg: f64) -> f64 {
    let ra1 = ra1_deg.to_radians();
    let dec1 = dec1_deg.to_radians();
    let ra2 = ra2_deg.to_radians();
    let dec2 = dec2_deg.to_radians();

    let d_ra = ra2 - ra1;

    let num = ((dec2.cos() * d_ra.sin()).powi(2)
        + (dec1.cos() * dec2.sin() - dec1.sin() * dec2.cos() * d_ra.cos()).powi(2))
    .sqrt();
    let den = dec1.sin() * dec2.sin() + dec1.cos() * dec2.cos() * d_ra.cos();

    let angle_rad = num.atan2(den);
    angle_rad.to_degrees() * 60.0 // convert degrees to arcminutes
}

/// Angular separation between two sky positions as a dimension-checked
/// [`Angle`] (the typed sibling of [`angular_separation_arcmin`]).
pub fn angular_separation(ra1_deg: f64, dec1_deg: f64, ra2_deg: f64, dec2_deg: f64) -> Angle {
    arcseconds(angular_separation_arcmin(ra1_deg, dec1_deg, ra2_deg, dec2_deg) * 60.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::minute;
    use uom::si::time::day;

    #[test]
    fn typed_siblings_match_f64_sources() {
        let (iras, akari) = (IrasSurvey::default(), AkariFisSurvey::default());
        let yr_d = julian_year().get::<day>();
        assert_relative_eq!(
            epoch_baseline_time(&iras, &akari).get::<day>(),
            epoch_baseline(&iras, &akari) * yr_d,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            angular_separation(10.0, 20.0, 11.0, 21.0).get::<minute>(),
            angular_separation_arcmin(10.0, 20.0, 11.0, 21.0),
            max_relative = 1e-12
        );
    }

    #[test]
    fn iras_default_matches_paper() {
        let iras = IrasSurvey::default();
        assert!((iras.wavelength_um - 60.0).abs() < 1e-10);
        assert!((iras.sensitivity_jy - 0.2).abs() < 1e-10);
        assert!((iras.sky_coverage - 0.96).abs() < 1e-10);
    }

    #[test]
    fn akari_default_matches_paper() {
        let akari = AkariFisSurvey::default();
        assert!((akari.wavelength_um - 90.0).abs() < 1e-10);
        assert!((akari.sensitivity_jy - 0.21).abs() < 1e-10);
        assert!(akari.sensitivity_jy < AKARI_BSC_V2_LIMIT_JY);
    }

    #[test]
    fn epoch_baseline_is_about_23_years() {
        // Real mid-epoch baseline 1983-06-24 → 2006-12-31 is 23.52 yr
        // (the paper's "23 years" quote uses the nominal 1983.5/2006.5).
        let iras = IrasSurvey::default();
        let akari = AkariFisSurvey::default();
        let baseline = epoch_baseline(&iras, &akari);
        assert!(
            (baseline - 23.52).abs() < 0.05,
            "baseline = {baseline} years"
        );
    }

    #[test]
    fn mid_epochs_inside_observing_windows() {
        let ts = Timescale::default();
        let iras = IrasSurvey::default();
        let akari = AkariFisSurvey::default();
        for (mid, (start, end)) in [
            (iras.mid_epoch(&ts), iras.window(&ts)),
            (akari.mid_epoch(&ts), akari.window(&ts)),
        ] {
            assert!(start.tdb() < mid.tdb() && mid.tdb() < end.tdb());
            // Mid-epoch within a week of the window's midpoint.
            let center = 0.5 * (start.tdb() + end.tdb());
            assert!((mid.tdb() - center).abs() < 7.0);
        }
    }

    #[test]
    fn angular_separation_zero() {
        let sep = angular_separation_arcmin(180.0, 45.0, 180.0, 45.0);
        assert!(sep.abs() < 1e-10, "same position: sep = {sep}");
    }

    #[test]
    fn angular_separation_one_degree() {
        // Points separated by 1° in declination
        let sep = angular_separation_arcmin(180.0, 45.0, 180.0, 46.0);
        assert!(
            (sep - 60.0).abs() < 0.1,
            "1° dec separation = {sep}', expected 60'"
        );
    }

    #[test]
    fn angular_separation_symmetric() {
        let sep1 = angular_separation_arcmin(10.0, 20.0, 11.0, 21.0);
        let sep2 = angular_separation_arcmin(11.0, 21.0, 10.0, 20.0);
        assert!(
            (sep1 - sep2).abs() < 1e-10,
            "not symmetric: {sep1} vs {sep2}"
        );
    }

    #[test]
    fn iras_detection_threshold() {
        let iras = IrasSurvey::default();
        assert!(iras.is_detectable(0.3));
        assert!(iras.is_detectable(0.2));
        assert!(!iras.is_detectable(0.19));
    }

    #[test]
    fn akari_detection_threshold() {
        let akari = AkariFisSurvey::default();
        assert!(akari.is_detectable(0.6));
        assert!(akari.is_detectable(0.21));
        assert!(!akari.is_detectable(0.20));
    }
}
