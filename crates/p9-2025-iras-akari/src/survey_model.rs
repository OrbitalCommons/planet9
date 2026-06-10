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
//! - Point source sensitivity ~0.55 Jy at 90 µm
//! - All-sky coverage ~94%
//! - Positional accuracy ~30" at 90 µm

use serde::{Deserialize, Serialize};

/// A far-infrared point source detection from either survey.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FirSource {
    /// Right ascension in degrees
    pub ra_deg: f64,
    /// Declination in degrees
    pub dec_deg: f64,
    /// Flux density in Janskys at the primary detection band
    pub flux_jy: f64,
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
    /// Survey epoch (midpoint)
    pub epoch_year: f64,
}

impl Default for IrasSurvey {
    fn default() -> Self {
        Self {
            wavelength_um: 60.0,
            sensitivity_jy: 0.2,
            sky_coverage: 0.96,
            position_error_arcsec: 20.0,
            epoch_year: 1983.5,
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
    /// Survey epoch (midpoint)
    pub epoch_year: f64,
}

impl Default for AkariFisSurvey {
    fn default() -> Self {
        Self {
            wavelength_um: 90.0,
            sensitivity_jy: 0.55,
            sky_coverage: 0.94,
            position_error_arcsec: 30.0,
            epoch_year: 2006.5,
        }
    }
}

impl IrasSurvey {
    /// Returns true if a source with the given flux would be catalogued.
    pub fn is_detectable(&self, flux_jy: f64) -> bool {
        flux_jy >= self.sensitivity_jy
    }
}

impl AkariFisSurvey {
    /// Returns true if a source with the given flux would be catalogued.
    pub fn is_detectable(&self, flux_jy: f64) -> bool {
        flux_jy >= self.sensitivity_jy
    }
}

/// The epoch baseline between IRAS and AKARI in years.
pub fn epoch_baseline(iras: &IrasSurvey, akari: &AkariFisSurvey) -> f64 {
    akari.epoch_year - iras.epoch_year
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

#[cfg(test)]
mod tests {
    use super::*;

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
        assert!((akari.sensitivity_jy - 0.55).abs() < 1e-10);
    }

    #[test]
    fn epoch_baseline_is_23_years() {
        let iras = IrasSurvey::default();
        let akari = AkariFisSurvey::default();
        let baseline = epoch_baseline(&iras, &akari);
        assert!((baseline - 23.0).abs() < 0.5, "baseline = {baseline} years");
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
        assert!(akari.is_detectable(0.55));
        assert!(!akari.is_detectable(0.54));
    }
}
