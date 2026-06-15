//! Apparent sky motion of Planet Nine between the IRAS (1983) and AKARI
//! (2006) epochs: heliocentric orbital motion plus geocentric parallax.
//!
//! The generic two-epoch candidate-pair geometry lives in
//! [`p9_core::coords::candidate_pair`]. The survey-window-bound ephemeris scans
//! below bind the IRAS/AKARI observing windows
//! (`survey_model::IrasSurvey`/`AkariFisSurvey`) to the core scans, keeping the
//! paper-specific epoch tables in this crate.

use p9_core::coords::candidate_pair;
use p9_core::coords::observer::{EarthProvider, Time, Timescale};

use crate::survey_model::{AkariFisSurvey, IrasSurvey};

/// The IRAS and AKARI observing windows that the survey-bound ephemeris
/// scans sample over.
fn survey_windows() -> ((Time, Time), (Time, Time)) {
    let ts = Timescale::default();
    (
        IrasSurvey::default().window(&ts),
        AkariFisSurvey::default().window(&ts),
    )
}

/// Range of apparent two-epoch separations (arcmin) achievable at
/// heliocentric distance `d_au` with real Earth ephemeris states, scanning
/// detection epochs inside the real IRAS and AKARI observing windows.
///
/// Returns `(sep_min, sep_max)`.
pub fn separation_window_ephemeris(
    earth: &mut impl EarthProvider,
    d_au: f64,
    e_max: f64,
) -> (f64, f64) {
    let (w1, w2) = survey_windows();
    candidate_pair::separation_window_ephemeris(earth, d_au, e_max, &w1, &w2)
}

/// Ephemeris-grade counterpart of `implied_distance` over the real IRAS and
/// AKARI observing windows.
pub fn implied_distance_ephemeris(earth: &mut impl EarthProvider, sep_arcmin: f64) -> f64 {
    let (w1, w2) = survey_windows();
    candidate_pair::implied_distance_ephemeris(earth, sep_arcmin, &w1, &w2)
}
