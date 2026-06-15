//! Apparent sky motion of Planet Nine between the IRAS (1983) and AKARI
//! (2006) epochs: heliocentric orbital motion plus geocentric parallax.
//!
//! The generic two-epoch candidate-pair geometry now lives in
//! [`p9_core::coords::candidate_pair`] and is re-exported here. The
//! survey-window-bound ephemeris scans below bind the IRAS/AKARI observing
//! windows (`survey_model::IrasSurvey`/`AkariFisSurvey`) to the core scans,
//! keeping the paper-specific epoch tables in this crate.

pub use p9_core::coords::candidate_pair::{
    annual_proper_motion_circular, apparent_separation_arcmin,
    apparent_separation_ephemeris_arcmin, circular_velocity, distance_range_for_separation,
    heliocentric_separation_arcmin, implied_distance, implied_distance_circular,
    orbital_period_years, proper_motion_table, separation_window, transverse_rate_arcmin_yr,
    transverse_rate_rad_day, vis_viva_velocity, E_MAX,
};

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
