//! Film export for `p9-2025-parallax-search`: the numbers its scene and ledger entry draw.

use p9_2025_parallax_search::discrimination::{contrast_ratio, max_tolerable_error_arcsec};
use p9_2025_parallax_search::parallax::{
    epoch_parallax_arcsec, half_parallax_arcmin, projected_baseline_au,
};
use p9_2025_parallax_search::published;
use p9_core::analysis::photometry::SOLAR_V_MINUS_R;
use p9_core::analysis::stacking::orbit_metric::apparent_sky_rate_at_opposition;
use p9_core::constants::YEAR_DAYS;
use p9_core::coords::candidate_pair::annual_proper_motion_circular;
use p9_core::data::reference_population::generate_reference_population;
use rand::SeedableRng;
use serde_json::{Value, json};

/// Size of the reference population whose brightness is compared with the
/// search depth.
const N_POP: usize = 3000;
/// Nights between the two visits of the 2023 campaign (24 and 25 September).
const NIGHTS_APART: f64 = 1.0;
/// Whole sky in square degrees.
const SKY_DEG2: f64 = 41_252.96;

/// Figures quoted by Socas-Navarro & Trujillo (2025) that the crate does not
/// carry: the nightly displacement range, the parallax-to-orbital-motion
/// ratio, the depth range and the field dimensions.
const PUBLISHED_NIGHTLY_ARCSEC: (f64, f64) = (4.0, 7.0);
const PUBLISHED_RATIO: (f64, f64) = (20.0, 30.0);
const PUBLISHED_DEPTH_R: (f64, f64) = (21.0, 21.4);
const PUBLISHED_FIELD_DEG: (f64, f64) = (17.9, 5.5);
const PUBLISHED_SEEING_ARCSEC: f64 = 1.0;

/// Orbital (heliocentric, circular) drift over `days`, in arcsec.
fn orbital_drift_arcsec(distance_au: f64, days: f64) -> f64 {
    annual_proper_motion_circular(distance_au) * 60.0 * days / YEAR_DAYS
}

pub fn export() -> Value {
    let distances: Vec<f64> = (0..=70).map(|k| 300.0 + 10.0 * k as f64).collect();
    let parallax: Vec<f64> = distances
        .iter()
        .map(|&d| epoch_parallax_arcsec(d, NIGHTS_APART))
        .collect();
    let orbital: Vec<f64> = distances
        .iter()
        .map(|&d| orbital_drift_arcsec(d, NIGHTS_APART))
        .collect();
    let net: Vec<f64> = distances
        .iter()
        .map(|&d| apparent_sky_rate_at_opposition(d) * NIGHTS_APART)
        .collect();
    let ratio: Vec<f64> = parallax.iter().zip(&orbital).map(|(p, o)| p / o).collect();

    // Distances at which the net nightly shift equals the published bounds.
    let distance_for = |shift: f64| {
        let (mut lo, mut hi) = (100.0_f64, 3000.0_f64);
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            if apparent_sky_rate_at_opposition(mid) * NIGHTS_APART > shift {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    };

    // How the shift builds up with the time between visits.
    let days: Vec<f64> = (0..=90).map(|k| 2.0 * k as f64 + 1.0).collect();
    let growth_500: Vec<f64> = days
        .iter()
        .map(|&t| epoch_parallax_arcsec(500.0, t) / 60.0)
        .collect();

    // Predicted brightness against the depth reached.
    let mut rng = rand::rngs::StdRng::seed_from_u64(2025);
    let r_mags: Vec<f64> = generate_reference_population(N_POP, &mut rng)
        .iter()
        .map(|p| p.v_magnitude - SOLAR_V_MINUS_R)
        .collect();
    let bright = r_mags
        .iter()
        .filter(|&&m| m < published::LIMITING_MAG_R)
        .count() as f64
        / N_POP as f64;

    let snr = 5.0;
    json!({
        "nights_apart": NIGHTS_APART,
        "baseline_au": projected_baseline_au(NIGHTS_APART),
        "distance_au": distances,
        "parallax_arcsec": parallax,
        "orbital_arcsec": orbital,
        "net_arcsec": net,
        "ratio": ratio,
        "shift_500": apparent_sky_rate_at_opposition(500.0) * NIGHTS_APART,
        "shift_700": apparent_sky_rate_at_opposition(700.0) * NIGHTS_APART,
        "ratio_500": epoch_parallax_arcsec(500.0, NIGHTS_APART)
            / orbital_drift_arcsec(500.0, NIGHTS_APART),
        "distance_for_shift": [
            distance_for(PUBLISHED_NIGHTLY_ARCSEC.1),
            distance_for(PUBLISHED_NIGHTLY_ARCSEC.0),
        ],
        "growth": {"days": days, "arcmin_500": growth_500},
        "half_parallax_ref_arcmin": half_parallax_arcmin(published::REFERENCE_DISTANCE_AU),
        "star_contrast_500": contrast_ratio(
            500.0,
            published::TYPICAL_FIELD_STAR_PARALLAX_ARCSEC,
        ),
        "tolerable_error_arcsec": {
            "snr": snr,
            "at_500": max_tolerable_error_arcsec(500.0, snr, NIGHTS_APART),
            "at_700": max_tolerable_error_arcsec(700.0, snr, NIGHTS_APART),
        },
        "r_mags": r_mags,
        "bright_fraction": bright,
        "sky_fraction": published::SURVEY_AREA_DEG2 / SKY_DEG2,
        "published": {
            "area_deg2": published::SURVEY_AREA_DEG2,
            "field_deg": [PUBLISHED_FIELD_DEG.0, PUBLISHED_FIELD_DEG.1],
            "depth_r": published::LIMITING_MAG_R,
            "depth_r_range": [PUBLISHED_DEPTH_R.0, PUBLISHED_DEPTH_R.1],
            "nightly_arcsec": [PUBLISHED_NIGHTLY_ARCSEC.0, PUBLISHED_NIGHTLY_ARCSEC.1],
            "ratio": [PUBLISHED_RATIO.0, PUBLISHED_RATIO.1],
            "seeing_arcsec": PUBLISHED_SEEING_ARCSEC,
            "reference_distance_au": published::REFERENCE_DISTANCE_AU,
            "half_parallax_ref_arcmin": published::HALF_PARALLAX_AT_REF_ARCMIN,
        },
    })
}
