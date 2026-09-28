//! Film export for `p9-2020-tess-shiftstack`: the numbers its scene and ledger entry draw.

use p9_2020_tess_shiftstack::depth::{
    p9_apparent_magnitude_at, sectors_to_reach_depth, stacked_depth_over_sectors,
};
use p9_2020_tess_shiftstack::published::{
    SENSITIVITY_DISTANCE_AU, SENSITIVITY_V_LIMIT, TESS_SINGLE_DEPTH, TESS_SINGLE_DEPTH_EFFECTIVE,
    V_BP519, V_SEDNA, V_TG422,
};
use p9_2020_tess_shiftstack::tess::{PIXEL_SCALE_ARCSEC, SECTOR_DAYS, frames_per_sector};
use p9_2020_tess_shiftstack::tracks::{n_trial_tracks, rate_cell_arcsec_per_day};
use p9_core::analysis::stacking::orbit_metric::apparent_sky_rate_at_opposition;
use p9_core::types::P9Params;
use serde_json::{Value, json};

/// Geometric albedo of the reflected-light comparison planet.
const ALBEDO: f64 = 0.4;

/// Distance range searched (AU) and the number of shift vectors per sector
/// (Rice & Laughlin 2020).
const SEARCH_DISTANCE_AU: (f64, f64) = (70.0, 800.0);
const PUBLISHED_SHIFT_VECTORS: f64 = 748.0;

pub fn export() -> Value {
    let p9 = P9Params::mcmc_2021();

    let log_sectors: Vec<f64> = (0..=60).map(|k| -2.0 + 3.2 * k as f64 / 60.0).collect();
    let depth: Vec<f64> = log_sectors
        .iter()
        .map(|&lg| stacked_depth_over_sectors(TESS_SINGLE_DEPTH_EFFECTIVE, 10f64.powf(lg)))
        .collect();

    let tnos: Vec<Value> = [
        ("Sedna", V_SEDNA),
        ("2015 BP519", V_BP519),
        ("2007 TG422", V_TG422),
    ]
    .iter()
    .map(|&(name, v)| {
        json!({
            "name": name,
            "v_mag": v,
            "sectors_needed": sectors_to_reach_depth(TESS_SINGLE_DEPTH_EFFECTIVE, v),
        })
    })
    .collect();

    let distances: Vec<f64> = (0..=73)
        .map(|k| SEARCH_DISTANCE_AU.0 + 10.0 * k as f64)
        .collect();
    let rate: Vec<f64> = distances
        .iter()
        .map(|&d| apparent_sky_rate_at_opposition(d))
        .collect();
    let drift_px: Vec<f64> = rate
        .iter()
        .map(|&r| r * SECTOR_DAYS / PIXEL_SCALE_ARCSEC)
        .collect();
    let v_p9: Vec<f64> = distances
        .iter()
        .map(|&d| p9_apparent_magnitude_at(&p9, ALBEDO, d))
        .collect();

    let perihelion = p9.a * (1.0 - p9.e);
    let aphelion = p9.a * (1.0 + p9.e);
    let one_sector = stacked_depth_over_sectors(TESS_SINGLE_DEPTH_EFFECTIVE, 1.0);

    json!({
        "frames_per_sector": frames_per_sector(),
        "pixel_scale_arcsec": PIXEL_SCALE_ARCSEC,
        "sector_days": SECTOR_DAYS,
        "single_depth_catalog": TESS_SINGLE_DEPTH,
        "single_depth_effective": TESS_SINGLE_DEPTH_EFFECTIVE,
        "depth_one_sector": one_sector,
        "depth_two_sectors": stacked_depth_over_sectors(TESS_SINGLE_DEPTH_EFFECTIVE, 2.0),
        "depth_vs_sectors": {"log_sectors": log_sectors, "depth": depth},
        "published_blind_limit_v": SENSITIVITY_V_LIMIT,
        "published_blind_distance_au": SENSITIVITY_DISTANCE_AU,
        "tnos": tnos,
        "distance_au": distances,
        "rate_arcsec_per_day": rate,
        "drift_px_per_sector": drift_px,
        "v_p9": v_p9,
        "p9": {
            "mass_earth": p9.mass_earth,
            "a_au": p9.a,
            "perihelion_au": perihelion,
            "aphelion_au": aphelion,
            "v_perihelion": p9_apparent_magnitude_at(&p9, ALBEDO, perihelion),
            "v_semi_major": p9_apparent_magnitude_at(&p9, ALBEDO, p9.a),
            "v_aphelion": p9_apparent_magnitude_at(&p9, ALBEDO, aphelion),
            "sectors_needed_at_a": sectors_to_reach_depth(
                TESS_SINGLE_DEPTH_EFFECTIVE,
                p9_apparent_magnitude_at(&p9, ALBEDO, p9.a),
            ),
        },
        "trial_tracks": n_trial_tracks(1.0),
        "published_trial_tracks": PUBLISHED_SHIFT_VECTORS,
        "rate_cell_arcsec_per_day": rate_cell_arcsec_per_day(1.0),
        "search_distance_au": [SEARCH_DISTANCE_AU.0, SEARCH_DISTANCE_AU.1],
    })
}
