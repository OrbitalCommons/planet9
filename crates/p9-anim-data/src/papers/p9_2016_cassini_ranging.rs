//! Film export for `p9-2016-cassini-ranging`: the numbers its scene and ledger entry draw.

use p9_2016_cassini_ranging::perturbation::{
    INPOP_RESIDUAL_FLOOR_KM, favored_true_anomaly, prefit_amplitude, range_perturbation_amplitude,
};
use p9_2016_cassini_ranging::published::PREFERRED_HELIO_DISTANCE_AU;
use p9_core::data::ephemeris_constraint::{
    FAVORED_INTERVAL_DEG, PREFERRED_TRUE_ANOMALY_DEG, brown_batygin_orbit,
};
use p9_core::types::helio_distance_at_true_anomaly;
use serde_json::{Value, json};

/// True-anomaly intervals Fienga et al. (2016) exclude (degrees, abstract),
/// mapped from their (-180, 180] convention onto [0, 360).
const PUBLISHED_EXCLUDED_DEG: [(f64, f64); 3] = [(0.0, 85.0), (230.0, 260.0), (295.0, 360.0)];

/// True-anomaly samples scanned for the favoured position.
const N_SCAN: usize = 1440;

const M_PER_KM: f64 = 1.0e3;

/// Contiguous runs of `flag == true` along a closed ring of samples, as
/// (start, end) in degrees; a run through 360 is split at the seam.
fn runs(nu_deg: &[f64], flag: &[bool]) -> Vec<(f64, f64)> {
    let step = nu_deg[1] - nu_deg[0];
    let mut out = Vec::new();
    let mut start: Option<f64> = None;
    for (k, &on) in flag.iter().enumerate() {
        match (on, start) {
            (true, None) => start = Some(nu_deg[k] - 0.5 * step),
            (false, Some(s)) => {
                out.push((s.max(0.0), nu_deg[k] - 0.5 * step));
                start = None;
            }
            _ => {}
        }
    }
    if let Some(s) = start {
        out.push((s.max(0.0), 360.0));
    }
    out
}

pub fn export() -> Value {
    let orbit = brown_batygin_orbit();

    let nu_deg: Vec<f64> = (0..180).map(|k| 2.0 * k as f64).collect();
    let post_m: Vec<f64> = nu_deg
        .iter()
        .map(|nu| M_PER_KM * range_perturbation_amplitude(&orbit, nu.to_radians()))
        .collect();
    let pre_m: Vec<f64> = nu_deg
        .iter()
        .map(|nu| M_PER_KM * prefit_amplitude(&orbit, nu.to_radians()))
        .collect();
    let r_au: Vec<f64> = nu_deg
        .iter()
        .map(|nu| helio_distance_at_true_anomaly(&orbit, nu.to_radians()))
        .collect();

    let floor_m = M_PER_KM * INPOP_RESIDUAL_FLOOR_KM;
    let excluded: Vec<bool> = post_m.iter().map(|&m| m > floor_m).collect();
    let blind: Vec<bool> = pre_m.iter().map(|&m| m < floor_m).collect();

    let favored = favored_true_anomaly(&orbit, N_SCAN);

    json!({
        "orbit": {
            "mass_earth": orbit.mass_earth,
            "a_au": orbit.a,
            "e": orbit.e,
            "perihelion_au": orbit.a * (1.0 - orbit.e),
            "aphelion_au": orbit.a * (1.0 + orbit.e),
        },
        "curve": {
            "nu_deg": nu_deg,
            "postfit_m": post_m,
            "prefit_m": pre_m,
            "r_au": r_au,
        },
        "floor_m": floor_m,
        "excluded_deg": runs(&nu_deg, &excluded),
        "blind_deg": runs(&nu_deg, &blind),
        "favored_deg": favored.to_degrees(),
        "favored_distance_au": helio_distance_at_true_anomaly(&orbit, favored),
        "favored_postfit_m": M_PER_KM * range_perturbation_amplitude(&orbit, favored),
        "published": {
            "favored_deg": PREFERRED_TRUE_ANOMALY_DEG,
            "favored_interval_deg": [FAVORED_INTERVAL_DEG.0, FAVORED_INTERVAL_DEG.1],
            "favored_distance_au": PREFERRED_HELIO_DISTANCE_AU,
            "excluded_deg": PUBLISHED_EXCLUDED_DEG,
        },
    })
}
