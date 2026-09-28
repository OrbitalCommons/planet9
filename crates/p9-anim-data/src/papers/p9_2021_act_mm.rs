//! Film export for `p9-2021-act-mm`: the numbers its scene and ledger entry draw.

use p9_2021_act_mm::detectability::{compare_reach, max_detectable_distance};
use p9_2021_act_mm::survey_model::{ActSurvey, NU_98_GHZ, NU_150_GHZ, NU_229_GHZ};
use p9_2021_act_mm::thermal_model::P9Millimeter;
use p9_core::analysis::thermal::R_EARTH_M;
use serde_json::{Value, json};
use std::f64::consts::PI;

/// Distance at which the paper tabulates the expected flux (AU).
const TABLE_DISTANCE_AU: f64 = 500.0;

/// One mass case of Naess et al. (2021): the Fortney et al. (2016) nominal
/// radius and temperature, the 150 GHz flux at 500 AU they give, and the
/// published range of the detection limit over the footprint.
struct Case {
    mass_earth: f64,
    radius_earth: f64,
    temp_k: f64,
    flux_150_mjy: f64,
    reach_au: (f64, f64),
    eliminated: f64,
}

const CASES: [Case; 2] = [
    Case {
        mass_earth: 5.0,
        radius_earth: 2.94,
        temp_k: 42.2,
        flux_150_mjy: 5.3,
        reach_au: (325.0, 625.0),
        eliminated: 0.17,
    },
    Case {
        mass_earth: 10.0,
        radius_earth: 3.46,
        temp_k: 48.3,
        flux_150_mjy: 8.5,
        reach_au: (425.0, 775.0),
        eliminated: 0.09,
    },
];

pub fn export() -> Value {
    let nominal = ActSurvey::default();
    let deepest = ActSurvey::deepest();
    let shallowest = ActSurvey::shallowest();

    let distances: Vec<f64> = (0..=90).map(|k| 200.0 + 10.0 * k as f64).collect();

    let cases: Vec<Value> = CASES
        .iter()
        .map(|c| {
            let body = P9Millimeter::new(c.mass_earth, TABLE_DISTANCE_AU);
            let flux: Vec<f64> = distances
                .iter()
                .map(|&d| P9Millimeter::new(c.mass_earth, d).flux_density_mjy(NU_150_GHZ))
                .collect();
            let wise = compare_reach(&nominal, c.mass_earth);
            json!({
                "mass_earth": c.mass_earth,
                "radius_earth": body.radius_m() / R_EARTH_M,
                "temp_k": body.effective_temp(),
                "flux_150_mjy_at_500au": body.flux_density_mjy(NU_150_GHZ),
                "flux_by_band_mjy": [
                    body.flux_density_mjy(NU_98_GHZ),
                    body.flux_density_mjy(NU_150_GHZ),
                    body.flux_density_mjy(NU_229_GHZ),
                ],
                "flux_150_mjy": flux,
                "reach_shallow_au": max_detectable_distance(&shallowest, c.mass_earth),
                "reach_nominal_au": max_detectable_distance(&nominal, c.mass_earth),
                "reach_deep_au": max_detectable_distance(&deepest, c.mass_earth),
                "wise_w1_reach_au": wise.wise_w1_au,
                "published": {
                    "radius_earth": c.radius_earth,
                    "temp_k": c.temp_k,
                    "flux_150_mjy_at_500au": c.flux_150_mjy,
                    "reach_lo_au": c.reach_au.0,
                    "reach_hi_au": c.reach_au.1,
                    "eliminated": c.eliminated,
                },
            })
        })
        .collect();

    let masses: Vec<f64> = (0..=24).map(|k| 3.0 + 0.5 * k as f64).collect();
    let reach = |s: &ActSurvey| -> Vec<f64> {
        masses
            .iter()
            .map(|&m| max_detectable_distance(s, m))
            .collect()
    };

    json!({
        "limit_lo_mjy": deepest.flux_limit_mjy,
        "limit_nominal_mjy": nominal.flux_limit_mjy,
        "limit_hi_mjy": shallowest.flux_limit_mjy,
        "area_deg2": nominal.area_deg2,
        "sky_fraction": nominal.area_sr() / (4.0 * PI),
        "table_distance_au": TABLE_DISTANCE_AU,
        "distance_au": distances,
        "cases": cases,
        "flux_5me_150_mjy": P9Millimeter::new(CASES[0].mass_earth, TABLE_DISTANCE_AU)
            .flux_density_mjy(NU_150_GHZ),
        "reach_vs_mass": {
            "mass_earth": masses,
            "deep_au": reach(&deepest),
            "nominal_au": reach(&nominal),
            "shallow_au": reach(&shallowest),
        },
    })
}
