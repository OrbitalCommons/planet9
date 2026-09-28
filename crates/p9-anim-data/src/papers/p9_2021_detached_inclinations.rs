//! Film export for `p9-2021-detached-inclinations`: the numbers its scene and ledger entry draw.

use p9_2021_detached_inclinations::forcing::{
    forced_inclination_typed, giant_planet_frequency_typed, planet_nine_frequency_typed,
};
use p9_2021_detached_inclinations::nominal_inclined_p9;
use p9_2021_detached_inclinations::reference::KOZAI_LIDOV_FLOOR_DEG;
use p9_core::constants::DEG2RAD;
use p9_core::units::{julian_year, radians};
use serde_json::{Value, json};

/// Perihelion and free inclination of the orbit the forcing curve is drawn for.
const Q_CURVE_AU: f64 = 45.0;
const I_CURVE_DEG: f64 = 10.0;

pub fn export() -> Value {
    let p9 = nominal_inclined_p9();

    // The tug of war along the belt: the giant planets hold the orbit plane
    // down, Planet Nine pulls it toward its own.
    let per_yr = radians(1.0) / julian_year();
    let a_grid: Vec<f64> = (0..=55).map(|k| 150.0 + 10.0 * k as f64).collect();
    let curve: Vec<Value> = a_grid
        .iter()
        .map(|&a| {
            let e = 1.0 - Q_CURVE_AU / a;
            let i = I_CURVE_DEG * DEG2RAD;
            let b_gp = (giant_planet_frequency_typed(a, e, i) / per_yr).value;
            let b_p9 = (planet_nine_frequency_typed(a, &p9) / per_yr).value;
            let forced = (forced_inclination_typed(a, e, i, &p9) / radians(1.0)).value;
            json!({
                "a_au": a,
                "giant_period_myr": std::f64::consts::TAU / b_gp / 1.0e6,
                "p9_period_myr": std::f64::consts::TAU / b_p9 / 1.0e6,
                "forced_deg": forced / DEG2RAD,
            })
        })
        .collect();
    let a_equal = curve
        .iter()
        .find(|c| c["p9_period_myr"].as_f64() <= c["giant_period_myr"].as_f64())
        .map(|c| c["a_au"].clone());

    json!({
        "p9": {
            "mass_earth": p9.mass_earth,
            "a_au": p9.a,
            "e": p9.e,
            "i_deg": p9.i / DEG2RAD,
        },
        "curve_q_au": Q_CURVE_AU,
        "curve": curve,
        "a_equal_au": a_equal,
        "floor_deg": KOZAI_LIDOV_FLOOR_DEG,
    })
}
