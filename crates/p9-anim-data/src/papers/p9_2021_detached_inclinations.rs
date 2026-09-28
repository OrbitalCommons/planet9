//! Film export for `p9-2021-detached-inclinations`: the numbers its scene and ledger entry draw.

use p9_2021_detached_inclinations::forcing::{
    forced_inclination_typed, giant_planet_frequency_typed, planet_nine_frequency_typed,
};
use p9_2021_detached_inclinations::nominal_inclined_p9;
use p9_2021_detached_inclinations::population::{
    PopulationConfig, build_population, dispersion, fraction_above, mean, mean_forced_inclination,
    observed_inclinations, observed_inclinations_no_p9,
};
use p9_2021_detached_inclinations::reference::{
    DETACHED_A_MIN_AU, DETACHED_Q_MIN_AU, KOZAI_LIDOV_FLOOR_DEG,
};
use p9_core::constants::DEG2RAD;
use p9_core::units::{julian_year, radians};
use serde_json::{Value, json};

/// Perihelion and free inclination of the orbit the forcing curve is drawn for.
const Q_CURVE_AU: f64 = 45.0;
const I_CURVE_DEG: f64 = 10.0;
/// Histogram: 3° bins from 0° to 60°.
const BIN_DEG: f64 = 3.0;
const N_BINS: usize = 20;

fn histogram(inclinations_rad: &[f64]) -> Vec<f64> {
    let mut counts = vec![0.0; N_BINS];
    for &i in inclinations_rad {
        let k = (i / DEG2RAD / BIN_DEG) as usize;
        if k < N_BINS {
            counts[k] += 1.0;
        }
    }
    let n = inclinations_rad.len() as f64;
    counts.iter().map(|c| c / n).collect()
}

fn summary(inclinations_rad: &[f64]) -> Value {
    let floor = KOZAI_LIDOV_FLOOR_DEG * DEG2RAD;
    json!({
        "fraction": histogram(inclinations_rad),
        "mean_deg": mean(inclinations_rad) / DEG2RAD,
        "dispersion_deg": dispersion(inclinations_rad) / DEG2RAD,
        "fraction_above_20": fraction_above(inclinations_rad, floor),
        "fraction_below_20": 1.0 - fraction_above(inclinations_rad, floor),
    })
}

pub fn export() -> Value {
    let p9 = nominal_inclined_p9();
    let cfg = PopulationConfig::default();
    let population = build_population(&cfg);
    let with_p9 = observed_inclinations(&population, &p9);
    let without = observed_inclinations_no_p9(&population);

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

    let edges: Vec<f64> = (0..=N_BINS).map(|k| k as f64 * BIN_DEG).collect();
    let with_summary = summary(&with_p9);
    let without_summary = summary(&without);

    json!({
        "p9": {
            "mass_earth": p9.mass_earth,
            "a_au": p9.a,
            "e": p9.e,
            "i_deg": p9.i / DEG2RAD,
        },
        "detached": {"a_min_au": DETACHED_A_MIN_AU, "q_min_au": DETACHED_Q_MIN_AU},
        "population": {
            "n": cfg.n,
            "a_min_au": cfg.a_min,
            "a_max_au": cfg.a_max,
            "i_free_sigma_deg": cfg.i_free_sigma / DEG2RAD,
        },
        "curve_q_au": Q_CURVE_AU,
        "curve": curve,
        "a_equal_au": a_equal,
        "edges_deg": edges,
        "floor_deg": KOZAI_LIDOV_FLOOR_DEG,
        "mean_forced_deg": mean_forced_inclination(&population, &p9) / DEG2RAD,
        "dispersion_with_p9_deg": with_summary["dispersion_deg"],
        "dispersion_without_p9_deg": without_summary["dispersion_deg"],
        "with_p9": with_summary,
        "without_p9": without_summary,
    })
}
