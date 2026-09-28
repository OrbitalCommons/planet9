//! Film export for `p9-2018-bp519`: the numbers its scene and ledger entry draw.

use p9_2018_bp519::bp519::{discovery_2018, jpl_current};
use p9_2018_bp519::clustering::{
    concentration_with_bp519, sample_mean_varpi, sample_varpi_concentration,
    varpi_offset_from_cluster_deg,
};
use p9_2018_bp519::pumping::{
    PumpConfig, PumpingHistory, integrate, no_planet_nine, scattered_particle,
};
use p9_core::constants::YEAR_DAYS;
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::types::P9Params;
use serde_json::{Value, json};
use std::thread;

use super::p9_2016_evidence::orbit_json;

/// Starting orbit of the pumped test particle (the crate's calibrated case).
const PARTICLE_A: f64 = 250.0;
const PARTICLE_E: f64 = 0.85;
/// Starting inclinations to Planet Nine's plane (deg).
const START_INCLINATIONS: [f64; 2] = [25.0, 40.0];
/// Integration steps between exported samples.
const RECORD_EVERY: usize = 400;

fn history_json(start_deg: f64, h: &PumpingHistory) -> Value {
    json!({
        "start_deg": start_deg,
        "t_myr": h.times.iter().map(|t| t / YEAR_DAYS / 1e6).collect::<Vec<_>>(),
        "i_deg": h.inclination_deg,
        "e": h.eccentricity,
        "max_i_deg": h.inclination_deg.iter().cloned().fold(f64::NEG_INFINITY, f64::max),
    })
}

pub fn export() -> Value {
    let bp = discovery_2018();
    let jpl = jpl_current();
    let p9 = P9Params::nominal_2016();

    let sample: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| orbit_json(o.name, &o.elements()))
        .collect();

    // Two starts under Planet Nine and one control without it, each its own
    // secular integration, run side by side.
    let with_p9 = PumpConfig {
        record_every: RECORD_EVERY,
        ..PumpConfig::fast()
    };
    let without_p9 = PumpConfig {
        include_j2: true,
        ..with_p9
    };
    let control_start = START_INCLINATIONS[1];
    let runs = [
        (START_INCLINATIONS[0], p9, with_p9),
        (START_INCLINATIONS[1], p9, with_p9),
        (control_start, no_planet_nine(), without_p9),
    ];
    let mut histories: Vec<Value> = thread::scope(|scope| {
        let handles: Vec<_> = runs
            .iter()
            .map(|&(i0, planet, cfg)| {
                scope.spawn(move || {
                    let particle = scattered_particle(PARTICLE_A, PARTICLE_E, i0);
                    history_json(i0, &integrate(particle, &planet, cfg))
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("secular integration finished"))
            .collect()
    });
    let control = histories.pop().expect("control run");
    let pumped = histories;
    let max_pumped = pumped
        .iter()
        .map(|h| h["max_i_deg"].as_f64().expect("finite inclination"))
        .fold(f64::INFINITY, f64::min);

    json!({
        "bp519": orbit_json("2015 BP519", &bp.elements()),
        "h_mag": bp.h_mag,
        "i_deg": bp.i_deg,
        "jpl": {"a": jpl.a, "e": jpl.e, "i_deg": jpl.i_deg},
        "sample": sample,
        "max_sample_i_deg": BROWN_2017_SAMPLE.iter().map(|o| o.i_deg).fold(0.0, f64::max),
        "n_sample": BROWN_2017_SAMPLE.len() + 1,
        "particle": {"a": PARTICLE_A, "e": PARTICLE_E},
        "p9": {"mass": p9.mass_earth, "a": p9.a, "e": p9.e, "i_deg": p9.i.to_degrees()},
        "pumped": pumped,
        "control": control,
        "pumped_i_deg": max_pumped,
        "cluster_varpi_deg": sample_mean_varpi().to_degrees().rem_euclid(360.0),
        "varpi_offset_deg": varpi_offset_from_cluster_deg(&bp),
        "r_bar_before": sample_varpi_concentration(),
        "r_bar_after": concentration_with_bp519(&bp),
    })
}
