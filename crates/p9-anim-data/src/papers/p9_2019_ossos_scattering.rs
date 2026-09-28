//! Film export for `p9-2019-ossos-scattering`: the numbers its scene and ledger entry draw.

use p9_2019_ossos_scattering::boundary::{DynamicalClass, boundary_curve, classify};
use p9_2019_ossos_scattering::planet_nine::{PlanetNine, SECULAR_SCULPTING_DAYS, perihelion_lift};
use p9_2019_ossos_scattering::population::{ScatteringTno, synthetic_scattering_population};
use p9_2019_ossos_scattering::published::SCATTERING_Q_BOUNDARY_AU;
use p9_core::constants::YEAR_DAYS;
use p9_core::types::P9Params;
use serde_json::{Value, json};
use std::thread;

/// Synthetic scattering objects followed (each costs a full secular level-curve search).
const N_OBJECTS: usize = 48;
/// Planet masses (Earth masses) swept at the planet's orbit.
const MASSES: [f64; 4] = [1.0, 2.5, 5.0, 10.0];

/// Perihelion lift of every object under one planet, spread over the available cores.
fn lifts(p9: &PlanetNine, population: &[ScatteringTno]) -> Vec<f64> {
    let workers = thread::available_parallelism().map_or(4, |n| n.get());
    let chunk = population.len().div_ceil(workers);
    thread::scope(|scope| {
        let handles: Vec<_> = population
            .chunks(chunk)
            .map(|part| {
                scope.spawn(move || {
                    part.iter()
                        .map(|t| perihelion_lift(p9, t.a_au, t.q_au, SECULAR_SCULPTING_DAYS))
                        .collect::<Vec<f64>>()
                })
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().expect("lift worker finished"))
            .collect()
    })
}

pub fn export() -> Value {
    let review = P9Params::revised_2019();
    let population = synthetic_scattering_population(1905, N_OBJECTS, 150.0, 900.0, 31.0, 40.0);

    let boundary: Vec<Value> = boundary_curve(100.0, 1000.0, 91)
        .iter()
        .map(|b| json!({"a": b.a_au, "q_crit": b.q_crit_au}))
        .collect();

    let mut sweep = Vec::new();
    let mut objects = Vec::new();
    let mut detached_fraction = 0.0;
    for mass in MASSES {
        let p9 = PlanetNine {
            mass_earth: mass,
            a_au: review.a,
            e: review.e,
        };
        let dq = lifts(&p9, &population);
        let detached = population
            .iter()
            .zip(&dq)
            .filter(|(t, dq)| classify(t.a_au, t.q_au + **dq) == DynamicalClass::Detached)
            .count();
        let fraction = detached as f64 / population.len() as f64;
        let mean_lift = dq.iter().sum::<f64>() / dq.len() as f64;
        sweep.push(json!({"mass": mass, "detached": fraction, "mean_lift_au": mean_lift}));
        if (mass - review.mass_earth).abs() < 1e-9 {
            detached_fraction = fraction;
            objects = population
                .iter()
                .zip(&dq)
                .map(|(t, dq)| {
                    json!({
                        "a": t.a_au,
                        "q": t.q_au,
                        "q_lifted": t.q_au + dq,
                        "detached": classify(t.a_au, t.q_au + dq) == DynamicalClass::Detached,
                    })
                })
                .collect();
        }
    }

    json!({
        "p9": {"mass": review.mass_earth, "a": review.a, "e": review.e},
        "sculpting_myr": SECULAR_SCULPTING_DAYS / YEAR_DAYS / 1e6,
        "boundary": boundary,
        "published_q_boundary": SCATTERING_Q_BOUNDARY_AU,
        "objects": objects,
        "n_objects": population.len(),
        "detached_fraction": detached_fraction,
        "sweep": sweep,
    })
}
