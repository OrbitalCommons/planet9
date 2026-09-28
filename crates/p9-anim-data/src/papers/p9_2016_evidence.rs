//! Film export for `p9-2016-evidence`: the numbers its scene and ledger entry draw.

use p9_2016_evidence::kbo_elements::{joint_clustering_significance, observed_clustering_stats};
use p9_core::analysis::circular::{mean_resultant_length, wrap_to_pi};
use p9_core::analysis::poles::{pole_resultant_length, pole_vector};
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::{GM_SUN, TWO_PI};
use p9_core::coords::sky::ecliptic_vec_to_equatorial_deg;
use p9_core::data::stable_kbos::stable_kbos;
use p9_core::types::{OrbitalElements, P9Params, elements_to_cartesian};
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

/// Points along each exported orbit.
const N_TRACK: usize = 96;
/// Null trials drawn for the scatter panel.
const N_NULL_SHOWN: usize = 1500;

fn round1(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}

/// Heliocentric ecliptic positions (AU) around one orbit, evenly spaced in
/// eccentric anomaly so both the perihelion and the aphelion arcs are smooth.
pub fn orbit_track(elements: &OrbitalElements, n: usize) -> Vec<[f64; 3]> {
    (0..=n)
        .map(|k| {
            let ecc_anomaly = TWO_PI * k as f64 / n as f64;
            let at = OrbitalElements {
                mean_anomaly: ecc_anomaly - elements.e * ecc_anomaly.sin(),
                ..*elements
            };
            let pos = elements_to_cartesian(&at, GM_SUN).pos;
            [round1(pos.x), round1(pos.y), round1(pos.z)]
        })
        .collect()
}

/// One orbit as the scenes draw it: elements in degrees and AU, the track in
/// the ecliptic frame, and where its perihelion sits on the sky.
pub fn orbit_json(name: &str, elements: &OrbitalElements) -> Value {
    let at_perihelion = OrbitalElements {
        mean_anomaly: 0.0,
        ..*elements
    };
    let peri = elements_to_cartesian(&at_perihelion, GM_SUN).pos;
    let (ra, dec) = ecliptic_vec_to_equatorial_deg(&peri);
    json!({
        "name": name,
        "a": elements.a,
        "e": elements.e,
        "q": elements.a * (1.0 - elements.e),
        "i_deg": elements.i.to_degrees(),
        "omega_deg": elements.omega.to_degrees().rem_euclid(360.0),
        "node_deg": elements.omega_big.to_degrees().rem_euclid(360.0),
        "varpi_deg": elements.longitude_of_perihelion().to_degrees().rem_euclid(360.0),
        "peri_xyz": [peri.x, peri.y, peri.z],
        "peri_ra_deg": ra,
        "peri_dec_deg": dec,
        "track": orbit_track(elements, N_TRACK),
    })
}

/// A Planet Nine parameter set as an orbit the scenes can draw.
pub fn p9_orbit_json(name: &str, p9: &P9Params) -> Value {
    let elements = OrbitalElements {
        a: p9.a,
        e: p9.e,
        i: p9.i,
        omega: p9.omega,
        omega_big: p9.omega_big,
        mean_anomaly: 0.0,
    };
    let mut orbit = orbit_json(name, &elements);
    orbit["mass_earth"] = json!(p9.mass_earth);
    orbit
}

pub fn export() -> Value {
    let joint = joint_clustering_significance(2_000_000, 2016);
    let (mean_varpi, mean_node, mean_omega) = observed_clustering_stats();
    let p9 = P9Params::nominal_2016();
    let kbos = stable_kbos();

    let objects: Vec<Value> = kbos
        .iter()
        .map(|k| orbit_json(k.name, &k.elements))
        .collect();

    // The same null the crate's Monte Carlo draws (uniform perihelion
    // longitudes and nodes, observed inclinations kept), sampled here so the
    // scene can show what chance alignments look like.
    let mut rng = rand::rngs::StdRng::seed_from_u64(2016);
    let null_trials: Vec<[f64; 2]> = (0..N_NULL_SHOWN)
        .map(|_| {
            let varpis: Vec<f64> = kbos.iter().map(|_| rng.gen_range(0.0..TWO_PI)).collect();
            let poles: Vec<[f64; 3]> = kbos
                .iter()
                .map(|k| pole_vector(k.elements.i, rng.gen_range(0.0..TWO_PI)))
                .collect();
            [
                mean_resultant_length(&varpis),
                pole_resultant_length(&poles),
            ]
        })
        .collect();

    let p9_varpi = (p9.omega + p9.omega_big).rem_euclid(TWO_PI);

    json!({
        "p_joint": joint.p_joint,
        "p_varpi": joint.p_varpi,
        "p_pole": joint.p_pole,
        "r_bar_varpi": joint.r_bar_varpi_obs,
        "r_pole": joint.r_pole_obs,
        "n_trials": joint.n_trials,
        "sigma": p_value_to_sigma(joint.p_joint),
        "n_sample": kbos.len(),
        "mean_varpi_deg": mean_varpi.to_degrees().rem_euclid(360.0),
        "mean_node_deg": mean_node.to_degrees().rem_euclid(360.0),
        "mean_omega_deg": mean_omega.to_degrees().rem_euclid(360.0),
        "objects": objects,
        "null_trials": null_trials,
        "mass": p9.mass_earth,
        "a": p9.a,
        "e": p9.e,
        "i": p9.i.to_degrees(),
        "p9": p9_orbit_json("Planet Nine", &p9),
        "p9_offset_deg": wrap_to_pi(p9_varpi - mean_varpi).to_degrees().abs(),
    })
}
