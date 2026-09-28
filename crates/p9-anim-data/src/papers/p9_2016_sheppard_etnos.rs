//! Film export for `p9-2016-sheppard-etnos`: the numbers its scene and ledger entry draw.

use p9_2016_sheppard_etnos::clustering::{ClusterStats, summarize};
use p9_2016_sheppard_etnos::discoveries::{
    A_THRESHOLD_AU, ANTI_ALIGNED_NAME, OUTER_OORT_NAME, Q_THRESHOLD_AU, SHEPPARD_2016_DISCOVERIES,
    clustering_subset,
};
use p9_core::analysis::circular::wrap_to_pi;
use p9_core::data::stable_kbos::{longitude_of_perihelion, stable_kbos};
use serde_json::{Value, json};

use super::p9_2016_evidence::orbit_json;

fn stats_json(s: &ClusterStats) -> Value {
    json!({
        "n": s.n,
        "mean_deg": s.mean_deg,
        "r_bar": s.r_bar,
        "std_deg": s.std_deg,
        "rayleigh_p": s.rayleigh_p,
    })
}

pub fn export() -> Value {
    let known = stable_kbos();
    let subset = clustering_subset();

    let known_orbits: Vec<Value> = known
        .iter()
        .map(|k| orbit_json(k.name, &k.elements))
        .collect();
    let new_orbits: Vec<Value> = SHEPPARD_2016_DISCOVERIES
        .iter()
        .map(|o| {
            let mut orbit = orbit_json(o.name, &o.elements());
            orbit["in_sample"] = json!(subset.iter().any(|s| s.name == o.name));
            orbit["anti_aligned"] = json!(o.name == ANTI_ALIGNED_NAME);
            orbit["outer_oort"] = json!(o.name == OUTER_OORT_NAME);
            orbit
        })
        .collect();

    // Arguments and longitudes of perihelion: the six objects of Batygin &
    // Brown (2016), then the same six plus the two new sample members.
    let known_omega: Vec<f64> = known.iter().map(|k| k.elements.omega).collect();
    let known_varpi: Vec<f64> = known
        .iter()
        .map(|k| longitude_of_perihelion(&k.elements))
        .collect();
    let mut after_omega = known_omega.clone();
    let mut after_varpi = known_varpi.clone();
    for o in &subset {
        after_omega.push(o.omega_deg.to_radians());
        after_varpi.push(o.longitude_of_perihelion());
    }

    let before_varpi_stats = summarize(&known_varpi);
    let cluster_varpi = before_varpi_stats
        .mean_deg
        .expect("clustered sample has a mean direction")
        .to_radians();
    let offset = |name: &str| -> f64 {
        let o = SHEPPARD_2016_DISCOVERIES
            .iter()
            .find(|o| o.name == name)
            .expect("object is in the discovery table");
        wrap_to_pi(o.longitude_of_perihelion() - cluster_varpi)
            .to_degrees()
            .abs()
    };

    json!({
        "known": known_orbits,
        "new": new_orbits,
        "a_threshold": A_THRESHOLD_AU,
        "q_threshold": Q_THRESHOLD_AU,
        "n_before": known.len(),
        "n_after": known.len() + subset.len(),
        "omega_before": stats_json(&summarize(&known_omega)),
        "omega_after": stats_json(&summarize(&after_omega)),
        "varpi_before": stats_json(&before_varpi_stats),
        "varpi_after": stats_json(&summarize(&after_varpi)),
        "sr349_offset_deg": offset("2014 SR349"),
        "ft28_offset_deg": offset(ANTI_ALIGNED_NAME),
    })
}
