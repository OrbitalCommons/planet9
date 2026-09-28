//! Film export for `p9-2019-clustering`: the numbers its scene and ledger entry draw.

use p9_2019_clustering::clustering_analysis::{
    NullModel, compute_poincare_states, monte_carlo_clustering, perihelion_direction,
    survey_bias_weight,
};
use p9_2019_clustering::kbo_sample::{DistantKbo, paper_sample_a230};
use p9_2019_clustering::ossos_comparison::{
    detectable_resultant_threshold, ossos_sample, sensitivity_analysis,
};
use p9_2019_clustering::poincare_variables::{
    PoincareState, mean_state, perihelion_clustering, pole_clustering,
};
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::TWO_PI;
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::types::OrbitalElements;
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

use super::p9_2016_evidence::orbit_json;

/// Monte Carlo iterations for the significance test.
const N_TEST: usize = 300_000;
/// Synthetic samples drawn for the panels.
const N_SHOWN: usize = 1200;
/// Sample sizes over which the detection threshold is traced.
const THRESHOLD_N: (usize, usize) = (3, 30);

pub fn export() -> Value {
    let kbos = paper_sample_a230();
    let biased = monte_carlo_clustering(&kbos, N_TEST, 2019, NullModel::SurveyBias);
    let uniform = monte_carlo_clustering(&kbos, N_TEST, 2019, NullModel::Uniform);

    let states = compute_poincare_states(&kbos);
    let mean = mean_state(&states);
    let ossos = ossos_sample();
    let objects: Vec<Value> = kbos
        .iter()
        .zip(&states)
        .map(|(k, s)| {
            let mut orbit = orbit_json(k.name, &k.elements);
            orbit["x"] = json!(s.x);
            orbit["y"] = json!(s.y);
            orbit["p"] = json!(s.p);
            orbit["q_var"] = json!(s.q_var);
            orbit["ossos"] = json!(ossos.iter().any(|o| o.name == k.name));
            orbit["in_2017"] = json!(BROWN_2017_SAMPLE.iter().any(|o| o.name == k.name));
            orbit
        })
        .collect();

    // Mean positions of synthetic fourteen-object samples drawn from the
    // survey-bias null: each object keeps (a, e, i), its angles are drawn
    // where the surveys could have found it.
    let mut rng = rand::rngs::StdRng::seed_from_u64(2019);
    let null_means: Vec<Value> = (0..N_SHOWN)
        .map(|_| {
            let drawn: Vec<PoincareState> = kbos
                .iter()
                .map(|k| {
                    loop {
                        let omega = rng.gen_range(0.0..TWO_PI);
                        let omega_big = rng.gen_range(0.0..TWO_PI);
                        let (lambda, beta) = perihelion_direction(k.elements.i, omega, omega_big);
                        if rng.gen_range(0.0..1.0) < survey_bias_weight(lambda, beta) {
                            return PoincareState::from_elements(&OrbitalElements {
                                omega,
                                omega_big,
                                ..k.elements
                            });
                        }
                    }
                })
                .collect();
            let m = mean_state(&drawn);
            json!({
                "x": m.x,
                "y": m.y,
                "p": m.p,
                "q_var": m.q_var,
                "perihelion": perihelion_clustering(&m),
                "pole": pole_clustering(&m),
            })
        })
        .collect();

    let sensitivity = sensitivity_analysis(&kbos, &ossos);

    // How aligned a sample of n perihelion directions must be before the
    // Rayleigh test flags it at 95%, against the alignment the full sample
    // and the OSSOS subsample actually show.
    let varpis = |set: &[DistantKbo]| -> Vec<f64> {
        set.iter()
            .map(|k| k.elements.longitude_of_perihelion())
            .collect()
    };
    let threshold_n: Vec<usize> = (THRESHOLD_N.0..=THRESHOLD_N.1).collect();
    let threshold_r: Vec<f64> = threshold_n
        .iter()
        .map(|&n| detectable_resultant_threshold(n, 0.05))
        .collect();

    json!({
        "objects": objects,
        "n_sample": kbos.len(),
        "mean": {"x": mean.x, "y": mean.y, "p": mean.p, "q_var": mean.q_var},
        "observed_perihelion": biased.observed_perihelion,
        "observed_pole": biased.observed_pole,
        "mean_varpi_deg": biased.mean_varpi.to_degrees().rem_euclid(360.0),
        "mean_node_deg": biased.mean_omega.to_degrees().rem_euclid(360.0),
        "p_perihelion": biased.p_perihelion,
        "p_pole": biased.p_pole,
        "p_combined": biased.p_combined,
        "p_combined_uniform": uniform.p_combined,
        "p_perihelion_uniform": uniform.p_perihelion,
        "p_pole_uniform": uniform.p_pole,
        "n_iterations": biased.n_iterations,
        "sigma": p_value_to_sigma(biased.p_combined),
        "null_means": null_means,
        "ossos": {
            "n": sensitivity.n_ossos,
            "perihelion": sensitivity.ossos_perihelion,
            "pole": sensitivity.ossos_pole,
            "r_bar_needed": detectable_resultant_threshold(sensitivity.n_ossos, 0.05),
            "r_bar_varpi": mean_resultant_length(&varpis(&ossos)),
        },
        "r_bar_varpi": mean_resultant_length(&varpis(&kbos)),
        "r_bar_needed": detectable_resultant_threshold(kbos.len(), 0.05),
        "threshold": {"n": threshold_n, "r_bar": threshold_r},
    })
}
