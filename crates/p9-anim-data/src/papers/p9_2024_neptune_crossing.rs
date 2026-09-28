//! Film export for `p9-2024-neptune-crossing`: the numbers its scene and ledger entry draw.

use p9_2024_neptune_crossing::hypothesis_test::{
    HypothesisResult, compute_zeta, discovery_distances, ks_p_value, ks_test, xi_values,
    zeta_null_samples, zeta_sigma_deviation,
};
use p9_2024_neptune_crossing::observed_tnos::{observed_sample, selection_criteria};
use p9_2024_neptune_crossing::simulation::{
    SimulationOptions, SimulationResult, quick_test_simulation,
};
use p9_core::constants::YEAR_DAYS;
use serde_json::{Value, json};

/// Draws of the ζ null distribution.
const N_NULL: usize = 200_000;
/// Footprints kept per model for the (a, q) panel.
const N_SHOWN: usize = 400;

fn model(result: &SimulationResult, xi: &[f64], null: &[f64]) -> Value {
    let step = (result.n_selected / N_SHOWN).max(1);
    let footprints: Vec<Value> = result
        .selected_footprints()
        .iter()
        .step_by(step)
        .map(|f| json!({"a_au": f.a, "q_au": f.q, "i_deg": f.i_deg}))
        .collect();
    let mut sorted = xi.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let zeta = compute_zeta(xi);
    let ks = ks_test(xi);
    json!({
        "n_footprints": result.n_footprints,
        "n_crossing": result.n_selected,
        "crossing_fraction": result.selection_fraction(),
        "footprints": footprints,
        "perihelia": result.selected_perihelia.iter().step_by(step).collect::<Vec<_>>(),
        "xi_sorted": sorted,
        "n_xi_zero": xi.iter().filter(|&&x| x <= 0.0).count(),
        "zeta": zeta,
        "zeta_sigma": zeta_sigma_deviation(zeta, null),
        "ks_d": ks,
        "ks_p": ks_p_value(ks, xi.len()),
    })
}

pub fn export() -> Value {
    let sample = observed_sample();
    let r_disc = discovery_distances(&sample);
    let cuts = selection_criteria();
    let null = zeta_null_samples(sample.len(), N_NULL, 2404);
    let null_mean = null.iter().sum::<f64>() / null.len() as f64;

    let with_p9 = quick_test_simulation(true);
    let without = quick_test_simulation(false);
    let xi_p9 = xi_values(&sample, &r_disc, &with_p9.selected_footprints());
    let xi_free = xi_values(&sample, &r_disc, &without.selected_footprints());

    let objects: Vec<Value> = sample
        .iter()
        .zip(&r_disc)
        .zip(xi_p9.iter().zip(&xi_free))
        .map(|((o, &r), (&xp, &xf))| {
            json!({
                "name": o.name,
                "a_au": o.a,
                "q_au": o.q,
                "i_deg": o.i,
                "r_discovery_au": r,
                "xi_p9": xp,
                "xi_free": xf,
            })
        })
        .collect();

    // Null distribution of ζ in unit bins from −30 to 0.
    let edges: Vec<f64> = (0..=60).map(|k| -30.0 + 0.5 * k as f64).collect();
    let density: Vec<f64> = edges
        .windows(2)
        .map(|w| null.iter().filter(|&&z| z >= w[0] && z < w[1]).count() as f64 / null.len() as f64)
        .collect();

    let paper = HypothesisResult::paper_values();
    let opts = SimulationOptions::reduced_scale(true);
    let boost = opts.p9.as_ref().map(|p| p.mass_earth).unwrap_or(0.0);

    json!({
        "cuts": {"a_min_au": cuts.a_min, "q_max_au": cuts.q_max, "i_max_deg": cuts.i_max},
        "n_objects": sample.len(),
        "objects": objects,
        "run": {
            "n_particles": opts.n_particles,
            "t_myr": opts.t_days / YEAR_DAYS / 1.0e6,
            "p9_mass_earth": boost,
        },
        "null": {"edges": edges, "fraction": density, "mean": null_mean, "n_draws": N_NULL},
        "null_sigma": (null_mean - paper.zeta_null)
            / zeta_sigma_deviation(paper.zeta_null, &null),
        "with_p9": model(&with_p9, &xi_p9, &null),
        "without_p9": model(&without, &xi_free, &null),
        "paper": {
            "zeta_p9": paper.zeta_p9,
            "zeta_free": paper.zeta_null,
            "ks_p_p9": paper.p_value_p9,
            "ks_p_free": paper.p_value_null,
            "zeta_free_sigma": zeta_sigma_deviation(paper.zeta_null, &null),
            "zeta_p9_sigma": zeta_sigma_deviation(paper.zeta_p9, &null),
        },
        "paper_zeta_free_sigma": zeta_sigma_deviation(paper.zeta_null, &null),
    })
}
