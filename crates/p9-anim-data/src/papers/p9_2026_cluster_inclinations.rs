//! Film export for `p9-2026-cluster-inclinations`: the numbers its scene and ledger entry draw.

use p9_2026_cluster_inclinations::debias::{
    CDF_GRID, ConditionalCdf, P_ONE_SIGMA, P_THREE_SIGMA, default_w_grid, fit_intrinsic_width,
};
use p9_2026_cluster_inclinations::population::BirthEnvironment;
use p9_2026_cluster_inclinations::reference::{
    P_REJECT_W26, W_CLUSTER_FREE_DEG, W_CLUSTER_INFLUENCED_DEG, W_OBS_1SIGMA_DEG, W_OBS_DEG,
};
use p9_2026_cluster_inclinations::sample::{HIGH_Q_SAMPLE, discovery_positions, sample_elements};
use p9_2026_cluster_inclinations::simulation::{
    CLUSTER_FREE_P9, CLUSTER_INFLUENCED_P9, P9Config, RunConfig, run,
};
use p9_2026_cluster_inclinations::width::{
    SelectionCuts, mean_pole, relative_inclinations, rotation_to_pole,
};
use p9_core::constants::{DEG2RAD, YEAR_DAYS};
use serde_json::{Value, json};

/// Film-scale integration, far smaller than the crate's reduced-scale test
/// (about 1.5 CPU-seconds per particle per 100 Myr). The widths are measured
/// over the whole distant population, which keeps enough orbits for a fit.
const N_CLUSTER_INFLUENCED: usize = 96;
const N_CLUSTER_FREE: usize = 48;
const T_MYR: f64 = 100.0;
const ANALYSIS_FROM_MYR: f64 = 50.0;
const SNAPSHOT_MYR: f64 = 5.0;

/// The observed sample in its own mean-plane frame: (discovery latitudes,
/// inclinations), radians. Same construction as the crate's headline test.
fn observed_in_mean_plane() -> (Vec<f64>, Vec<f64>) {
    let elements = sample_elements();
    let pole = mean_pole(&elements);
    let rot = rotation_to_pole(&pole);
    let incs = relative_inclinations(&elements, &pole);
    let betas = discovery_positions()
        .iter()
        .zip(&incs)
        .map(|(p, &i)| {
            let p = rot * p;
            (p.z / p.norm()).asin().abs().min(i - 1e-6)
        })
        .collect();
    (betas, incs)
}

/// Sorted Brown (2001) probabilities P_j of the sample for width `w_deg`.
fn sorted_probabilities(w_deg: f64, betas: &[f64], incs: &[f64]) -> Vec<f64> {
    let mut p: Vec<f64> = betas
        .iter()
        .zip(incs)
        .map(|(&b, &i)| ConditionalCdf::new(w_deg * DEG2RAD, b, CDF_GRID).eval(i))
        .collect();
    p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    p
}

fn population(environment: BirthEnvironment, p9: P9Config, n_particles: usize) -> Value {
    let mut cfg = RunConfig::reduced_scale(environment, p9);
    cfg.n_particles = n_particles;
    cfg.t_days = T_MYR * 1.0e6 * YEAR_DAYS;
    cfg.analysis_from_days = ANALYSIS_FROM_MYR * 1.0e6 * YEAR_DAYS;
    cfg.snapshot_every_days = SNAPSHOT_MYR * 1.0e6 * YEAR_DAYS;
    let result = run(&cfg);
    let whole = SelectionCuts::whole_population();
    let history = result.width_history(&whole);
    let window = result.width_history(&SelectionCuts::width_sample());
    json!({
        "p9": {"mass_earth": p9.mass_earth, "a_au": p9.a, "e": p9.e, "i_deg": p9.i_deg},
        "n_particles": n_particles,
        "t_myr": history.iter().map(|h| h.0).collect::<Vec<_>>(),
        "w_deg": history.iter().map(|h| h.1).collect::<Vec<_>>(),
        "n_selected": history.iter().map(|h| h.2).collect::<Vec<_>>(),
        "w_window_deg": window.iter().map(|h| h.1).collect::<Vec<_>>(),
        "n_window": window.iter().map(|h| h.2).collect::<Vec<_>>(),
        "w_initial_deg": result.initial_width_with(&whole),
        "w_pooled_deg": result.pooled_width(&whole),
        "w_window_pooled_deg": result.final_width(),
        "kappa": result.final_kappa(),
        "survival": result.survival(),
    })
}

pub fn export() -> Value {
    // The two integrations are independent: run them side by side.
    let (cluster_influenced, cluster_free) = std::thread::scope(|s| {
        let stirred = s.spawn(|| {
            population(
                BirthEnvironment::ClusterInfluenced,
                CLUSTER_INFLUENCED_P9[0],
                N_CLUSTER_INFLUENCED,
            )
        });
        let quiet = population(
            BirthEnvironment::ClusterFree,
            CLUSTER_FREE_P9[6],
            N_CLUSTER_FREE,
        );
        (stirred.join().unwrap(), quiet)
    });

    let (betas, incs) = observed_in_mean_plane();
    let fit = fit_intrinsic_width(&betas, &incs, &default_w_grid(), 2000, 1);
    let (lo, hi) = fit.interval(P_ONE_SIGMA);
    let rejection = fit_intrinsic_width(&betas, &incs, &[26.0], 20_000, 2);

    let objects: Vec<Value> = HIGH_Q_SAMPLE
        .iter()
        .zip(betas.iter().zip(&incs))
        .map(|(o, (&b, &i))| {
            json!({
                "name": o.name,
                "a_au": o.a,
                "q_au": o.perihelion(),
                "i_ecliptic_deg": o.i_deg,
                "i_deg": i / DEG2RAD,
                "beta_deg": b / DEG2RAD,
            })
        })
        .collect();

    json!({
        "n_observed": HIGH_Q_SAMPLE.len(),
        "objects": objects,
        "scan": {
            "w_deg": fit.w_grid_deg,
            "statistic": fit.statistic,
            "p_value": fit.p_value,
        },
        "w_best_deg": fit.w_best_deg,
        "w_lo_deg": lo,
        "w_hi_deg": hi,
        "p_one_sigma": P_ONE_SIGMA,
        "p_three_sigma": P_THREE_SIGMA,
        "p_at_26": rejection.p_value[0],
        "probabilities_best": sorted_probabilities(fit.w_best_deg, &betas, &incs),
        "probabilities_26": sorted_probabilities(26.0, &betas, &incs),
        "published": {
            "w_obs_deg": W_OBS_DEG,
            "w_obs_1sigma_deg": [W_OBS_1SIGMA_DEG.0, W_OBS_1SIGMA_DEG.1],
            "w_cluster_influenced_deg": [W_CLUSTER_INFLUENCED_DEG.0, W_CLUSTER_INFLUENCED_DEG.1],
            "w_cluster_free_deg": [W_CLUSTER_FREE_DEG.0, W_CLUSTER_FREE_DEG.1],
            "p_reject_w26": P_REJECT_W26,
        },
        "t_myr": T_MYR,
        "analysis_from_myr": ANALYSIS_FROM_MYR,
        "cluster_influenced": cluster_influenced,
        "cluster_free": cluster_free,
    })
}
