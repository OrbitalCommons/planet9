//! Film export for `p9-2021-orbit`: the numbers its scene and ledger entry draw.

use p9_2021_orbit::statistical_measures::{
    NullModel, clustering_significance_mc, paper_clustering_confidence, survey_bias_weight,
};
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::{DEG2RAD, GM_SUN, RAD2DEG, TWO_PI};
use p9_core::data::etno::{BROWN_2017_SAMPLE, longitudes_of_perihelion};
use p9_core::data::posterior::{
    A_Q_CORRELATION, AsymmetricGaussian, mcmc_2021_posterior, sample_from_posterior,
};
use p9_core::types::{OrbitalElements, P9Params, elements_to_cartesian, true_to_mean_anomaly};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

const SEED: u64 = 2021;
const N_MC: usize = 200_000;
const N_NULL: usize = 20_000;
const R_BINS: usize = 25;
/// Posterior draws shown in the scene.
const N_POSTERIOR: usize = 1500;
/// Points along each drawn orbit.
const N_PATH: usize = 120;
/// Objects in the paper's own clustering sample.
const PAPER_N_SAMPLE: usize = 11;
/// Published clustering confidence.
const PAPER_CONFIDENCE: f64 = 0.996;

/// The orbit projected on the ecliptic plane, as (x, y) in AU.
fn orbit_path(elements: &OrbitalElements) -> Vec<(f64, f64)> {
    (0..=N_PATH)
        .map(|k| {
            let nu = TWO_PI * k as f64 / N_PATH as f64;
            let at = OrbitalElements {
                mean_anomaly: true_to_mean_anomaly(elements.e, nu),
                ..*elements
            };
            let pos = elements_to_cartesian(&at, GM_SUN).pos;
            (pos.x, pos.y)
        })
        .collect()
}

/// Density (per unit R-bar) of the mean resultant length of `n` longitudes
/// drawn from the crate's survey-bias weight.
fn null_r_bar_density(n: usize, biased: bool, seed: u64) -> Vec<f64> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut counts = vec![0usize; R_BINS];
    let mut sample = vec![0.0; n];
    for _ in 0..N_NULL {
        for s in sample.iter_mut() {
            *s = loop {
                let lon = rng.gen_range(0.0..TWO_PI);
                if !biased || rng.gen_range(0.0..1.0) < survey_bias_weight(lon) {
                    break lon;
                }
            };
        }
        let bin = (mean_resultant_length(&sample) * R_BINS as f64) as usize;
        counts[bin.min(R_BINS - 1)] += 1;
    }
    counts
        .iter()
        .map(|&c| c as f64 * R_BINS as f64 / N_NULL as f64)
        .collect()
}

fn interval(g: &AsymmetricGaussian) -> Value {
    json!({"median": g.median, "plus": g.sigma_upper, "minus": g.sigma_lower})
}

pub fn export() -> Value {
    let varpis = longitudes_of_perihelion();
    let n = varpis.len();
    let rayleigh = paper_clustering_confidence();
    let mc_bias = clustering_significance_mc(&varpis, N_MC, SEED, NullModel::SurveyBias);
    let mc_flat = clustering_significance_mc(&varpis, N_MC, SEED, NullModel::Uniform);

    let lon_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
    let bias_weight: Vec<f64> = lon_deg
        .iter()
        .map(|&l| survey_bias_weight(l * DEG2RAD))
        .collect();
    let edges: Vec<f64> = (0..=R_BINS).map(|k| k as f64 / R_BINS as f64).collect();

    // The published posterior, resampled by the workspace's emulator.
    let post = mcmc_2021_posterior();
    let mut rng = StdRng::seed_from_u64(SEED);
    let samples: Vec<Value> = (0..N_POSTERIOR)
        .map(|_| {
            let p = sample_from_posterior(&post, &mut rng);
            json!({
                "mass": p.mass_earth,
                "a": p.a,
                "e": p.e,
                "q": p.a * (1.0 - p.e),
                "i_deg": p.i * RAD2DEG,
            })
        })
        .collect();

    let best = P9Params::mcmc_2021();
    let best_elements = OrbitalElements {
        a: best.a,
        e: best.e,
        i: best.i,
        omega: best.omega,
        omega_big: best.omega_big,
        mean_anomaly: 0.0,
    };

    let etnos: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "q": o.perihelion(),
                "varpi_deg": o.longitude_of_perihelion() * RAD2DEG,
                "path": orbit_path(&o.elements()),
            })
        })
        .collect();

    json!({
        "n_sample": n,
        "paper_n_sample": PAPER_N_SAMPLE,
        "r_bar": rayleigh.mean_resultant_length,
        "mean_varpi_deg": rayleigh.mean_direction * RAD2DEG,
        "confidence": rayleigh.confidence,
        "confidence_mc_bias": mc_bias.confidence,
        "confidence_mc_flat": mc_flat.confidence,
        "paper_confidence": PAPER_CONFIDENCE,
        "sigma": p_value_to_sigma(rayleigh.p_value),
        "sigma_mc_bias": p_value_to_sigma(mc_bias.p_value),
        "paper_sigma": p_value_to_sigma(1.0 - PAPER_CONFIDENCE),
        "bias": {"lon_deg": lon_deg, "weight": bias_weight},
        "null_r_bar": {
            "edges": edges,
            "flat": null_r_bar_density(n, false, SEED + 1),
            "biased": null_r_bar_density(n, true, SEED + 2),
        },
        "posterior": {
            "mass": interval(&post.mass),
            "a": interval(&post.a),
            "i": interval(&post.i),
            "q": interval(&post.perihelion),
            "a_q_correlation": A_Q_CORRELATION,
            "samples": samples,
        },
        "mass": best.mass_earth,
        "a": best.a,
        "e": best.e,
        "i": best.i * RAD2DEG,
        "q": best.a * (1.0 - best.e),
        "varpi_deg": ((best.omega + best.omega_big) * RAD2DEG).rem_euclid(360.0),
        "p9_path": orbit_path(&best_elements),
        "etnos": etnos,
    })
}
