//! Film export for `p9-2017-bias`: the numbers its scene and ledger entry draw.

use p9_2017_bias::bias_function::{BiasParams, angular_bias};
use p9_2017_bias::clustering_test::monte_carlo_clustering_test;
use p9_2017_bias::kbo_sample::{longitudes_of_perihelion, paper_sample_a230};
use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::TWO_PI;
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

use super::p9_2016_evidence::orbit_json;

/// Monte Carlo iterations for the significance test.
const N_TEST: usize = 400_000;
/// Synthetic samples drawn for the histogram of chance alignments.
const N_HIST: usize = 40_000;
/// Histogram bins in mean resultant length.
const N_BINS: usize = 40;
/// Longitude steps of the exported bias curve.
const N_LON: usize = 72;

fn histogram(values: &[f64]) -> Vec<f64> {
    let mut counts = vec![0.0; N_BINS];
    for &v in values {
        let k = ((v * N_BINS as f64) as usize).min(N_BINS - 1);
        counts[k] += 1.0 / values.len() as f64;
    }
    counts
}

pub fn export() -> Value {
    let kbos = paper_sample_a230();
    let params = BiasParams::default();
    let test = monte_carlo_clustering_test(&kbos, N_TEST, 2017);

    let varpis = longitudes_of_perihelion(&kbos);
    let r_bar = mean_resultant_length(&varpis);
    let mean_varpi = circular_mean(&varpis).expect("clustered sample has a mean direction");

    // Relative discovery probability against longitude of perihelion: the
    // crate's angular bias averaged over argument of perihelion and over the
    // inclinations of the ten objects.
    let lon_deg: Vec<f64> = (0..=N_LON)
        .map(|k| 360.0 * k as f64 / N_LON as f64)
        .collect();
    let bias: Vec<f64> = lon_deg
        .iter()
        .map(|&lon| {
            let mut sum = 0.0;
            for kbo in &kbos {
                for j in 0..N_LON {
                    let omega = TWO_PI * (j as f64 + 0.5) / N_LON as f64;
                    sum += angular_bias(lon.to_radians(), omega, kbo.i_deg.to_radians(), &params);
                }
            }
            sum / (kbos.len() * N_LON) as f64
        })
        .collect();
    let bias_peak = bias.iter().cloned().fold(0.0, f64::max);

    // Chance alignments of ten objects: drawn uniformly, and drawn from the
    // bias (rejection sampling against the crate's angular bias, as its own
    // Monte Carlo does).
    let mut rng = rand::rngs::StdRng::seed_from_u64(2017);
    let mut uniform = Vec::with_capacity(N_HIST);
    let mut biased = Vec::with_capacity(N_HIST);
    for _ in 0..N_HIST {
        let flat: Vec<f64> = kbos.iter().map(|_| rng.gen_range(0.0..TWO_PI)).collect();
        uniform.push(mean_resultant_length(&flat));
        let drawn: Vec<f64> = kbos
            .iter()
            .map(|kbo| {
                loop {
                    let varpi = rng.gen_range(0.0..TWO_PI);
                    let omega = rng.gen_range(0.0..TWO_PI);
                    let weight = angular_bias(varpi, omega, kbo.i_deg.to_radians(), &params);
                    if rng.gen_range(0.0..1.0) < weight {
                        return varpi;
                    }
                }
            })
            .collect();
        biased.push(mean_resultant_length(&drawn));
    }
    let exceed = |set: &[f64]| set.iter().filter(|&&r| r >= r_bar).count() as f64 / N_HIST as f64;

    let objects: Vec<Value> = kbos
        .iter()
        .map(|k| orbit_json(k.name, &k.elements()))
        .collect();
    let edges: Vec<f64> = (0..=N_BINS).map(|k| k as f64 / N_BINS as f64).collect();

    json!({
        "objects": objects,
        "n_sample": kbos.len(),
        "r_bar": r_bar,
        "mean_varpi_deg": mean_varpi.to_degrees().rem_euclid(360.0),
        "p_varpi": test.p_varpi,
        "p_omega": test.p_omega,
        "p_combined": test.p_combined,
        "n_iterations": test.n_iterations,
        "sigma_varpi": p_value_to_sigma(test.p_varpi),
        "sigma": p_value_to_sigma(test.p_combined),
        "bias": {
            "lon_deg": lon_deg,
            "weight": bias.iter().map(|b| b / bias_peak).collect::<Vec<_>>(),
        },
        "chance": {
            "edges": edges,
            "uniform": histogram(&uniform),
            "biased": histogram(&biased),
            "p_uniform": exceed(&uniform),
            "p_biased": exceed(&biased),
        },
    })
}
