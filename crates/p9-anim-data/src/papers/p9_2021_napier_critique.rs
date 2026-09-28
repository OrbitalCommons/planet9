//! Film export for `p9-2021-napier-critique`: the numbers its scene and ledger entry draw.

use p9_2021_napier_critique::critique::{
    NAPIER_2021_CONSISTENCY_BAND, SelectionFunction, directional_consistency_p_value, run_critique,
};
use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::analysis::stats::p_value_to_sigma;
use p9_core::constants::{DEG2RAD, RAD2DEG, TWO_PI};
use p9_core::data::etno::{BROWN_2017_SAMPLE, longitudes_of_perihelion};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

const SEED: u64 = 2021;
const N_MC: usize = 20_000;
/// Synthetic samples drawn per null for the R-bar histograms.
const N_NULL: usize = 20_000;
const R_BINS: usize = 25;
/// Number of ETNOs in the paper's own sample (DES + OSSOS + Sheppard & Trujillo).
const PAPER_N_ETNOS: usize = 14;

/// Histogram (density per unit R-bar) of the mean resultant length of `n`
/// angles drawn from the density proportional to `weight`.
fn null_r_bar_density(n: usize, weight: &dyn Fn(f64) -> f64, w_max: f64, seed: u64) -> Vec<f64> {
    let mut rng = StdRng::seed_from_u64(seed);
    let mut counts = vec![0usize; R_BINS];
    let mut sample = vec![0.0; n];
    for _ in 0..N_NULL {
        for s in sample.iter_mut() {
            *s = loop {
                let th = rng.gen_range(0.0..TWO_PI);
                if rng.gen_range(0.0..w_max) <= weight(th) {
                    break th;
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

pub fn export() -> Value {
    let varpis = longitudes_of_perihelion();
    let sel = SelectionFunction::default();
    let n = varpis.len();

    let mut rng = StdRng::seed_from_u64(SEED);
    let res = run_critique(&varpis, &sel, N_MC, &mut rng);

    // Direction-aware check, with the lobe where the crate assumes it and with
    // the same lobe turned a quarter of the way round the sky.
    let mut rng = StdRng::seed_from_u64(SEED);
    let p_dir_aligned = directional_consistency_p_value(&varpis, &sel, N_MC, &mut rng);
    let rotated = SelectionFunction {
        phi1: sel.phi1 + 90.0 * DEG2RAD,
        phi2: sel.phi2 + 90.0 * DEG2RAD,
        ..sel
    };
    let mut rng = StdRng::seed_from_u64(SEED);
    let p_dir_rotated = directional_consistency_p_value(&varpis, &rotated, N_MC, &mut rng);

    let lon_deg: Vec<f64> = (0..=180).map(|k| 2.0 * k as f64).collect();
    let weight: Vec<f64> = lon_deg.iter().map(|&l| sel.weight(l * DEG2RAD)).collect();
    let w_max = sel.weight_max();
    let w_min = weight.iter().cloned().fold(f64::MAX, f64::min);

    let edges: Vec<f64> = (0..=R_BINS).map(|k| k as f64 / R_BINS as f64).collect();
    let null_flat = null_r_bar_density(n, &|_| 1.0, 1.0, SEED + 1);
    let null_sel = null_r_bar_density(n, &|v| sel.weight(v), w_max, SEED + 2);

    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a": o.a,
                "e": o.e,
                "varpi_deg": o.longitude_of_perihelion() * RAD2DEG,
            })
        })
        .collect();

    json!({
        "objects": objects,
        "n_sample": n,
        "paper_n_sample": PAPER_N_ETNOS,
        "mean_varpi_deg": circular_mean(&varpis).unwrap_or(0.0).rem_euclid(TWO_PI) * RAD2DEG,
        "r_bar": res.r_bar_observed,
        "rayleigh_p": res.rayleigh_p,
        "consistency_p": res.consistency_p,
        "sigma_naive": p_value_to_sigma(res.rayleigh_p),
        "sigma": p_value_to_sigma(res.consistency_p),
        "p_directional_aligned": p_dir_aligned,
        "p_directional_rotated": p_dir_rotated,
        "paper_band": [NAPIER_2021_CONSISTENCY_BAND.0, NAPIER_2021_CONSISTENCY_BAND.1],
        "selection": {
            "a1": sel.a1,
            "phi1_deg": sel.phi1 * RAD2DEG,
            "a2": sel.a2,
            "lon_deg": lon_deg,
            "weight": weight,
            "contrast": w_max / w_min,
        },
        "null_r_bar": {
            "edges": edges,
            "flat": null_flat,
            "selection": null_sel,
        },
    })
}
