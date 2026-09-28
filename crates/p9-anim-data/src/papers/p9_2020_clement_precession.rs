//! Film export for `p9-2020-clement-precession`: the numbers its scene and ledger entry draw.

use p9_2020_clement_precession::precession::{
    giant_planet_precession_rate, p9_precession_rate, precession_period_gyr,
};
use p9_core::analysis::circular::mean_resultant_length;
use p9_core::constants::{GYR_DAYS, TWO_PI};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::types::P9Params;
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

/// Sample size of the paper's observed ETNO + IOCO set.
const N_PAPER: usize = 17;
/// Random samples drawn for each null distribution.
const N_DRAWS: usize = 20_000;
/// Perihelion of the reference period curve (AU).
const Q_CURVE: f64 = 40.0;

/// Mean resultant length of `n` uniformly random apsidal angles, `N_DRAWS` times.
fn uniform_r_bar(n: usize, seed: u64) -> Vec<f64> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut out: Vec<f64> = (0..N_DRAWS)
        .map(|_| {
            let angles: Vec<f64> = (0..n).map(|_| rng.gen_range(0.0..TWO_PI)).collect();
            mean_resultant_length(&angles)
        })
        .collect();
    out.sort_by(|a, b| a.partial_cmp(b).unwrap());
    out
}

/// A few individual draws of `n` uniformly random apsidal angles (degrees),
/// each with its mean resultant length, to show what one null sample looks like.
fn example_draws(n: usize, count: usize, seed: u64) -> Vec<Value> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    (0..count)
        .map(|_| {
            let angles: Vec<f64> = (0..n).map(|_| rng.gen_range(0.0..TWO_PI)).collect();
            json!({
                "varpi_deg": angles.iter().map(|a| a.to_degrees()).collect::<Vec<_>>(),
                "r_bar": mean_resultant_length(&angles),
            })
        })
        .collect()
}

fn quantile(sorted: &[f64], f: f64) -> f64 {
    sorted[((sorted.len() - 1) as f64 * f).round() as usize]
}

fn histogram(sorted: &[f64], edges: &[f64]) -> Vec<f64> {
    edges
        .windows(2)
        .map(|w| {
            sorted.iter().filter(|&&r| r >= w[0] && r < w[1]).count() as f64 / sorted.len() as f64
        })
        .collect()
}

pub fn export() -> Value {
    let p9 = P9Params::revised_2019();

    // Each observed object's apsidal precession under the giant planets alone.
    let varpi0: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| o.longitude_of_perihelion())
        .collect();
    let rates: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| giant_planet_precession_rate(o.a, o.e))
        .collect();
    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .zip(&rates)
        .map(|(o, &rate)| {
            json!({
                "name": o.name,
                "a_au": o.a,
                "e": o.e,
                "q_au": o.perihelion(),
                "varpi_deg": o.longitude_of_perihelion().to_degrees(),
                "giant_period_gyr": precession_period_gyr(rate),
                "p9_period_gyr": precession_period_gyr(p9_precession_rate(o.a, o.e, &p9)),
            })
        })
        .collect();
    let periods: Vec<f64> = rates.iter().map(|&r| precession_period_gyr(r)).collect();
    let fastest = periods.iter().cloned().fold(f64::INFINITY, f64::min);
    let slowest = periods.iter().cloned().fold(0.0, f64::max);

    // Let the observed apsides drift at those rates: how long does the
    // alignment survive with nothing holding it together?
    let t_gyr: Vec<f64> = (0..=200).map(|k| 0.02 * k as f64).collect();
    let r_bar_t: Vec<f64> = t_gyr
        .iter()
        .map(|&t| {
            let now: Vec<f64> = varpi0
                .iter()
                .zip(&rates)
                .map(|(&w, &rate)| w + rate * t * GYR_DAYS)
                .collect();
            mean_resultant_length(&now)
        })
        .collect();
    let r_bar_observed = mean_resultant_length(&varpi0);

    // Period against semi-major axis at fixed perihelion.
    let a_grid: Vec<f64> = (0..=65).map(|k| 150.0 + 10.0 * k as f64).collect();
    let giant_curve: Vec<f64> = a_grid
        .iter()
        .map(|&a| precession_period_gyr(giant_planet_precession_rate(a, 1.0 - Q_CURVE / a)))
        .collect();
    let p9_curve: Vec<f64> = a_grid
        .iter()
        .map(|&a| precession_period_gyr(p9_precession_rate(a, 1.0 - Q_CURVE / a, &p9)))
        .collect();

    // The paper's caution: a small sample of uniformly random orbits still
    // shows some clustering.
    let null_paper = uniform_r_bar(N_PAPER, 2005);
    let null_sample = uniform_r_bar(BROWN_2017_SAMPLE.len(), 5326);
    let edges: Vec<f64> = (0..=25).map(|k| 0.04 * k as f64).collect();
    let exceed = null_sample.iter().filter(|&&r| r >= r_bar_observed).count();

    json!({
        "p9": {"mass_earth": p9.mass_earth, "a_au": p9.a, "e": p9.e},
        "objects": objects,
        "n_objects": BROWN_2017_SAMPLE.len(),
        "giant_period_fastest_gyr": fastest,
        "giant_period_slowest_gyr": slowest,
        "drift": {"t_gyr": t_gyr, "r_bar": r_bar_t},
        "r_bar_observed": r_bar_observed,
        "r_bar_after_4gyr": r_bar_t.last().copied(),
        "curve": {
            "q_au": Q_CURVE,
            "a_au": a_grid,
            "giant_period_gyr": giant_curve,
            "p9_period_gyr": p9_curve,
        },
        "null": {
            "n_paper": N_PAPER,
            "n_draws": N_DRAWS,
            "edges": edges,
            "fraction_paper": histogram(&null_paper, &edges),
            "fraction_sample": histogram(&null_sample, &edges),
        },
        "examples_17": example_draws(N_PAPER, 4, 17),
        "r_bar_null_mean_17":null_paper.iter().sum::<f64>() / null_paper.len() as f64,
        "r_bar_null_p95_17": quantile(&null_paper, 0.95),
        "r_bar_null_mean_sample": null_sample.iter().sum::<f64>() / null_sample.len() as f64,
        "r_bar_null_p95_sample": quantile(&null_sample, 0.95),
        "p_observed_by_chance": exceed as f64 / null_sample.len() as f64,
    })
}
