//! Film export for `p9-2021-perihelion-gap`: the numbers its scene and ledger entry draw.

use p9_2021_perihelion_gap::boundaries::scattering_boundary;
use p9_2021_perihelion_gap::distribution::{Histogram, count_in_window, dip_statistic};
use p9_2021_perihelion_gap::published::{
    GAP_ECCENTRICITY_FLOOR, GAP_Q_HIGH_AU, GAP_Q_LOW_AU, GAP_RELATIVE_ABUNDANCE,
};
use p9_2021_perihelion_gap::sample::{
    EXTENDED_DISTANT_TNOS, observed_perihelia, paper_epoch_perihelia,
};
use p9_2021_perihelion_gap::synthetic::continuous_null_dip_p_value;
use p9_core::data::etno::BROWN_2017_SAMPLE;
use serde_json::{Value, json};

/// Width of each flank window beside the gap (AU), as in the crate's tests.
const FLANK_AU: f64 = 12.0;
/// Seed of the crate's significance test.
const SEED: u64 = 2021;
/// Histogram range and bin count: 5 AU bins from 30 to 90 AU.
const Q_LO: f64 = 30.0;
const Q_HI: f64 = 90.0;
const N_BINS: usize = 12;
/// Monte Carlo draws from the single-population null.
const N_MC: usize = 20_000;
/// The one object of the table announced after the paper.
const POST_PAPER: &str = "2021 RR205";

/// Expected counts per histogram bin under the single continuous population
/// the crate's null test fits: a shifted exponential with the sample minimum
/// as origin and the mean excess as scale.
fn null_expectation(perihelia: &[f64], edges: &[f64]) -> (f64, f64, Vec<f64>) {
    let n = perihelia.len() as f64;
    let q0 = perihelia.iter().cloned().fold(f64::INFINITY, f64::min);
    let lambda = perihelia.iter().map(|q| q - q0).sum::<f64>() / n;
    let cdf = |q: f64| {
        if q <= q0 {
            0.0
        } else {
            1.0 - (-(q - q0) / lambda).exp()
        }
    };
    let expected = edges
        .windows(2)
        .map(|w| n * (cdf(w[1]) - cdf(w[0])))
        .collect();
    (q0, lambda, expected)
}

fn epoch(perihelia: &[f64], seed: u64) -> Value {
    let hist = Histogram::new(perihelia, Q_LO, Q_HI, N_BINS);
    let dip = dip_statistic(perihelia, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, FLANK_AU);
    let (q0, lambda, expected) = null_expectation(perihelia, &hist.edges);
    json!({
        "n": perihelia.len(),
        "edges": hist.edges,
        "counts": hist.counts,
        "null_expected": expected,
        "null_q0_au": q0,
        "null_scale_au": lambda,
        "n_in_gap": count_in_window(perihelia, GAP_Q_LOW_AU, GAP_Q_HIGH_AU),
        "n_low_flank": count_in_window(perihelia, GAP_Q_LOW_AU - FLANK_AU, GAP_Q_LOW_AU),
        "n_high_flank": count_in_window(perihelia, GAP_Q_HIGH_AU, GAP_Q_HIGH_AU + FLANK_AU),
        "dip_ratio": dip.dip_ratio,
        "p_single_population": continuous_null_dip_p_value(
            perihelia, GAP_Q_LOW_AU, GAP_Q_HIGH_AU, FLANK_AU, N_MC, seed,
        ),
    })
}

pub fn export() -> Value {
    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| (o.name, o.a, o.e))
        .chain(EXTENDED_DISTANT_TNOS.iter().map(|o| (o.name, o.a, o.e)))
        .map(|(name, a, e)| {
            json!({
                "name": name,
                "a_au": a,
                "e": e,
                "q_au": a * (1.0 - e),
                "post_paper": name == POST_PAPER,
            })
        })
        .collect();

    let a_grid: Vec<f64> = (0..=60).map(|k| 250.0 + 40.0 * k as f64).collect();
    let boundary: Vec<f64> = a_grid.iter().map(|&a| scattering_boundary(a)).collect();

    let paper = epoch(&paper_epoch_perihelia(), SEED);
    let today = epoch(&observed_perihelia(), SEED);

    json!({
        "gap_au": [GAP_Q_LOW_AU, GAP_Q_HIGH_AU],
        "flank_au": FLANK_AU,
        "e_floor": GAP_ECCENTRICITY_FLOOR,
        "objects": objects,
        "post_paper_name": POST_PAPER,
        "scattering_boundary": {"a_au": a_grid, "q_au": boundary},
        "dip_ratio_paper_epoch": paper["dip_ratio"],
        "p_paper_epoch": paper["p_single_population"],
        "p_today": today["p_single_population"],
        "paper_epoch": paper,
        "today": today,
        "published_gap_relative_abundance": GAP_RELATIVE_ABUNDANCE,
    })
}
