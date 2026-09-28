//! Film export for `p9-2020-resonance-hopping`: the numbers its scene and ledger entry draw.

use p9_2020_resonance_hopping::classification::{
    Class, ResonanceLandscape, classify, resonance_half_width, resonance_overlap_k,
};
use p9_2020_resonance_hopping::population::{
    class_fractions, default_landscape, default_m_p9_solar, synthetic_population,
};
use p9_2020_resonance_hopping::published::{A_P9_AU, E_P9};
use p9_core::analysis::resonance::{critical_perihelion, neptune_diffusion_coefficient};
use p9_core::constants::EARTH_MASS_SOLAR;
use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

const A_MIN: f64 = 150.0;
const A_MAX: f64 = 499.0;
/// The crate's labelled O(1) pendulum-width constant (its test value).
const WIDTH_COEFF: f64 = 1.2;
/// Perihelion range of the synthetic population (AU), the crate's test range.
const Q_MIN: f64 = 30.0;
const Q_MAX: f64 = 120.0;
/// Perihelia at which Neptune's kicks are evaluated (AU).
const Q_KICKS: [f64; 3] = [30.0, 33.0, 36.0];
const N_BINS: usize = 7;
/// Random walk under Neptune's kicks: perihelion (AU), step (yr), steps.
const WALK_Q: f64 = 33.0;
const WALK_DT_YR: f64 = 5_000.0;
const WALK_STEPS: usize = 600;

fn class_name(c: Class) -> &'static str {
    match c {
        Class::Resonant => "resonant",
        Class::Hopping => "hopping",
        Class::NonResonant => "non_resonant",
    }
}

pub fn export() -> Value {
    let landscape: ResonanceLandscape = default_landscape(A_MIN, A_MAX);
    let m9 = default_m_p9_solar();
    let population = synthetic_population(2020, 6000, A_MIN, A_MAX, Q_MIN, Q_MAX);
    let all = class_fractions(&population, &landscape, m9, WIDTH_COEFF);

    // Class fractions in equal-width semi-major-axis bins.
    let width = (A_MAX - A_MIN) / N_BINS as f64;
    let bins: Vec<Value> = (0..N_BINS)
        .map(|k| {
            let (lo, hi) = (A_MIN + k as f64 * width, A_MIN + (k + 1) as f64 * width);
            let members: Vec<_> = population
                .iter()
                .copied()
                .filter(|t| t.a >= lo && t.a < hi)
                .collect();
            let f = class_fractions(&members, &landscape, m9, WIDTH_COEFF);
            json!({
                "a_lo": lo,
                "a_hi": hi,
                "n": f.n,
                "resonant": f.resonant,
                "hopping": f.hopping,
                "non_resonant": f.non_resonant,
            })
        })
        .collect();

    // A thinned view of the classified population for the (a, q) panel.
    let sample: Vec<Value> = population
        .iter()
        .take(700)
        .map(|t| {
            let c = classify(&landscape, t.a, t.e, m9, WIDTH_COEFF);
            json!({"a_au": t.a, "q_au": t.a * (1.0 - t.e), "class": class_name(c.class)})
        })
        .collect();

    // The low-order commensurabilities with their libration widths at the
    // eccentricity of a q = 40 AU orbit.
    let simple: Vec<Value> = landscape
        .simple
        .iter()
        .map(|&(mmr, a_res)| {
            let e = 1.0 - 40.0 / a_res;
            json!({
                "label": format!("{}:{}", mmr.p, mmr.q),
                "a_au": a_res,
                "half_width_au": resonance_half_width(a_res, mmr, e, m9, WIDTH_COEFF),
            })
        })
        .collect();

    // Overlap of neighbouring Planet Nine resonances along the belt.
    let a_grid: Vec<f64> = (0..=349).map(|k| A_MIN + k as f64).collect();
    let overlap: Vec<Option<f64>> = a_grid
        .iter()
        .map(|&a| resonance_overlap_k(&landscape, a, 1.0 - 40.0 / a, m9, WIDTH_COEFF))
        .collect();

    let median_spacing = {
        let mut gaps: Vec<f64> = landscape
            .simple
            .windows(2)
            .map(|w| w[1].1 - w[0].1)
            .collect();
        gaps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        gaps[gaps.len() / 2]
    };

    // Neptune's kicks: the time for the semi-major-axis random walk to cover
    // the spacing of the low-order Planet Nine resonances, from the Batygin,
    // Mardling & Nesvorný (2021) diffusion coefficient. That coefficient holds
    // inside Neptune's chaotic layer, q < q_crit(a), so the smallest
    // semi-major axis at which it applies is exported with it.
    let kicks: Vec<Value> = Q_KICKS
        .iter()
        .map(|&q| {
            let d = neptune_diffusion_coefficient(q);
            let a_valid = (150..=1500)
                .map(|a| a as f64)
                .find(|&a| critical_perihelion(a) >= q);
            json!({
                "q_au": q,
                "d_au2_per_yr": d,
                "crossing_time_myr": median_spacing * median_spacing / d / 1.0e6,
                "applies_beyond_a_au": a_valid,
            })
        })
        .collect();

    // One random walk in semi-major axis under Neptune's kicks at q = 33 AU,
    // started on the 3:2 resonance with Planet Nine: steps of √(2 D dt).
    let d_walk = neptune_diffusion_coefficient(WALK_Q);
    let a_start = landscape.simple.iter().map(|&(_, a)| a).fold(0.0, f64::max);
    let mut rng = rand::rngs::StdRng::seed_from_u64(2234);
    let step = (2.0 * d_walk * WALK_DT_YR).sqrt();
    let mut a_walk = vec![a_start];
    for _ in 0..WALK_STEPS {
        let last = *a_walk.last().unwrap();
        let kick: f64 = if rng.gen_bool(0.5) { step } else { -step };
        a_walk.push((last + kick).clamp(A_MIN, A_MAX));
    }
    let walk = json!({
        "q_au": WALK_Q,
        "dt_myr": WALK_DT_YR / 1.0e6,
        "a_au": a_walk,
        "labels": a_walk
            .iter()
            .map(|&a| {
                let c = classify(&landscape, a, 1.0 - WALK_Q / a, m9, WIDTH_COEFF);
                class_name(c.class)
            })
            .collect::<Vec<_>>(),
    });

    json!({
        "p9": {"mass_earth": m9 / EARTH_MASS_SOLAR, "a_au": A_P9_AU, "e": E_P9},
        "walk": walk,
        "a_range_au": [A_MIN, A_MAX],
        "q_range_au": [Q_MIN, Q_MAX],
        "n_population": all.n,
        "resonant": all.resonant,
        "hopping": all.hopping,
        "non_resonant": all.non_resonant,
        "bins": bins,
        "sample": sample,
        "simple_resonances": simple,
        "n_spectrum": landscape.spectrum.len(),
        "overlap": {"a_au": a_grid, "k": overlap},
        "kicks": kicks,
        "kick_crossing_time_myr": kicks[1]["crossing_time_myr"],
        "median_simple_spacing_au": median_spacing,
    })
}
