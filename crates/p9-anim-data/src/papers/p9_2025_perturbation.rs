//! Film export for `p9-2025-perturbation`: the numbers its scene and ledger entry draw.

use p9_2025_perturbation::chirikov::overlap_parameter;
use p9_2025_perturbation::resonance::{
    Resonance, build_1j_chain, build_2j_chain, build_3j_chain, build_4j_chain,
};
use p9_2025_perturbation::stability_boundary::{
    BoundaryCurve, OrbitClass, classify_orbit, effective_boundary,
};
use p9_core::analysis::resonance::critical_perihelion;
use p9_core::constants::A_NEPTUNE_AU;
use p9_core::data::etno::BROWN_2017_SAMPLE;
use serde_json::{Value, json};

/// Semi-major-axis window of the resonance-comb panel (AU).
const A_LO: f64 = 185.0;
const A_HI: f64 = 215.0;
/// Perihelia at which the comb is drawn (AU).
const Q_COMB: [f64; 2] = [33.0, 45.0];
/// Semi-major axes at which the overlap onset is solved (AU).
const A_ONSET: [f64; 10] = [
    100.0, 125.0, 150.0, 175.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0,
];
/// Semi-major axis of the headline comparison (AU).
const A_HEADLINE: f64 = 200.0;

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 { a } else { gcd(b, a % b) }
}

/// The m:j chain (m = 1..=`m_max`) between `a_lo` and `a_hi` at perihelion
/// `q`, lowest terms only, sorted by semi-major axis.
fn comb(m_max: u32, a_lo: f64, a_hi: f64, q: f64) -> Vec<Resonance> {
    let mut all: Vec<Resonance> = Vec::new();
    for m in 1..=m_max {
        // a = a_N (j/m)^(2/3), so the window maps to a range of j.
        let j_of = |a: f64| m as f64 * (a / A_NEPTUNE_AU).powf(1.5);
        let (j_min, j_max) = (j_of(a_lo).floor() as i64, j_of(a_hi).ceil() as i64);
        let chain = match m {
            1 => build_1j_chain(j_min, j_max, q),
            2 => build_2j_chain(j_min, j_max, q),
            3 => build_3j_chain(j_min, j_max, q),
            _ => build_4j_chain(j_min, j_max, q),
        };
        all.extend(
            chain.into_iter().filter(|r| {
                gcd(r.j, r.m as i64) == 1 && r.a_nominal >= a_lo && r.a_nominal <= a_hi
            }),
        );
    }
    all.sort_by(|x, y| x.a_nominal.partial_cmp(&y.a_nominal).unwrap());
    all
}

/// Median overlap parameter of neighbouring resonances within ±10% of `a`.
fn local_overlap(m_max: u32, a: f64, q: f64) -> f64 {
    let members = comb(m_max, 0.9 * a, 1.1 * a, q);
    let mut k: Vec<f64> = members
        .windows(2)
        .map(|w| overlap_parameter(&w[0], &w[1]))
        .collect();
    if k.is_empty() {
        return 0.0;
    }
    k.sort_by(|x, y| x.partial_cmp(y).unwrap());
    k[k.len() / 2]
}

/// Perihelion at which the local overlap parameter reaches 1 (bisection).
fn overlap_onset(m_max: u32, a: f64) -> Option<f64> {
    let (mut lo, mut hi) = (25.0, 80.0);
    if local_overlap(m_max, a, lo) < 1.0 {
        return None;
    }
    for _ in 0..12 {
        let mid = 0.5 * (lo + hi);
        if local_overlap(m_max, a, mid) >= 1.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    Some(0.5 * (lo + hi))
}

fn class_name(c: OrbitClass) -> &'static str {
    match c {
        OrbitClass::Unstable => "unstable",
        OrbitClass::BoundedChaos => "bounded_chaos",
        OrbitClass::Stable => "stable",
    }
}

pub fn export() -> Value {
    let fits = BoundaryCurve::compute(50.0, 600.0, 111);
    let quadrupole: Vec<f64> = fits
        .a_values
        .iter()
        .map(|&a| critical_perihelion(a))
        .collect();

    let combs: Vec<Value> = Q_COMB
        .iter()
        .map(|&q| {
            let members: Vec<Value> = comb(4, A_LO, A_HI, q)
                .iter()
                .map(|r| {
                    json!({
                        "m": r.m,
                        "j": r.j,
                        "a_au": r.a_nominal,
                        "width_au": r.delta_a,
                    })
                })
                .collect();
            json!({
                "q_au": q,
                "resonances": members,
                "overlap_2j_only": local_overlap(2, 0.5 * (A_LO + A_HI), q),
                "overlap_all": local_overlap(4, 0.5 * (A_LO + A_HI), q),
            })
        })
        .collect();

    let solved: Vec<(Option<f64>, Option<f64>)> = std::thread::scope(|s| {
        let handles: Vec<_> = A_ONSET
            .iter()
            .map(|&a| s.spawn(move || (overlap_onset(2, a), overlap_onset(4, a))))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    let onset: Vec<Value> = A_ONSET
        .iter()
        .zip(&solved)
        .map(|(&a, &(quadrupole_only, all_chains))| {
            json!({
                "a_au": a,
                "q_quadrupole_chain_au": quadrupole_only,
                "q_all_chains_au": all_chains,
                "q_published_fit_au": effective_boundary(a).0,
            })
        })
        .collect();
    let headline = A_ONSET
        .iter()
        .zip(&solved)
        .find(|(a, _)| **a == A_HEADLINE)
        .and_then(|(_, s)| s.1);

    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a_au": o.a,
                "q_au": o.perihelion(),
                "class": class_name(classify_orbit(o.a, o.perihelion())),
            })
        })
        .collect();

    let crossing = fits
        .a_values
        .iter()
        .zip(fits.q_chaotic.iter().zip(&fits.q_diffusion))
        .find(|(_, (c, d))| d > c)
        .map(|(&a, _)| a);

    json!({
        "fits": {
            "a_au": fits.a_values,
            "q_chaotic_au": fits.q_chaotic,
            "q_diffusion_au": fits.q_diffusion,
            "q_comb_au": fits.q_comb,
            "q_instability_au": fits.q_effective,
            "q_quadrupole_2021_au": quadrupole,
        },
        "a_regime_change_au": crossing,
        "comb_window_au": [A_LO, A_HI],
        "combs": combs,
        "onset": onset,
        "a_headline_au": A_HEADLINE,
        "q_onset_headline_au": headline,
        "q_published_headline_au": effective_boundary(A_HEADLINE).0,
        "objects": objects,
    })
}
