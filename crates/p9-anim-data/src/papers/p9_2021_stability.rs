//! Film export for `p9-2021-stability`: the numbers its scene and ledger entry draw.

use p9_2021_stability::chirikov::OVERLAP_GLOBAL_CHAOS;
use p9_2021_stability::diffusion::measurement_from_ensemble;
use p9_2021_stability::nbody_validation::neptune_scattering_series;
use p9_2021_stability::resonance_chain::{build_resonance_chain, resonance_spacing};
use p9_2021_stability::stability::lyapunov_time_days;
use p9_core::analysis::resonance::{
    chirikov_overlap_parameter, critical_perihelion, neptune_diffusion_coefficient,
};
use p9_core::constants::YEAR_DAYS;
use p9_core::data::etno::BROWN_2017_SAMPLE;
use serde_json::{Value, json};

/// Semi-major axis at which the headline critical perihelion is quoted, and
/// at which the resonance chain and the N-body particles are placed (AU).
const A_REF: f64 = 500.0;
/// Perihelia on either side of the boundary at `A_REF` (AU).
const Q_INSIDE: f64 = 33.0;
const Q_OUTSIDE: f64 = 50.0;
/// Integrator step (days) and snapshot cadence (steps) of the N-body check.
const DT_DAYS: f64 = 2500.0;
const SNAPSHOT_EVERY: usize = 40;
/// Snapshots skipped between drawn points of a track.
const THIN: usize = 12;
/// 2:j resonances bracketing `A_REF` (2:136 sits at 500.6 AU).
const J_MIN: i64 = 131;
const J_MAX: i64 = 141;

/// The 2:j chain around `A_REF` at perihelion `q`: centres and half-widths.
fn chain(q: f64) -> Value {
    let members: Vec<Value> = build_resonance_chain(J_MIN, J_MAX, q)
        .iter()
        .map(|r| json!({"j": r.j, "a_au": r.a_nominal, "half_width_au": r.delta_a}))
        .collect();
    json!({
        "q_au": q,
        "overlap": chirikov_overlap_parameter(A_REF, q),
        "resonances": members,
    })
}

/// a(t) of a few test particles started at (`A_REF`, q) under Sun + Neptune,
/// with the measured and the analytic diffusion coefficient.
fn nbody(q: f64, seed: u64) -> Value {
    let t_total = 1.0e6 * YEAR_DAYS;
    let series = neptune_scattering_series(A_REF, q, 8, t_total, DT_DAYS, SNAPSHOT_EVERY, seed);
    let dt_yr = DT_DAYS * SNAPSHOT_EVERY as f64 / YEAR_DAYS;
    let measured = measurement_from_ensemble(A_REF, q, &series, dt_yr).map(|m| m.d_measured);
    let tracks: Vec<Vec<f64>> = series
        .iter()
        .map(|s| s.iter().step_by(THIN).copied().collect())
        .collect();
    json!({
        "q_au": q,
        "dt_myr": THIN as f64 * dt_yr / 1.0e6,
        "a_au": tracks,
        "d_measured": measured,
        "d_analytic": neptune_diffusion_coefficient(q),
    })
}

pub fn export() -> Value {
    let a: Vec<f64> = (0..=90).map(|k| 250.0 + 15.0 * k as f64).collect();
    let q_crit: Vec<f64> = a.iter().map(|&x| critical_perihelion(x)).collect();

    // q at which the overlap parameter reaches the standard-map global-chaos
    // threshold (s = 0.63), found by bisection on the crate's own parameter.
    let q_onset: Vec<f64> = a
        .iter()
        .map(|&x| {
            let (mut lo, mut hi) = (0.0, 120.0);
            for _ in 0..50 {
                let mid = 0.5 * (lo + hi);
                if chirikov_overlap_parameter(x, mid) > OVERLAP_GLOBAL_CHAOS {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            0.5 * (lo + hi)
        })
        .collect();

    let objects: Vec<Value> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a_au": o.a,
                "q_au": o.perihelion(),
                "chaotic": o.perihelion() < critical_perihelion(o.a),
            })
        })
        .collect();
    let n_chaotic = BROWN_2017_SAMPLE
        .iter()
        .filter(|o| o.perihelion() < critical_perihelion(o.a))
        .count();

    json!({
        "a_ref_au": A_REF,
        "q_crit_headline_au": critical_perihelion(A_REF),
        "spacing_au": resonance_spacing((J_MIN + J_MAX) / 2),
        "lyapunov_time_yr": lyapunov_time_days(A_REF) / YEAR_DAYS,
        "boundary": {"a_au": a, "q_crit_au": q_crit, "q_onset_au": q_onset},
        "chains": [chain(Q_INSIDE), chain(Q_OUTSIDE)],
        "objects": objects,
        "n_objects": BROWN_2017_SAMPLE.len(),
        "n_chaotic": n_chaotic,
        "nbody": [nbody(Q_INSIDE, 11), nbody(Q_OUTSIDE, 11)],
    })
}
