//! Film export for `p9-2018-resonance`: the numbers its scene and ledger entry draw.
//!
//! Bailey, Brown & Batygin (2018) ask which commensurabilities a scattered
//! disk actually occupies under an eccentric Planet Nine, and whether the
//! observed objects can then pin down Planet Nine's semimajor axis. This export
//! runs the crate's planar N-body census at reduced scale (seven Planet Nine
//! eccentricities in parallel, 0.4 Myr instead of 4 Gyr), classifies each
//! particle by libration of its resonant angle, and computes the crate's
//! implied-a9 distributions for the six Batygin & Brown (2016) objects with
//! the Farey F5 catalogue and with the full catalogue.

use std::thread;

use p9_2018_resonance::probability_analysis::{
    compare_distributions, observed_kbo_axes, p_all_simple,
};
use p9_2018_resonance::resonance_catalog::{
    Resonance, extended_catalog, farey_f5, resonant_angle_from_angles,
};
use p9_2018_resonance::simulation::{
    ResonanceRunResult, ResonanceSimConfig, classify_particle, run_planar_simulation,
};
use p9_core::constants::YEAR_DAYS;
use p9_core::units::au;
use serde_json::{Value, json};

/// Planet Nine eccentricities run in parallel.
const E9: [f64; 7] = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7];
/// Particles per run.
const N_PARTICLES: usize = 100;
/// Run length (years).
const T_YR: f64 = 4.0e5;
/// Resonant-angle samples per run.
const N_SAMPLES: f64 = 150.0;

struct Run {
    e9: f64,
    a_p9: f64,
    result: ResonanceRunResult,
}

fn run(e9: f64, seed: u64) -> Run {
    let config = ResonanceSimConfig {
        n_particles: N_PARTICLES,
        t_total: T_YR * YEAR_DAYS,
        angle_sample_interval: T_YR / N_SAMPLES * YEAR_DAYS,
        ..ResonanceSimConfig::quick_test(e9)
    };
    let result = run_planar_simulation(&config, seed);
    Run {
        e9,
        a_p9: config.a_p9,
        result,
    }
}

/// Resonant angle history (degrees) of one particle in one resonance.
fn angle_series(run: &Run, idx: usize, res: &Resonance) -> (Vec<f64>, Vec<f64>) {
    let n = run.result.p9_samples.len();
    let mut t = Vec::new();
    let mut phi = Vec::new();
    for (k, s) in run.result.particle_series[idx].iter().enumerate() {
        if let Some(s) = s {
            let p9 = &run.result.p9_samples[k];
            let angle = resonant_angle_from_angles(s.lambda, p9.lambda, s.varpi, res);
            t.push(T_YR * k as f64 / (n - 1) as f64 / 1e3);
            phi.push(angle.to_degrees());
        }
    }
    (t, phi)
}

pub fn export() -> Value {
    let runs: Vec<Run> = thread::scope(|s| {
        let handles: Vec<_> = E9
            .iter()
            .enumerate()
            .map(|(k, &e)| s.spawn(move || run(e, 2018 + k as u64)))
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });

    let catalog = extended_catalog();
    let mut resonant = Vec::new();
    let mut n_total = 0usize;
    let mut n_simple = 0usize;
    let mut per_e = Vec::new();
    // One librating example of each kind, for the resonant-angle panel.
    let mut example_simple: Option<(usize, usize, Resonance, f64)> = None;
    let mut example_high: Option<(usize, usize, Resonance, f64)> = None;
    for (ri, r) in runs.iter().enumerate() {
        let (mut simple_e, mut res_e) = (0usize, 0usize);
        for (idx, series) in r.result.particle_series.iter().enumerate() {
            n_total += 1;
            let Some((res, amp)) =
                classify_particle(series, &r.result.p9_samples, r.a_p9, &catalog)
            else {
                continue;
            };
            let live: Vec<f64> = series.iter().flatten().map(|s| s.a).collect();
            let a_mean = live.iter().sum::<f64>() / live.len() as f64;
            res_e += 1;
            if res.is_simple() {
                simple_e += 1;
            }
            let slot = if res.is_simple() {
                &mut example_simple
            } else {
                &mut example_high
            };
            if slot.is_none_or(|(_, _, _, best)| amp < best) {
                *slot = Some((ri, idx, res, amp));
            }
            resonant.push(json!({
                "e9": r.e9,
                "p": res.p,
                "q": res.q,
                "a_mean_au": a_mean,
                "a_res_au": (res.semimajor_axis_typed(r.a_p9) / au(1.0)).value,
                "amplitude_deg": amp.to_degrees(),
                "simple": res.is_simple(),
            }));
        }
        n_simple += simple_e;
        per_e.push(json!({
            "e9": r.e9,
            "n_particles": r.result.total_count,
            "n_survivors": r.result.active_count,
            "n_resonant": res_e,
            "n_simple": simple_e,
        }));
    }
    let n_resonant = resonant.len();
    let p_simple = n_simple as f64 / n_resonant.max(1) as f64;

    let example = |slot: Option<(usize, usize, Resonance, f64)>| -> Value {
        match slot {
            Some((ri, idx, res, amp)) => {
                let (t, phi) = angle_series(&runs[ri], idx, &res);
                json!({"p": res.p, "q": res.q, "e9": runs[ri].e9,
                       "amplitude_deg": amp.to_degrees(), "t_kyr": t, "phi_deg": phi})
            }
            None => Value::Null,
        }
    };

    // Resonance locations for Planet Nine at 600 AU inside the particle range.
    let a_p9 = runs[0].a_p9;
    let locations = |cat: &[Resonance]| -> Vec<Value> {
        cat.iter()
            .filter(|r| r.period_ratio() < 1.0)
            .map(|r| (r, (r.semimajor_axis_typed(a_p9) / au(1.0)).value))
            .filter(|(_, a)| (100.0..=600.0).contains(a))
            .map(|(r, a)| json!({"p": r.p, "q": r.q, "a_au": a, "simple": r.is_simple()}))
            .collect()
    };

    let cmp = compare_distributions(&observed_kbo_axes());

    json!({
        "a_p9": a_p9,
        "t_kyr": T_YR / 1e3,
        "per_e9": per_e,
        "n_total": n_total,
        "n_resonant": n_resonant,
        "n_simple": n_simple,
        "p_simple": p_simple,
        "p_all6": p_all_simple(p_simple, 6),
        "resonant": resonant,
        "example_simple": example(example_simple),
        "example_high": example(example_high),
        "catalog_f5": locations(&farey_f5()),
        "catalog_full": locations(&catalog),
        "kbo_axes": observed_kbo_axes(),
        "a9_dist": {
            "a9_au": cmp.f5_distribution.bins,
            "f5": cmp.f5_distribution.density,
            "full": cmp.ext_distribution.density,
        },
        "f5_peak_a9": cmp.f5_peak_a9,
        "f5_peak_to_mean": cmp.f5_peak_to_mean,
        "full_peak_a9": cmp.ext_peak_a9,
        "full_peak_to_mean": cmp.ext_peak_to_mean,
    })
}
