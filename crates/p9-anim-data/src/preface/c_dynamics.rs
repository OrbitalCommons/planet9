//! Film export for `scenes/preface/preface_c_dynamics.py`: the numbers its scenes draw.
//!
//! - `clustering`: the Brown (2017) ETNO sample, its mean resultant length and
//!   a Monte Carlo of how often ten random orbit directions line up as well.
//! - `precession`: the giant planets' J2 apsidal precession, as a curve against
//!   semi-major axis and applied to the real ETNOs started in perfect alignment.
//! - `resonance`: Pluto in Neptune's rotating frame, at the centre of the 3:2
//!   resonance and at the opposite (unprotected) phase.

use std::f64::consts::{PI, TAU};

use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

use p9_2017_dynamics::hamiltonian::{compute_j2_effective, precession_rate_j2_per_year};
use p9_core::analysis::circular::{mean_resultant_length, rayleigh_p_value};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::types::solve_kepler;

const A_NEPTUNE: f64 = 30.07;
const A_PLUTO: f64 = 39.48;
const E_PLUTO: f64 = 0.2488;
const AGE_MYR: f64 = 4500.0;

fn round(v: f64, k: i32) -> f64 {
    let s = 10f64.powi(k);
    (v * s).round() / s
}

fn orbital_period_yr(a: f64) -> f64 {
    a.powf(1.5)
}

/// Giant-planet J2 precession period (Myr) of an orbit (a, e).
fn precession_period_myr(a: f64, e: f64, j2: f64) -> f64 {
    TAU / precession_rate_j2_per_year(a, e, j2) / 1e6
}

fn clustering() -> Value {
    let varpis: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| o.longitude_of_perihelion())
        .collect();
    let n = varpis.len();
    let r_obs = mean_resultant_length(&varpis);

    // How often do N uniformly random directions reach the observed R-bar?
    let trials = 200_000;
    let bins = 40;
    let mut hist = vec![0u32; bins];
    let mut beat = 0u32;
    let mut rng = rand::rngs::StdRng::seed_from_u64(2017);
    let mut draw = vec![0.0; n];
    for _ in 0..trials {
        for d in draw.iter_mut() {
            *d = rng.gen_range(0.0..TAU);
        }
        let r = mean_resultant_length(&draw);
        hist[((r * bins as f64) as usize).min(bins - 1)] += 1;
        if r >= r_obs {
            beat += 1;
        }
    }
    // Three example "random skies" for the head-to-tail demonstration.
    let examples: Vec<Value> = (0..3)
        .map(|_| {
            let v: Vec<f64> = (0..n).map(|_| rng.gen_range(0.0..TAU)).collect();
            json!({
                "varpi_deg": v.iter().map(|x| round(x.to_degrees(), 2)).collect::<Vec<_>>(),
                "r_bar": round(mean_resultant_length(&v), 3),
            })
        })
        .collect();

    json!({
        "objects": BROWN_2017_SAMPLE.iter().map(|o| json!({
            "name": o.name,
            "a": o.a,
            "e": o.e,
            "q": round(o.perihelion(), 1),
            "varpi_deg": round(o.longitude_of_perihelion().to_degrees(), 2),
        })).collect::<Vec<_>>(),
        "r_bar": round(r_obs, 3),
        "rayleigh_p": rayleigh_p_value(&varpis),
        "mc_trials": trials,
        "mc_hist": hist,
        "mc_frac_beat": beat as f64 / trials as f64,
        "random_examples": examples,
    })
}

fn precession() -> Value {
    let j2 = compute_j2_effective();
    let q_ref = 40.0;
    let a_curve: Vec<f64> = (0..=80)
        .map(|k| 10f64.powf(2.0 + k as f64 / 80.0))
        .collect();
    let per_curve: Vec<f64> = a_curve
        .iter()
        .map(|&a| round(precession_period_myr(a, 1.0 - q_ref / a, j2), 2))
        .collect();

    let objs = &BROWN_2017_SAMPLE;
    let rates: Vec<f64> = objs
        .iter()
        .map(|o| precession_rate_j2_per_year(o.a, o.e, j2) * 1e6) // rad/Myr
        .collect();
    // Start every real ETNO at the same apse direction and let each precess
    // at its own rate: how fast does the alignment decohere?
    let n_t = 181;
    let t_myr: Vec<f64> = (0..n_t)
        .map(|k| AGE_MYR * k as f64 / (n_t - 1) as f64)
        .collect();
    let r_bar: Vec<f64> = t_myr
        .iter()
        .map(|&t| {
            let v: Vec<f64> = rates.iter().map(|r| r * t).collect();
            round(mean_resultant_length(&v), 3)
        })
        .collect();
    // First time the ensemble R-bar falls below 0.5 (fine scan).
    let t_half = (1..=45_000)
        .map(|k| k as f64 * 0.1)
        .find(|&t| {
            let v: Vec<f64> = rates.iter().map(|r| r * t).collect();
            mean_resultant_length(&v) < 0.5
        })
        .unwrap_or(AGE_MYR);
    // Afterwards, how much of the time does chance re-alignment reach the
    // observed sample's R-bar?
    let r_obs = mean_resultant_length(
        &objs
            .iter()
            .map(|o| o.longitude_of_perihelion())
            .collect::<Vec<_>>(),
    );
    let later: Vec<f64> = (0..=((AGE_MYR - t_half) / 0.5) as usize)
        .map(|k| t_half + 0.5 * k as f64)
        .collect();
    let above = later
        .iter()
        .filter(|&&t| {
            let v: Vec<f64> = rates.iter().map(|r| r * t).collect();
            mean_resultant_length(&v) >= r_obs
        })
        .count();

    let sedna = &objs[0];
    json!({
        "giants": [["Jupiter", 5.203], ["Saturn", 9.537], ["Uranus", 19.189], ["Neptune", A_NEPTUNE]],
        "curve_q_au": q_ref,
        "curve_a_au": a_curve.iter().map(|a| round(*a, 2)).collect::<Vec<_>>(),
        "curve_period_myr": per_curve,
        "objects": objs.iter().zip(&rates).map(|(o, r)| json!({
            "name": o.name,
            "a": o.a,
            "e": o.e,
            "rate_deg_per_myr": round(r.to_degrees(), 4),
            "period_myr": round(TAU / r, 1),
        })).collect::<Vec<_>>(),
        "decohere_t_myr": t_myr.iter().map(|t| round(*t, 1)).collect::<Vec<_>>(),
        "decohere_r_bar": r_bar,
        "t_half_myr": round(t_half, 1),
        "r_bar_observed": round(r_obs, 3),
        "frac_time_above_observed": round(above as f64 / later.len() as f64, 4),
        "age_myr": AGE_MYR,
        "sedna": {
            "orbit_yr": round(orbital_period_yr(sedna.a), 0),
            "precession_myr": round(precession_period_myr(sedna.a, sedna.e, j2), 1),
            "orbits_per_turn": round(precession_period_myr(sedna.a, sedna.e, j2) * 1e6 / orbital_period_yr(sedna.a), -2),
        },
    })
}

/// Pluto's path seen from a frame rotating with Neptune, over one full 3:2
/// cycle (2 Pluto orbits = 3 Neptune orbits), for resonant phase `phi`
/// (φ = 3λ_P − 2λ_N − ϖ_P). Pluto's perihelion is on +x at t = 0.
fn pluto_rotating(phi: f64, n: usize) -> Value {
    let p_n = orbital_period_yr(A_NEPTUNE);
    let n_n = TAU / p_n;
    let n_p = 2.0 / 3.0 * n_n; // exact commensurability at the resonance centre
    let lam_n0 = -phi / 2.0; // φ = 3·0 − 2λ_N0 − 0
    let t_end = 3.0 * p_n;
    let (mut xs, mut ys, mut d_min, mut r_at_conj) = (Vec::new(), Vec::new(), f64::MAX, 0.0);
    let mut sep_min = f64::MAX;
    for k in 0..n {
        let t = t_end * k as f64 / (n - 1) as f64;
        let ea = solve_kepler(E_PLUTO, n_p * t);
        let x = A_PLUTO * (ea.cos() - E_PLUTO);
        let y = A_PLUTO * (1.0 - E_PLUTO * E_PLUTO).sqrt() * ea.sin();
        let lam_n = lam_n0 + n_n * t;
        // rotate into Neptune's frame: Neptune fixed on +x at (a_N, 0)
        let (s, c) = (-lam_n).sin_cos();
        let (xr, yr) = (c * x - s * y, s * x + c * y);
        xs.push(round(xr, 3));
        ys.push(round(yr, 3));
        let d = ((xr - A_NEPTUNE).powi(2) + yr * yr).sqrt();
        if d < d_min {
            d_min = d;
        }
        // Conjunction: Pluto's longitude equals Neptune's (yr ≈ 0, xr > 0).
        let sep = yr.atan2(xr).abs();
        if sep < sep_min {
            sep_min = sep;
            r_at_conj = (xr * xr + yr * yr).sqrt();
        }
    }
    json!({
        "phi_deg": round(phi.to_degrees(), 1),
        "x": xs,
        "y": ys,
        "min_dist_au": round(d_min, 2),
        "r_at_conjunction_au": round(r_at_conj, 1),
    })
}

fn resonance() -> Value {
    json!({
        "neptune": { "a": A_NEPTUNE, "period_yr": round(orbital_period_yr(A_NEPTUNE), 1) },
        "pluto": {
            "a": A_PLUTO,
            "e": E_PLUTO,
            "q": round(A_PLUTO * (1.0 - E_PLUTO), 2),
            "Q": round(A_PLUTO * (1.0 + E_PLUTO), 2),
            "period_yr": round(orbital_period_yr(A_PLUTO), 1),
        },
        "period_ratio": round(orbital_period_yr(A_PLUTO) / orbital_period_yr(A_NEPTUNE), 4),
        "protected": pluto_rotating(PI, 721),
        "unprotected": pluto_rotating(0.0, 721),
    })
}

pub fn export() -> Value {
    json!({
        "clustering": clustering(),
        "precession": precession(),
        "resonance": resonance(),
    })
}
