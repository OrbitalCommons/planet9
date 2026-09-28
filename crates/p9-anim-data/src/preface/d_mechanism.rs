//! Film export for `scenes/preface/preface_d_mechanism.py`: the numbers its scenes draw.
//!
//! The mechanism scenes show how an eccentric Planet Nine confines the apsides
//! of the distant belt. Everything here is the coplanar, doubly-averaged
//! (secular) problem of Batygin & Brown (2016): Planet Nine smeared into a
//! Gauss ring (`numerical_secular_hamiltonian`, exact 1/Δ, no multipole
//! expansion) plus the four giant planets smeared into an effective J2. At a
//! fixed semi-major axis this has one degree of freedom, (Δϖ, G) with
//! G = √(GMa(1−e²)), and Hamilton's equations
//!
//!   dΔϖ/dt = ∂H/∂G,    dG/dt = −∂H/∂Δϖ
//!
//! are integrated with RK4 on a tabulated H(e, Δϖ). Neptune's scattering is
//! represented by removing any orbit whose perihelion dips inside 30 AU.

use std::f64::consts::{PI, TAU};

use rand::{Rng, SeedableRng};
use serde_json::{Value, json};

use p9_2017_dynamics::hamiltonian::compute_j2_effective;
use p9_core::analysis::circular::{circular_mean, mean_resultant_length};
use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::{GM_SUN, YEAR_DAYS};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::types::{P9Params, solve_kepler};

const DAYS_PER_MYR: f64 = 1e6 * YEAR_DAYS;
const A_NEPTUNE: f64 = 30.07;
/// Semi-major axis of the belt orbit used in the single-orbit beats (AU).
const A_BELT: f64 = 400.0;
/// Eccentricity of the belt orbit used in the ring-energy beat.
const E_BELT: f64 = 0.8;
/// Age of the solar system used for the swarm run (Myr).
const T_SWARM_MYR: f64 = 4000.0;

/// H(e, Δϖ) tabulated at one semi-major axis, with smooth interpolation.
struct Table {
    a: f64,
    e_lo: f64,
    de: f64,
    n_e: usize,
    n_w: usize,
    h: Vec<f64>, // [ie * n_w + iw]
}

impl Table {
    fn build(a: f64, p9: &P9Params, j2: f64) -> Self {
        let (e_lo, e_hi, n_e, n_w) = (0.01, 0.97, 97usize, 72usize);
        let de = (e_hi - e_lo) / (n_e - 1) as f64;
        let soft = 0.02 * p9.a;
        let mut h = vec![0.0; n_e * n_w];
        for ie in 0..n_e {
            let e = e_lo + ie as f64 * de;
            let h_j2 = giants_j2_energy(a, e, j2);
            // H is even in Δϖ (reflection symmetry about P9's apse line).
            for iw in 0..=n_w / 2 {
                let w = TAU * iw as f64 / n_w as f64;
                let ring = numerical_secular_hamiltonian(
                    a,
                    e,
                    0.0,
                    w,
                    0.0,
                    p9.a,
                    p9.e,
                    p9.gm(),
                    128,
                    soft,
                );
                h[ie * n_w + iw] = ring + h_j2;
                if iw > 0 && iw < n_w / 2 {
                    h[ie * n_w + (n_w - iw)] = ring + h_j2;
                }
            }
        }
        Self {
            a,
            e_lo,
            de,
            n_e,
            n_w,
            h,
        }
    }

    fn at(&self, ie: isize, iw: isize) -> f64 {
        let ie = ie.clamp(0, self.n_e as isize - 1) as usize;
        let iw = iw.rem_euclid(self.n_w as isize) as usize;
        self.h[ie * self.n_w + iw]
    }

    /// Catmull-Rom bicubic interpolation (periodic in Δϖ).
    fn eval(&self, e: f64, w: f64) -> f64 {
        let x = (e - self.e_lo) / self.de;
        let y = w.rem_euclid(TAU) / TAU * self.n_w as f64;
        let (ix, fx) = (x.floor() as isize, x - x.floor());
        let (iy, fy) = (y.floor() as isize, y - y.floor());
        let cr = |p0: f64, p1: f64, p2: f64, p3: f64, t: f64| {
            0.5 * (2.0 * p1
                + (-p0 + p2) * t
                + (2.0 * p0 - 5.0 * p1 + 4.0 * p2 - p3) * t * t
                + (-p0 + 3.0 * p1 - 3.0 * p2 + p3) * t * t * t)
        };
        let mut col = [0.0; 4];
        for (k, c) in col.iter_mut().enumerate() {
            let i = ix - 1 + k as isize;
            *c = cr(
                self.at(i, iy - 1),
                self.at(i, iy),
                self.at(i, iy + 1),
                self.at(i, iy + 2),
                fy,
            );
        }
        cr(col[0], col[1], col[2], col[3], fx)
    }

    /// Hamilton's equations in (Δϖ, e): (dΔϖ/dt, de/dt) per day.
    fn rhs(&self, w: f64, e: f64) -> (f64, f64) {
        let l = (GM_SUN * self.a).sqrt();
        let dg_de = -l * e / (1.0 - e * e).sqrt();
        let (he, hw) = (1e-4, 1e-3);
        let dh_de = (self.eval(e + he, w) - self.eval(e - he, w)) / (2.0 * he);
        let dh_dw = (self.eval(e, w + hw) - self.eval(e, w - hw)) / (2.0 * hw);
        (dh_de / dg_de, -dh_dw / dg_de)
    }
}

/// Orbit-averaged energy of the four giant planets smeared into a J2 ring,
/// −GM·J2/(2a³(1−e²)^{3/2}): it drives the prograde apsidal precession
/// (3/2)·n·J2/(a²(1−e²)²).
fn giants_j2_energy(a: f64, e: f64, j2: f64) -> f64 {
    -GM_SUN * j2 / (2.0 * a.powi(3) * (1.0 - e * e).powf(1.5))
}

/// Raw adaptive-step trajectory: (t days, unwrapped Δϖ, e).
struct Raw {
    t: Vec<f64>,
    w: Vec<f64>,
    e: Vec<f64>,
}

/// Adaptive RK4 (step limited to ~0.01 rad of Δϖ or ~0.0025 of e), stopping at
/// `t_end` days or, with `one_cycle`, as soon as the trajectory closes (Δϖ
/// crosses its starting value, mod 2π, moving the same way as at the start).
fn integrate(tab: &Table, w0: f64, e0: f64, t_end: f64, one_cycle: bool) -> Raw {
    let (mut t, mut w, mut e) = (0.0, w0, e0);
    let mut raw = Raw {
        t: vec![0.0],
        w: vec![w0],
        e: vec![e0],
    };
    let dir0 = tab.rhs(w0, e0).0.signum();
    let phase = |x: f64| (x - w0 + PI).rem_euclid(TAU) - PI;
    while t < t_end {
        let (k1w, k1e) = tab.rhs(w, e);
        let rate = k1w.abs() + 4.0 * k1e.abs();
        let dt = (0.01 / rate.max(1e-30)).min(t_end - t);
        let (k2w, k2e) = tab.rhs(w + 0.5 * dt * k1w, e + 0.5 * dt * k1e);
        let (k3w, k3e) = tab.rhs(w + 0.5 * dt * k2w, e + 0.5 * dt * k2e);
        let (k4w, k4e) = tab.rhs(w + dt * k3w, e + dt * k3e);
        let w_new = w + dt / 6.0 * (k1w + 2.0 * k2w + 2.0 * k3w + k4w);
        let e_new = (e + dt / 6.0 * (k1e + 2.0 * k2e + 2.0 * k3e + k4e)).clamp(0.012, 0.968);
        let crossed = one_cycle
            && raw.t.len() > 20
            && phase(w).signum() == -dir0
            && phase(w_new).signum() != -dir0
            && (phase(w_new) - phase(w)).abs() < PI;
        if crossed {
            // Close the loop exactly at the crossing.
            let f = phase(w).abs() / (phase(w).abs() + phase(w_new).abs());
            raw.t.push(t + f * dt);
            raw.w.push(w + f * (w_new - w));
            raw.e.push(e + f * (e_new - e));
            break;
        }
        (w, e) = (w_new, e_new);
        t += dt;
        raw.t.push(t);
        raw.w.push(w);
        raw.e.push(e);
    }
    raw
}

/// Linear resample of a raw trajectory on `n` uniform times over its span.
fn resample(raw: &Raw, n: usize) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let t_end = *raw.t.last().unwrap();
    let (mut ts, mut ws, mut es) = (Vec::new(), Vec::new(), Vec::new());
    let mut j = 0;
    for k in 0..n {
        let tk = t_end * k as f64 / (n - 1) as f64;
        while j + 1 < raw.t.len() - 1 && raw.t[j + 1] < tk {
            j += 1;
        }
        let span = (raw.t[j + 1] - raw.t[j]).max(1e-30);
        let f = ((tk - raw.t[j]) / span).clamp(0.0, 1.0);
        ts.push(tk / DAYS_PER_MYR);
        ws.push(raw.w[j] + f * (raw.w[j + 1] - raw.w[j]));
        es.push(raw.e[j] + f * (raw.e[j + 1] - raw.e[j]));
    }
    (ts, ws, es)
}

fn round(v: f64, k: i32) -> f64 {
    let s = 10f64.powi(k);
    (v * s).round() / s
}

fn deg(v: &[f64]) -> Vec<f64> {
    v.iter().map(|x| round(x.to_degrees(), 2)).collect()
}

fn r4(v: &[f64]) -> Vec<f64> {
    v.iter().map(|x| round(*x, 4)).collect()
}

fn span(v: &[f64]) -> (f64, f64) {
    v.iter()
        .fold((f64::MAX, f64::MIN), |(a, b), &x| (a.min(x), b.max(x)))
}

/// Planet Nine positions at uniform steps in time (mean anomaly), apse on +x:
/// the "time-density" of the smeared ring.
fn p9_ring(p9: &P9Params, n: usize) -> Value {
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for k in 0..n {
        let m = TAU * k as f64 / n as f64;
        let ea = solve_kepler(p9.e, m);
        xs.push(round(p9.a * (ea.cos() - p9.e), 2));
        ys.push(round(p9.a * (1.0 - p9.e * p9.e).sqrt() * ea.sin(), 2));
    }
    json!({ "x": xs, "y": ys })
}

/// The ring-ring interaction energy vs orientation for the belt orbit, and
/// the torque −∂H/∂Δϖ, evaluated directly with a fine quadrature (the giants'
/// J2 term is independent of Δϖ and drops out of both).
fn ring_energy(p9: &P9Params) -> Value {
    let h_at = |w: f64| {
        numerical_secular_hamiltonian(
            A_BELT,
            E_BELT,
            0.0,
            w,
            0.0,
            p9.a,
            p9.e,
            p9.gm(),
            192,
            0.02 * p9.a,
        )
    };
    let ws: Vec<f64> = (0..=240).map(|k| TAU * k as f64 / 240.0).collect();
    let h: Vec<f64> = ws.iter().map(|&w| h_at(w)).collect();
    let mean = h.iter().sum::<f64>() / h.len() as f64;
    let amp = h.iter().map(|x| (x - mean).abs()).fold(0.0, f64::max);
    let hw = 0.5_f64.to_radians();
    let torque: Vec<f64> = ws
        .iter()
        .map(|&w| -(h_at(w + hw) - h_at(w - hw)) / (2.0 * hw))
        .collect();
    let tamp = torque.iter().map(|x| x.abs()).fold(0.0, f64::max);
    json!({
        "a": A_BELT,
        "e": E_BELT,
        "dvarpi_deg": deg(&ws),
        "h_norm": r4(&h.iter().map(|x| (x - mean) / amp).collect::<Vec<_>>()),
        "torque_norm": r4(&torque.iter().map(|x| x / tamp).collect::<Vec<_>>()),
    })
}

/// A family of closed secular trajectories at `A_BELT`: the phase portrait.
fn portrait(tab: &Table) -> Value {
    let starts = [
        (180.0, 0.62),
        (180.0, 0.70),
        (180.0, 0.76),
        (180.0, 0.86),
        (180.0, 0.90),
        (180.0, 0.94),
        (0.0, 0.15),
        (0.0, 0.30),
        (0.0, 0.45),
        (0.0, 0.62),
        (0.0, 0.75),
        (0.0, 0.85),
    ];
    let t_cap = 20_000.0 * DAYS_PER_MYR;
    let tracks: Vec<Value> = starts
        .iter()
        .map(|&(w0, e0)| {
            let raw = integrate(tab, f64::to_radians(w0), e0, t_cap, true);
            let (ts, ws, es) = resample(&raw, 240);
            let (wlo, whi) = span(&ws);
            let (elo, ehi) = span(&es);
            let circulates = whi - wlo >= TAU - 0.05;
            json!({
                "start_deg": w0,
                "e0": e0,
                "kind": if circulates { "circulation" } else { "libration" },
                "center_deg": if circulates { Value::Null } else {
                    json!(round((0.5 * (wlo + whi)).to_degrees(), 1).rem_euclid(360.0) + 0.0)
                },
                "period_myr": round(*ts.last().unwrap(), 1),
                "q_min": round(A_BELT * (1.0 - ehi), 1),
                "q_max": round(A_BELT * (1.0 - elo), 1),
                "t_myr": r4(&ts),
                "dvarpi_deg": deg(&ws),
                "e": r4(&es),
            })
        })
        .collect();
    json!({
        "a": A_BELT,
        "e_neptune": 1.0 - A_NEPTUNE / A_BELT,
        "tracks": tracks,
    })
}

/// One swarm orbit: semi-major axis, sampled (Δϖ, e), and the sample at which
/// its perihelion first reached Neptune (if it did).
struct SwarmTrack {
    a: f64,
    w: Vec<f64>,
    e: Vec<f64>,
    removed: Option<usize>,
}

/// A scattered disk (q = 32–50 AU, random orientation) evolved for 4 Gyr at a
/// handful of semi-major axes. Orbits whose perihelion reaches Neptune are
/// marked removed at that sample.
fn swarm(p9: &P9Params, j2: f64) -> Value {
    let a_vals = [250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0];
    let per_a = 6;
    let n_out = 201;
    let mut rng = rand::rngs::StdRng::seed_from_u64(2016);
    let mut particles = Vec::new();
    let mut tracks: Vec<SwarmTrack> = Vec::new();
    for &a in &a_vals {
        let tab = Table::build(a, p9, j2);
        for _ in 0..per_a {
            let q0: f64 = rng.gen_range(32.0..50.0);
            let w0: f64 = rng.gen_range(0.0..TAU);
            let raw = integrate(&tab, w0, 1.0 - q0 / a, T_SWARM_MYR * DAYS_PER_MYR, false);
            let (_, ws, es) = resample(&raw, n_out);
            let removed = es.iter().position(|&e| a * (1.0 - e) < A_NEPTUNE);
            tracks.push(SwarmTrack {
                a,
                w: ws,
                e: es,
                removed,
            });
        }
    }
    let t_myr: Vec<f64> = (0..n_out)
        .map(|k| round(T_SWARM_MYR * k as f64 / (n_out - 1) as f64, 2))
        .collect();
    let alive_at = |k: usize, r: &Option<usize>| r.is_none_or(|i| k < i);
    let mut rbar = Vec::new();
    let mut mean_dw = Vec::new();
    let mut alive = Vec::new();
    for k in 0..n_out {
        let ws: Vec<f64> = tracks
            .iter()
            .filter(|t| alive_at(k, &t.removed))
            .map(|t| t.w[k])
            .collect();
        alive.push(ws.len());
        rbar.push(round(mean_resultant_length(&ws), 3));
        mean_dw.push(round(
            circular_mean(&ws)
                .unwrap_or(0.0)
                .to_degrees()
                .rem_euclid(360.0),
            1,
        ));
    }
    let q_at = |k: usize| -> Vec<f64> {
        let mut q: Vec<f64> = tracks
            .iter()
            .filter(|t| alive_at(k, &t.removed))
            .map(|t| t.a * (1.0 - t.e[k]))
            .collect();
        q.sort_by(|a, b| a.partial_cmp(b).unwrap());
        q
    };
    let median = |v: Vec<f64>| if v.is_empty() { 0.0 } else { v[v.len() / 2] };
    for t in &tracks {
        particles.push(json!({
            "a": t.a,
            "dvarpi_deg": deg(&t.w),
            "e": r4(&t.e),
            "removed_at": t.removed,
        }));
    }
    json!({
        "t_myr": t_myr,
        "particles": particles,
        "alive": alive,
        "r_bar": rbar,
        "mean_dvarpi_deg": mean_dw,
        "median_q_start": round(median(q_at(0)), 1),
        "median_q_end": round(median(q_at(n_out - 1)), 1),
        "q_neptune": A_NEPTUNE,
    })
}

pub fn export() -> Value {
    let p9 = P9Params::nominal_2016();
    let j2 = compute_j2_effective();
    let tab = Table::build(A_BELT, &p9, j2);
    let varpi9 = (p9.omega + p9.omega_big).rem_euclid(TAU);
    let etno_dvarpi: Vec<f64> = BROWN_2017_SAMPLE
        .iter()
        .map(|o| (o.longitude_of_perihelion() - varpi9).rem_euclid(TAU))
        .collect();
    json!({
        "p9": {
            "mass_earth": p9.mass_earth,
            "a": p9.a,
            "e": p9.e,
            "period_yr": round(p9.a.powf(1.5), 0),
        },
        "neptune_a": A_NEPTUNE,
        "ring": p9_ring(&p9, 90),
        "energy": ring_energy(&p9),
        "portrait": portrait(&tab),
        "swarm": swarm(&p9, j2),
        "etno_dvarpi_deg": deg(&etno_dvarpi),
        "etno_r_bar": round(mean_resultant_length(&etno_dvarpi), 3),
    })
}
