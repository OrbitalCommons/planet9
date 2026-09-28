//! Film export for `p9-2016-secular-resonance`: the numbers its scene and ledger entry draw.

use p9_2016_secular_resonance::SecularModel;
use p9_2016_secular_resonance::published::nominal_p9;
use p9_2017_dynamics::hamiltonian::{SecularHamiltonianParams, octupole_term};
use p9_core::constants::EARTH_MASS_SOLAR;
use serde_json::{Value, json};

/// Semi-major axis of the portrait (AU): a distant orbit in the middle of the
/// observed clustered sample.
const A_PORTRAIT: f64 = 400.0;

/// Portrait grid.
const N_E: usize = 95;
const N_W: usize = 121;
const E_MIN: f64 = 0.02;
const E_MAX: f64 = 0.96;

/// Number of level curves drawn per portrait.
const N_LEVELS: usize = 16;

/// Level curves of a gridded function by marching squares. `z[k][j]` is the
/// value at `(x[j], y[k])`. Returns the curves as polylines of `[x, y]`.
pub fn level_curves(x: &[f64], y: &[f64], z: &[Vec<f64>], level: f64) -> Vec<Vec<[f64; 2]>> {
    let mut segments: Vec<([f64; 2], [f64; 2])> = Vec::new();
    let cross = |xa: f64, ya: f64, za: f64, xb: f64, yb: f64, zb: f64| -> [f64; 2] {
        let t = (level - za) / (zb - za);
        [xa + t * (xb - xa), ya + t * (yb - ya)]
    };
    for k in 0..y.len() - 1 {
        for j in 0..x.len() - 1 {
            let corners = [
                (x[j], y[k], z[k][j]),
                (x[j + 1], y[k], z[k][j + 1]),
                (x[j + 1], y[k + 1], z[k + 1][j + 1]),
                (x[j], y[k + 1], z[k + 1][j]),
            ];
            let mut hits = Vec::with_capacity(4);
            for c in 0..4 {
                let (xa, ya, za) = corners[c];
                let (xb, yb, zb) = corners[(c + 1) % 4];
                if (za < level) != (zb < level) {
                    hits.push(cross(xa, ya, za, xb, yb, zb));
                }
            }
            match hits.len() {
                2 => segments.push((hits[0], hits[1])),
                4 => {
                    segments.push((hits[0], hits[1]));
                    segments.push((hits[2], hits[3]));
                }
                _ => {}
            }
        }
    }

    // Chain the segments into polylines by matching shared end points.
    let close = |p: [f64; 2], q: [f64; 2]| (p[0] - q[0]).abs() < 1e-9 && (p[1] - q[1]).abs() < 1e-9;
    let mut used = vec![false; segments.len()];
    let mut curves = Vec::new();
    for s in 0..segments.len() {
        if used[s] {
            continue;
        }
        used[s] = true;
        let mut line = vec![segments[s].0, segments[s].1];
        loop {
            let mut grew = false;
            for t in 0..segments.len() {
                if used[t] {
                    continue;
                }
                let (p, q) = segments[t];
                let tail = *line.last().unwrap();
                let head = line[0];
                if close(tail, p) {
                    line.push(q);
                } else if close(tail, q) {
                    line.push(p);
                } else if close(head, q) {
                    line.insert(0, p);
                } else if close(head, p) {
                    line.insert(0, q);
                } else {
                    continue;
                }
                used[t] = true;
                grew = true;
            }
            if !grew {
                break;
            }
        }
        curves.push(line);
    }
    curves
}

/// Evenly spaced quantile levels of a grid, so the curves fill the portrait.
pub fn quantile_levels(z: &[Vec<f64>], n: usize) -> Vec<f64> {
    let mut all: Vec<f64> = z.iter().flatten().copied().collect();
    all.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (0..n)
        .map(|k| all[((k as f64 + 0.5) / n as f64 * all.len() as f64) as usize])
        .collect()
}

fn portrait(a: f64) -> Value {
    let p9 = nominal_p9();
    let model = SecularModel::new(a, &p9).with_giants_j2();
    let es: Vec<f64> = (0..N_E)
        .map(|k| E_MIN + (E_MAX - E_MIN) * k as f64 / (N_E - 1) as f64)
        .collect();
    let ws: Vec<f64> = (0..N_W)
        .map(|j| 360.0 * j as f64 / (N_W - 1) as f64)
        .collect();
    let h: Vec<Vec<f64>> = es
        .iter()
        .map(|&e| {
            ws.iter()
                .map(|&w| model.hamiltonian(e, w.to_radians()))
                .collect()
        })
        .collect();
    // Each level curve is classified by the range of apsidal angle it covers:
    // spanning 180 degrees without coming within 45 degrees of alignment
    // (anti-aligned), confined near 0 (aligned), or spanning the whole circle
    // (circulation).
    let mut lines = Vec::new();
    let mut anti_halfwidth: f64 = 0.0;
    for level in quantile_levels(&h, N_LEVELS) {
        for line in level_curves(&ws, &es, &h, level) {
            if line.len() < 8 {
                continue;
            }
            let lo = line.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
            let hi = line.iter().map(|p| p[0]).fold(f64::NEG_INFINITY, f64::max);
            let kind = if lo > 45.0 && hi < 315.0 && lo < 180.0 && hi > 180.0 {
                anti_halfwidth = anti_halfwidth.max(0.5 * (hi - lo));
                "anti"
            } else if lo <= 0.0 && hi >= 360.0 {
                "circulating"
            } else if hi < 120.0 || lo > 240.0 {
                "aligned"
            } else {
                "other"
            };
            lines.push(json!({"kind": kind, "points": line}));
        }
    }

    json!({
        "a_au": a,
        "alpha": a / p9.a,
        "lines": lines,
        "anti_halfwidth_deg": anti_halfwidth,
    })
}

pub fn export() -> Value {
    let p9 = nominal_p9();

    // Amplitude of the cos(apsidal angle) term: the exact ring average against
    // the truncated (octupole) expansion.
    let a_cut = A_PORTRAIT;
    let model = SecularModel::new(a_cut, &p9);
    let params = SecularHamiltonianParams {
        a9: p9.a,
        e9: p9.e,
        m9_solar: p9.mass_earth * EARTH_MASS_SOLAR,
        j2_eff: 0.0,
        precession_rate_9: 0.0,
    };
    let es: Vec<f64> = (0..=45).map(|k| 0.05 + 0.02 * k as f64).collect();
    let exact: Vec<f64> = es.iter().map(|&e| model.apsidal_harmonics(e).1).collect();
    let truncated: Vec<f64> = es
        .iter()
        .map(|&e| octupole_term(a_cut, e, 0.0, &params))
        .collect();
    let unit = exact
        .iter()
        .chain(truncated.iter())
        .fold(0.0_f64, |m, v| m.max(v.abs()));

    let worst = exact
        .iter()
        .zip(&truncated)
        .map(|(x, t)| ((t - x) / x).abs())
        .fold(0.0_f64, f64::max);
    let portrait = portrait(A_PORTRAIT);

    json!({
        "planet": {"mass_earth": p9.mass_earth, "a_au": p9.a, "e": p9.e},
        "island_halfwidth_deg": portrait["anti_halfwidth_deg"],
        "truncation_error": worst,
        "harmonic": {
            "a_au": a_cut,
            "e": es,
            "exact": exact.iter().map(|v| v / unit).collect::<Vec<_>>(),
            "truncated": truncated.iter().map(|v| v / unit).collect::<Vec<_>>(),
        },
        "portrait": portrait,
    })
}
