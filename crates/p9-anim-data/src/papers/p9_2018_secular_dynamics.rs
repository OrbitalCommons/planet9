//! Film export for `p9-2018-secular-dynamics`: the numbers its scene and ledger entry draw.
//!
//! Everything comes from the crate's coplanar secular model: the conserved
//! Hamiltonian H(e, Δϖ; a) (giants' J2 + Planet Nine's Gauss-ring average),
//! the giants' free apsidal precession rate, and the critical semi-major axis
//! that separates circulation from libration.

use p9_2018_secular_dynamics::published::LIBRATION_A_THRESHOLD_AU;
use p9_2018_secular_dynamics::secular_model::{
    critical_semimajor_axis_typed, hamiltonian, librates, nominal_p9, precession_cross_check,
};
use p9_core::constants::GM_SUN;
use p9_core::forces::j2_secular::combined_j2_jsu;
use p9_core::types::P9Params;
use p9_core::units::au;
use serde_json::{Value, json};
use std::f64::consts::PI;

/// Days per million years.
const DAYS_PER_MYR: f64 = 365.25e6;
/// Largest eccentricity drawn in the (k, h) plane.
const E_EDGE: f64 = 0.92;
/// Grid points per side of the Hamiltonian map.
const N_GRID: usize = 61;
/// Arc-length step of the level-curve tracer in the (k, h) plane.
const DS: f64 = 0.01;

/// Secular velocity (dk/dt, dh/dt) in 1/day from Hamilton's equations for the
/// crate's Hamiltonian, with (k, h) = e (cos Δϖ, sin Δϖ):
/// de/dt = (√(1−e²)/(n a² e)) ∂H/∂Δϖ,  dΔϖ/dt = −(√(1−e²)/(n a² e)) ∂H/∂e.
fn velocity(a: f64, k: f64, h: f64, p9: &P9Params) -> (f64, f64) {
    let e = k.hypot(h);
    let phi = h.atan2(k);
    let (de, dphi) = (1e-4, 1e-3);
    let h_e = (hamiltonian(a, e + de, phi, p9) - hamiltonian(a, e - de, phi, p9)) / (2.0 * de);
    let h_phi =
        (hamiltonian(a, e, phi + dphi, p9) - hamiltonian(a, e, phi - dphi, p9)) / (2.0 * dphi);
    let n = ((GM_SUN + combined_j2_jsu().2) / a.powi(3)).sqrt();
    let f = (1.0 - e * e).sqrt() / (n * a * a * e);
    let e_dot = f * h_phi;
    let phi_dot = -f * h_e;
    (
        phi.cos() * e_dot - e * phi.sin() * phi_dot,
        phi.sin() * e_dot + e * phi.cos() * phi_dot,
    )
}

/// Follow one secular trajectory (a level curve of H) from (e0, Δϖ0) with RK4
/// in arc length, accumulating the physical time. Stops after one closed
/// cycle or when the orbit leaves the drawn disc.
fn trace(a: f64, e0: f64, dvarpi0: f64, p9: &P9Params) -> (Vec<[f64; 3]>, bool) {
    let (k0, h0) = (e0 * dvarpi0.cos(), e0 * dvarpi0.sin());
    let unit = |k: f64, h: f64| {
        let (vx, vy) = velocity(a, k, h, p9);
        let s = vx.hypot(vy);
        (vx / s, vy / s, 1.0 / s)
    };
    let (mut k, mut h, mut t) = (k0, h0, 0.0);
    let mut path = vec![[k, h, 0.0]];
    for step in 0..2000 {
        let (a1, b1, c1) = unit(k, h);
        let (a2, b2, c2) = unit(k + 0.5 * DS * a1, h + 0.5 * DS * b1);
        let (a3, b3, c3) = unit(k + 0.5 * DS * a2, h + 0.5 * DS * b2);
        let (a4, b4, c4) = unit(k + DS * a3, h + DS * b3);
        k += DS / 6.0 * (a1 + 2.0 * a2 + 2.0 * a3 + a4);
        h += DS / 6.0 * (b1 + 2.0 * b2 + 2.0 * b3 + b4);
        t += DS / 6.0 * (c1 + 2.0 * c2 + 2.0 * c3 + c4);
        path.push([k, h, t / DAYS_PER_MYR]);
        if k.hypot(h) > E_EDGE {
            return (path, false);
        }
        if step > 20 && (k - k0).hypot(h - h0) < 0.6 * DS {
            return (path, true);
        }
    }
    (path, false)
}

/// The Hamiltonian on a square (k, h) grid (null outside the drawn disc) and
/// one traced trajectory from (e0, Δϖ0).
fn portrait(a: f64, e0: f64, dvarpi0_deg: f64, p9: &P9Params) -> Value {
    let axis: Vec<f64> = (0..N_GRID)
        .map(|i| -E_EDGE + 2.0 * E_EDGE * i as f64 / (N_GRID - 1) as f64)
        .collect();
    let grid: Vec<Vec<Option<f64>>> = axis
        .iter()
        .map(|&h| {
            axis.iter()
                .map(|&k| {
                    let e = k.hypot(h);
                    (e <= E_EDGE && e > 1e-3).then(|| hamiltonian(a, e, h.atan2(k), p9))
                })
                .collect()
        })
        .collect();

    let (path, closed) = trace(a, e0, dvarpi0_deg.to_radians(), p9);
    let dvarpi: Vec<f64> = path
        .iter()
        .map(|p| p[1].atan2(p[0]).to_degrees().rem_euclid(360.0))
        .collect();
    let circulates = dvarpi.iter().any(|&d| !(10.0..=350.0).contains(&d));
    let e_path: Vec<f64> = path.iter().map(|p| p[0].hypot(p[1])).collect();
    json!({
        "a_au": a,
        "e0": e0,
        "q_au": a * (1.0 - e0),
        "dvarpi0_deg": dvarpi0_deg,
        "axis": axis,
        "h_grid": grid,
        "path": {
            "k": path.iter().map(|p| p[0]).collect::<Vec<_>>(),
            "h": path.iter().map(|p| p[1]).collect::<Vec<_>>(),
            "t_myr": path.iter().map(|p| p[2]).collect::<Vec<_>>(),
        },
        "closed": closed,
        "circulates": circulates,
        "period_myr": path.last().map(|p| p[2]),
        "dvarpi_min_deg": dvarpi.iter().copied().fold(f64::INFINITY, f64::min),
        "dvarpi_max_deg": dvarpi.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        "e_min": e_path.iter().copied().fold(f64::INFINITY, f64::min),
        "e_max": e_path.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        "crate_librates": librates(a, p9),
    })
}

pub fn export() -> Value {
    let p9 = nominal_p9();
    let a_crit = critical_semimajor_axis_typed(&p9, 50.0, 480.0).map(|l| (l / au(1.0)).value);

    let a_crit_by_mass: Vec<Value> = [5.0, 10.0, 15.0, 20.0]
        .iter()
        .map(|&m| {
            let mut p = nominal_p9();
            p.mass_earth = m;
            json!({
                "mass_earth": m,
                "a_crit_au": critical_semimajor_axis_typed(&p, 50.0, 480.0).map(|l| (l / au(1.0)).value),
            })
        })
        .collect();

    // The giants' free apsidal precession period against semi-major axis.
    let a_grid: Vec<f64> = (0..=76).map(|k| 80.0 + 5.0 * k as f64).collect();
    let period_myr: Vec<f64> = a_grid
        .iter()
        .map(|&a| 2.0 * PI / precession_cross_check(a, 0.01).0 / DAYS_PER_MYR)
        .collect();

    // Two ETNOs with the same 60 AU perihelion: one inside, one outside a_crit.
    let portraits = vec![
        portrait(150.0, 0.6, 180.0, &p9),
        portrait(300.0, 0.8, 140.0, &p9),
    ];

    json!({
        "p9": {"mass_earth": p9.mass_earth, "a_au": p9.a, "e": p9.e},
        "a_crit_au": a_crit,
        "published_a_threshold_au": LIBRATION_A_THRESHOLD_AU,
        "a_crit_by_mass": a_crit_by_mass,
        "precession": {"a_au": a_grid, "period_myr": period_myr},
        "portraits": portraits,
    })
}
