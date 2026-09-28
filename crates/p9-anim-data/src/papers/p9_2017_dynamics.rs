//! Film export for `p9-2017-dynamics`: the numbers its scene and ledger entry draw.

use p9_2017_dynamics::hamiltonian::{SecularHamiltonianParams, compute_j2_effective};
use p9_2017_dynamics::resonance::{
    FIDUCIAL_E9, Resonance, critical_period_ratio, critical_period_ratio_fiducial, period_ratio,
};
use p9_core::analysis::secular::phase_portrait;
use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};
use serde_json::{Value, json};

/// Shared perihelion of the scattered-disk chain (the crate's fiducial).
const Q_PERI: f64 = 33.0;

/// Semi-major axis of the secular portrait (AU): P/P9 = 0.28, well inside
/// the resonance-dominated regime the paper identifies.
const A_PORTRAIT: f64 = 300.0;

/// Whether two coplanar, confocal Kepler ellipses intersect: the particle
/// orbit (a, e) with its perihelion `dv` radians from Planet Nine's, and
/// Planet Nine's orbit (a9, e9). Sampled on a fine grid of true longitude.
pub fn orbits_cross(a: f64, e: f64, dv: f64, a9: f64, e9: f64) -> bool {
    let (p, p9) = (a * (1.0 - e * e), a9 * (1.0 - e9 * e9));
    let gap = |th: f64| p / (1.0 + e * (th - dv).cos()) - p9 / (1.0 + e9 * th.cos());
    let first = gap(0.0).signum();
    (1..720).any(|k| gap(std::f64::consts::TAU * k as f64 / 720.0).signum() != first)
}

/// The paper's coplanar secular portrait H(e, Δϖ) at semi-major axis `a`:
/// Planet Nine's exactly ring-averaged potential plus the giant planets'
/// orbit-averaged J2 field.
pub fn portrait(a: f64, params: &SecularHamiltonianParams) -> Value {
    let gm9 = params.m9_solar * GM_SUN;
    let (e_vals, dv_vals, grid) = phase_portrait(a, params.a9, params.e9, gm9, 40, 72);
    let j2 = compute_j2_effective();
    let h: Vec<Vec<f64>> = e_vals
        .iter()
        .zip(&grid)
        .map(|(&e, row)| {
            // Orbit-averaged J2 potential −GM J2R²/(2a³η³), the Hamiltonian
            // whose apsidal rate is the crate's `precession_rate_j2`.
            let eta3 = (1.0 - e * e).powf(1.5);
            let h_j2 = -0.5 * GM_SUN * j2 / (a * a * a * eta3);
            row.iter().map(|&h9| h9 + h_j2).collect()
        })
        .collect();
    let crossing: Vec<Vec<bool>> = e_vals
        .iter()
        .map(|&e| {
            dv_vals
                .iter()
                .map(|&dv| orbits_cross(a, e, dv, params.a9, params.e9))
                .collect()
        })
        .collect();
    json!({
        "a": a,
        "e": e_vals,
        "dvarpi_deg": dv_vals.iter().map(|v| v.to_degrees()).collect::<Vec<_>>(),
        "h": h,
        "crossing": crossing,
    })
}

pub fn export() -> Value {
    let params = SecularHamiltonianParams::default_paper();
    let (a9, e9, m9) = (params.a9, params.e9, params.m9_solar);
    let q9 = a9 * (1.0 - e9);

    // The N/1 resonance chain of scattered-disk objects sharing q = 33 AU:
    // pendulum widths against the spacing to the next member, and whether
    // the member's aphelion reaches Planet Nine's perihelion.
    let chain: Vec<Value> = (2..=16i64)
        .map(|j| {
            let a_j = Resonance::nominal_semimajor(j, 1, a9);
            let e_j = 1.0 - Q_PERI / a_j;
            let aphelion = 2.0 * a_j - Q_PERI;
            let half_width = Resonance::pendulum_half_width(j, 1, e_j, m9, a9, FIDUCIAL_E9);
            let spacing = a_j - Resonance::nominal_semimajor(j + 1, 1, a9);
            json!({
                "j": j,
                "a": a_j,
                "e": e_j,
                "period_ratio": period_ratio(a_j, a9),
                "aphelion": aphelion,
                "half_width": half_width,
                "spacing": spacing,
                "reaches_p9": aphelion >= q9,
            })
        })
        .collect();

    let critical = critical_period_ratio_fiducial();
    let alt = SecularHamiltonianParams::alternative_600();
    let critical_600 = critical_period_ratio(Q_PERI, alt.m9_solar, alt.a9, alt.e9);
    let by_mass: Vec<Value> = [5.0, 10.0, 20.0]
        .iter()
        .map(|&m| {
            json!({
                "mass_earth": m,
                "critical_ratio": critical_period_ratio(Q_PERI, m * EARTH_MASS_SOLAR, a9, e9),
            })
        })
        .collect();

    json!({
        "a9": a9,
        "e9": e9,
        "m9_earth": m9 / EARTH_MASS_SOLAR,
        "q9": q9,
        "q_peri": Q_PERI,
        "portrait": portrait(A_PORTRAIT, &params),
        "chain": chain,
        "critical_ratio": critical,
        "critical_a": a9 * critical.powf(2.0 / 3.0),
        "critical_ratio_600": critical_600,
        "critical_by_mass": by_mass,
    })
}
