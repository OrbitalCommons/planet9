//! Film export for `p9-2020-secular-octupole`: the numbers its scene and ledger entry draw.
//!
//! The crate reproduces a coplanar octupole-order secular model of the
//! Planet Nine interaction, not the N-body clone integrations of Köhne &
//! Batygin (2020); the export carries what the crate computes, plus the
//! published orbit of the retrograde Trojan as a labelled reference.

use p9_2020_secular_octupole::libration::{ApsidalMotion, Order, classify_trajectory};
use p9_2020_secular_octupole::octupole::{
    hamiltonian_octupole, hamiltonian_quadrupole, numerical_octupole_amplitude, octupole_coupling,
    octupole_strength_ratio,
};
use p9_2020_secular_octupole::published::{P9_A_AU, P9_E, P9_MASS_EARTH};
use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};
use serde_json::{Value, json};

/// Prototype distant orbit of the crate's headline test (AU).
const A_PROTOTYPE: f64 = 250.0;
/// Initial condition integrated under both truncations.
const E0: f64 = 0.3;
const DVARPI0: f64 = 0.3;

/// 514107 Ka'epaoka'awela (Wiegert et al. 2017), quoted by the paper.
const TROJAN_A_AU: f64 = 5.14;
const TROJAN_E: f64 = 0.38;
const TROJAN_I_DEG: f64 = 163.0;
const JUPITER_A_AU: f64 = 5.203;

/// Eccentricity grid of the level curves.
const N_E: usize = 400;

/// The level curve of the octupole Hamiltonian through (`e0`, `dvarpi0`), as a
/// polyline in (Δϖ [deg], e): the −Δϖ branch walked down in e, then the +Δϖ
/// branch walked back up. Librating curves close on themselves.
fn level_curve(e0: f64, dvarpi0: f64, gm: f64) -> (Vec<f64>, Vec<f64>) {
    let h0 = hamiltonian_octupole(A_PROTOTYPE, e0, dvarpi0, P9_A_AU, P9_E, gm);
    let c_oct = octupole_coupling(A_PROTOTYPE, P9_A_AU, P9_E, gm);
    let branch: Vec<(f64, f64)> = (1..N_E)
        .filter_map(|k| {
            let e = 0.95 * k as f64 / N_E as f64;
            let h_quad = hamiltonian_quadrupole(A_PROTOTYPE, e, P9_A_AU, P9_E, gm);
            let cos_dv = (h0 - h_quad) / (c_oct * e * P9_E);
            (cos_dv.abs() <= 1.0).then(|| (cos_dv.acos().to_degrees(), e))
        })
        .collect();
    let mut dv = Vec::new();
    let mut ecc = Vec::new();
    for &(d, e) in branch.iter().rev() {
        dv.push(-d);
        ecc.push(e);
    }
    for &(d, e) in &branch {
        dv.push(d);
        ecc.push(e);
    }
    (dv, ecc)
}

fn motion_name(m: ApsidalMotion) -> &'static str {
    match m {
        ApsidalMotion::Circulation => "circulation",
        ApsidalMotion::Libration => "libration",
    }
}

pub fn export() -> Value {
    let gm = P9_MASS_EARTH * EARTH_MASS_SOLAR * GM_SUN;
    let epsilon = octupole_strength_ratio(A_PROTOTYPE, P9_A_AU, P9_E);

    // The same orbit under the two truncations.
    let (dt, n_steps) = (2.0e6, 400_000);
    let (quad, _, _) = classify_trajectory(
        Order::Quadrupole,
        A_PROTOTYPE,
        E0,
        DVARPI0,
        P9_A_AU,
        P9_E,
        gm,
        dt,
        n_steps,
    );
    let (oct, dv_min, dv_max) = classify_trajectory(
        Order::Octupole,
        A_PROTOTYPE,
        E0,
        DVARPI0,
        P9_A_AU,
        P9_E,
        gm,
        dt,
        n_steps,
    );

    // Level curves of the octupole Hamiltonian, classified by the integrator.
    let starts: [(f64, f64); 9] = [
        (0.45, 0.0),
        (0.40, 0.0),
        (E0, DVARPI0),
        (0.20, 0.0),
        (0.10, 0.0),
        (0.15, std::f64::consts::PI),
        (0.35, std::f64::consts::PI),
        (0.55, std::f64::consts::PI),
        (0.75, std::f64::consts::PI),
    ];
    let curves: Vec<Value> = starts
        .iter()
        .map(|&(e0, dv0)| {
            let (dv, ecc) = level_curve(e0, dv0, gm);
            let (motion, _, _) = classify_trajectory(
                Order::Octupole,
                A_PROTOTYPE,
                e0,
                dv0,
                P9_A_AU,
                P9_E,
                gm,
                dt,
                200_000,
            );
            json!({
                "e0": e0,
                "dvarpi0_deg": dv0.to_degrees(),
                "motion": motion_name(motion),
                "dvarpi_deg": dv,
                "e": ecc,
            })
        })
        .collect();

    // Octupole strength across the distant belt.
    let a_grid: Vec<f64> = (0..=50).map(|k| 100.0 + 10.0 * k as f64).collect();
    let eps_curve: Vec<f64> = a_grid
        .iter()
        .map(|&a| octupole_strength_ratio(a, P9_A_AU, P9_E))
        .collect();

    // How far the analytic octupole can be trusted: its amplitude against the
    // cos(Δϖ) harmonic of the exact ring average.
    let e_check = 0.2;
    let check: Vec<Value> = [0.03, 0.06, 0.1, 0.15, 0.2, 0.25, 0.3]
        .iter()
        .map(|&alpha| {
            let a = alpha * P9_A_AU;
            let exact = numerical_octupole_amplitude(a, e_check, P9_A_AU, P9_E, gm, 128, 0.0, 48);
            let analytic = octupole_coupling(a, P9_A_AU, P9_E, gm) * e_check * P9_E;
            json!({"alpha": alpha, "a_au": a, "analytic_over_exact": analytic / exact})
        })
        .collect();

    json!({
        "p9": {"mass_earth": P9_MASS_EARTH, "a_au": P9_A_AU, "e": P9_E},
        "a_prototype_au": A_PROTOTYPE,
        "epsilon_oct": epsilon,
        "e_fixed_point": epsilon / 3.0,
        "start": {"e": E0, "dvarpi_deg": DVARPI0.to_degrees()},
        "quadrupole_motion": motion_name(quad),
        "octupole_motion": motion_name(oct),
        "libration_min_deg": dv_min.to_degrees(),
        "libration_max_deg": dv_max.to_degrees(),
        "curves": curves,
        "strength": {"a_au": a_grid, "epsilon_oct": eps_curve},
        "ring_check": check,
        "trojan": {
            "name": "Ka'epaoka'awela",
            "a_au": TROJAN_A_AU,
            "e": TROJAN_E,
            "i_deg": TROJAN_I_DEG,
            "jupiter_a_au": JUPITER_A_AU,
        },
    })
}
