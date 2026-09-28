//! Film export for `p9-2018-kuiper-belt`: the numbers its scene and ledger entry draw.

use p9_2017_dynamics::hamiltonian::SecularHamiltonianParams;
use p9_2018_kuiper_belt::simulation::{
    KuiperBeltConfig, delta_varpi_values, perihelion_distances, run_simulation,
};
use p9_core::constants::YEAR_DAYS;
use serde_json::{Value, json};

use super::p9_2017_dynamics::portrait;

/// Parallel single-particle integrations per primordial disk.
const N_NARROW: usize = 12;
const N_BROAD: usize = 16;
/// Integration span (Myr) and snapshot count: the first 1.5% of the paper's
/// 4 Gyr, all an export can afford.
const T_MYR: f64 = 60.0;
const N_SNAP: f64 = 30.0;

/// Semi-major axis of the secular portrait: the paper's Figure 4. The
/// integrated particles are drawn from a window around it so their tracks
/// can be read against the portrait.
const A_PORTRAIT: f64 = 345.0;
const A_WINDOW: (f64, f64) = (320.0, 370.0);

/// One particle of `base`, integrated for `T_MYR` with seed `seed`; returns
/// its (t, a, q, Δϖ) track while it survives.
fn track(base: &KuiperBeltConfig, seed: u64) -> Value {
    let config = KuiperBeltConfig {
        n_particles: 1,
        a_min: A_WINDOW.0,
        a_max: A_WINDOW.1,
        t_total: T_MYR * 1e6 * YEAR_DAYS,
        snapshot_interval: T_MYR * 1e6 * YEAR_DAYS / N_SNAP,
        ..base.clone()
    };
    let result = run_simulation(&config, seed);
    let mut t = Vec::new();
    let mut a = Vec::new();
    let mut e = Vec::new();
    let mut q = Vec::new();
    let mut dv = Vec::new();
    for snap in &result.snapshots {
        let (Some(el), Some(qi), Some(dvi)) = (
            snap.elements[0].as_ref(),
            perihelion_distances(snap)[0],
            delta_varpi_values(snap, result.varpi_p9)[0],
        ) else {
            break;
        };
        t.push(snap.t / YEAR_DAYS / 1e6);
        a.push(el.a);
        e.push(el.e);
        q.push(qi);
        dv.push(dvi.to_degrees());
    }
    let survived = result.snapshots.last().is_some_and(|s| s.active_count == 1);
    json!({"t_myr": t, "a": a, "e": e, "q": q, "dvarpi_deg": dv, "survived": survived})
}

/// The aligned orbits that can never meet Planet Nine: on the portrait, the
/// secular paths around Δϖ = 0 whose Hamiltonian exceeds every value it takes
/// on an orbit-crossing cell nearby. Returns (that level, the lowest perihelion
/// among those paths in AU).
fn aligned_safe(por: &Value) -> (f64, f64) {
    let a = por["a"].as_f64().expect("portrait a");
    let nums = |v: &Value| -> Vec<f64> {
        v.as_array()
            .expect("array")
            .iter()
            .map(|x| x.as_f64().expect("number"))
            .collect()
    };
    let e = nums(&por["e"]);
    let dv = nums(&por["dvarpi_deg"]);
    let h: Vec<Vec<f64>> = por["h"].as_array().expect("h").iter().map(nums).collect();
    let crossing = &por["crossing"];
    let near_aligned = |j: usize| dv[j].abs() < 90.0;
    let mut level = f64::NEG_INFINITY;
    for (i, row) in h.iter().enumerate() {
        for (j, &hij) in row.iter().enumerate() {
            if near_aligned(j) && crossing[i][j] == json!(true) {
                level = level.max(hij);
            }
        }
    }
    let mut q_min = f64::INFINITY;
    for (i, row) in h.iter().enumerate() {
        for (j, &hij) in row.iter().enumerate() {
            if near_aligned(j) && hij > level {
                q_min = q_min.min(a * (1.0 - e[i]));
            }
        }
    }
    (level, q_min)
}

fn run_disk(base: KuiperBeltConfig, n: usize, seed0: u64) -> Vec<Value> {
    std::thread::scope(|s| {
        let handles: Vec<_> = (0..n)
            .map(|k| {
                let c = base.clone();
                s.spawn(move || track(&c, seed0 + k as u64))
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("integration thread"))
            .collect()
    })
}

pub fn export() -> Value {
    let (narrow, broad) = std::thread::scope(|s| {
        let n = s.spawn(|| run_disk(KuiperBeltConfig::narrow(), N_NARROW, 100));
        let b = s.spawn(|| run_disk(KuiperBeltConfig::broad(), N_BROAD, 500));
        (n.join().expect("narrow"), b.join().expect("broad"))
    });
    let lost = |disk: &[Value]| {
        disk.iter()
            .filter(|p| p["survived"] == json!(false))
            .count() as f64
            / disk.len() as f64
    };
    let cfg = KuiperBeltConfig::broad();
    let por = portrait(A_PORTRAIT, &SecularHamiltonianParams::default_paper());
    let (safe_level, safe_q_min) = aligned_safe(&por);
    json!({
        "a9": cfg.p9.a,
        "e9": cfg.p9.e,
        "m9_earth": cfg.p9.mass_earth,
        "t_myr": T_MYR,
        "a_window": [A_WINDOW.0, A_WINDOW.1],
        "narrow_q": [KuiperBeltConfig::narrow().q_min, KuiperBeltConfig::narrow().q_max],
        "broad_q": [cfg.q_min, cfg.q_max],
        "narrow_lost": lost(&narrow),
        "broad_lost": lost(&broad),
        "narrow": narrow,
        "broad": broad,
        "aligned_safe_level": safe_level,
        "aligned_safe_q_min": safe_q_min,
        "portrait": por,
    })
}
