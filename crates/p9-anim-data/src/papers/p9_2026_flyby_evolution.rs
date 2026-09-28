//! Film export for `p9-2026-flyby-evolution`: the numbers its scene and ledger entry draw.

use p9_2026_flyby_evolution::colours::{E_SPLIT, I_SPLIT_DEG, colour_dichotomy};
use p9_2026_flyby_evolution::evolution::{EvolutionConfig, evolve, near_neptune};
use p9_2026_flyby_evolution::flyby::{FlybyConfig, run_flyby};
use p9_2026_flyby_evolution::groups::{GROUPS, classify};
use p9_2026_flyby_evolution::reference::{
    CAPTURED_BY_PERTURBER, HIGH_E_LOSS, LOW_E_LOSS, NEAR_NEPTUNE_LOSS_4_5GYR,
    NEAR_NEPTUNE_LOSS_10MYR, TABLE1_FRACTION_T0, TABLE1_FRACTION_T4_5GYR, TABLE1_RETENTION_T4_5GYR,
};
use p9_core::constants::{GM_SUN, RAD2DEG, TWO_PI, YEAR_DAYS};
use p9_core::types::{OrbitalElements, elements_to_cartesian};
use serde_json::{Value, json};

/// Tracers sent through the encounter.
const N_FLYBY: usize = 1500;
/// Near-Neptune tracers followed afterwards, and for how long. The paper
/// quotes the loss within 10 Myr; the film follows half of that to keep the
/// export to tens of CPU-seconds.
const N_EVOLVED: usize = 96;
const EVOLVE_MYR: f64 = 5.0;
const SNAPSHOT_MYR: f64 = 0.25;
/// Extent of the drawn perturber path: Sun-star distance at its ends (AU).
const PATH_EXTENT_AU: f64 = 480.0;
const N_PATH: usize = 121;
/// Post-flyby orbits drawn, points per orbit, and the largest drawn orbit (AU).
const N_DRAWN: usize = 90;
const N_ORBIT_POINTS: usize = 64;
const DRAWN_A_MAX_AU: f64 = 260.0;

fn perihelion(el: &OrbitalElements) -> f64 {
    el.a * (1.0 - el.e)
}

/// The orbit as a closed polyline (AU, disc frame), sampled evenly in
/// eccentric anomaly.
fn orbit_polyline(el: &OrbitalElements) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut out = (Vec::new(), Vec::new(), Vec::new());
    for k in 0..=N_ORBIT_POINTS {
        let ea = TWO_PI * k as f64 / N_ORBIT_POINTS as f64;
        let at = OrbitalElements {
            mean_anomaly: ea - el.e * ea.sin(),
            ..*el
        };
        let pos = elements_to_cartesian(&at, GM_SUN).pos;
        out.0.push(pos.x);
        out.1.push(pos.y);
        out.2.push(pos.z);
    }
    out
}

fn tracer_json(r0: f64, el: &OrbitalElements) -> Value {
    json!({
        "r0_au": r0,
        "a_au": el.a,
        "e": el.e,
        "i_deg": el.i * RAD2DEG,
        "q_au": perihelion(el),
        "group": classify(el).map(|g| g.label()),
    })
}

pub fn export() -> Value {
    let config = FlybyConfig::model_a1(N_FLYBY);
    let flyby = config.flyby;
    let out = run_flyby(&config);
    let bound = out.bound();

    // The perturber's parabola through the disc.
    let t_far = flyby.time_at_distance(PATH_EXTENT_AU);
    let path: Vec<Value> = (0..N_PATH)
        .map(|k| {
            let t = t_far * (2.0 * k as f64 / (N_PATH - 1) as f64 - 1.0);
            let s = flyby.relative_state(t);
            json!({"t_yr": t / YEAR_DAYS, "x": s.pos.x, "y": s.pos.y, "z": s.pos.z})
        })
        .collect();
    let peri = flyby.relative_state(0.0);

    let tracers: Vec<Value> = bound.iter().map(|(r0, el)| tracer_json(*r0, el)).collect();
    let unbound_radii: Vec<f64> = out
        .initial_radii
        .iter()
        .zip(&out.elements)
        .filter(|(_, e)| e.is_none())
        .map(|(&r, _)| r)
        .collect();

    let drawable: Vec<&(f64, OrbitalElements)> = bound
        .iter()
        .filter(|(_, el)| el.a < DRAWN_A_MAX_AU)
        .collect();
    let stride = (drawable.len() / N_DRAWN).max(1);
    let drawn: Vec<Value> = drawable
        .iter()
        .step_by(stride)
        .map(|(r0, el)| {
            let (x, y, z) = orbit_polyline(el);
            json!({
                "r0_au": r0,
                "e": el.e,
                "i_deg": el.i * RAD2DEG,
                "group": classify(el).map(|g| g.label()),
                "x": x,
                "y": y,
                "z": z,
            })
        })
        .collect();

    let census = out.census();
    let fractions = census.fractions();
    let families: Vec<Value> = GROUPS
        .iter()
        .enumerate()
        .map(|(k, g)| {
            json!({
                "label": g.label(),
                "fraction": fractions[k],
                "published_t0": TABLE1_FRACTION_T0[k],
                "published_4_5gyr": TABLE1_FRACTION_T4_5GYR[k],
                "published_retention": TABLE1_RETENTION_T4_5GYR[k],
            })
        })
        .collect();
    let deviation: f64 = fractions
        .iter()
        .zip(&TABLE1_FRACTION_T0)
        .map(|(f, p)| (f - p).abs())
        .sum();

    // What Neptune does to the tracers the flyby left near it.
    let near: Vec<(f64, OrbitalElements)> = bound
        .iter()
        .filter(|(_, el)| near_neptune(el))
        .copied()
        .collect();
    let evolution = evolve(
        &near,
        &EvolutionConfig {
            t_days: EVOLVE_MYR * 1e6 * YEAR_DAYS,
            n_particles: Some(N_EVOLVED),
            snapshot_every_days: SNAPSHOT_MYR * 1e6 * YEAR_DAYS,
            ..EvolutionConfig::reduced_scale()
        },
    );
    let last = evolution.final_snapshot();
    let t_myr: Vec<f64> = evolution.snapshots.iter().map(|s| s.t_myr()).collect();
    let loss: Vec<f64> = evolution
        .snapshots
        .iter()
        .map(|s| evolution.region_loss(s, near_neptune))
        .collect();
    let removed: Vec<f64> = evolution
        .snapshots
        .iter()
        .map(|s| evolution.ejected_fraction(s, near_neptune))
        .collect();
    let followed: Vec<Value> = evolution
        .initial()
        .elements
        .iter()
        .zip(&last.elements)
        .zip(&evolution.initial_radii)
        .filter_map(|((e0, e1), &r0)| {
            e0.as_ref().map(|e0| {
                json!({
                    "r0_au": r0,
                    "start": {"a_au": e0.a, "e": e0.e, "i_deg": e0.i * RAD2DEG, "q_au": perihelion(e0)},
                    "end": e1.as_ref().map(|e1| json!({
                        "a_au": e1.a, "e": e1.e, "i_deg": e1.i * RAD2DEG, "q_au": perihelion(e1),
                    })),
                    "still_near_neptune": e1.as_ref().is_some_and(near_neptune),
                })
            })
        })
        .collect();

    let (red_low_i, red_high_i, red_low_e, red_high_e) =
        colour_dichotomy(bound.iter().map(|(r, e)| (*r, e)));

    json!({
        "flyby": {
            "mass_solar": flyby.mass_solar,
            "q_au": flyby.q_au,
            "inclination_deg": flyby.inclination_deg,
            "arg_periastron_deg": flyby.arg_periastron_deg,
            "periastron": {"x": peri.pos.x, "y": peri.pos.y, "z": peri.pos.z},
            "periastron_speed_kms": peri.vel.norm() / p9_core::constants::KMS_TO_AUDAY,
            "path": path,
        },
        "disc": {"r_min_au": config.disc.r_min_au, "r_max_au": config.disc.r_max_au},
        "n_tracers": N_FLYBY,
        "n_bound": out.n_bound(),
        "unbound_fraction": out.unbound_fraction(),
        "published_captured": CAPTURED_BY_PERTURBER,
        "tracers": tracers,
        "unbound_radii_au": unbound_radii,
        "drawn_orbits": drawn,
        "families": families,
        "census_n": census.n,
        "census_deviation": deviation,
        "belt_fraction": fractions[0] + fractions[1],
        "published_belt_fraction": TABLE1_FRACTION_T0[0] + TABLE1_FRACTION_T0[1],
        "sedna_like_fraction": fractions[3],
        "evolution": {
            "n_followed": followed.len(),
            "span_myr": EVOLVE_MYR,
            "t_myr": t_myr,
            "loss": loss,
            "removed": removed,
            "tracers": followed,
            "loss_low_e": evolution.region_loss(last, |e| near_neptune(e) && e.e < 0.2),
            "loss_high_e": evolution.region_loss(last, |e| near_neptune(e) && e.e > 0.4),
            "published_loss_low_e_4_5gyr": LOW_E_LOSS,
            "published_loss_high_e_4_5gyr": HIGH_E_LOSS,
            "published_loss_10myr": NEAR_NEPTUNE_LOSS_10MYR,
            "published_loss_4_5gyr": NEAR_NEPTUNE_LOSS_4_5GYR,
        },
        "near_neptune_loss": evolution.region_loss(last, near_neptune),
        "near_neptune_removed": evolution.ejected_fraction(last, near_neptune),
        "colours": {
            "i_split_deg": I_SPLIT_DEG,
            "e_split": E_SPLIT,
            "very_red_low_i": red_low_i,
            "very_red_high_i": red_high_i,
            "very_red_low_e": red_low_e,
            "very_red_high_e": red_high_e,
        },
    })
}
