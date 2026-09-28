//! Film export for `p9-2016-inclined-tnos`: the numbers its scene and ledger entry draw.
//!
//! The paper's experiment is 3,200 particles for 4 Gyr with all four giant
//! planets integrated directly. The film runs the crate's own initial
//! conditions and Planet Nine at reduced scale: 48 particles for 50 Myr, one
//! per thread, Neptune direct and Jupiter, Saturn and Uranus as an averaged
//! quadrupole. The snapshots are then passed to the crate's own pathway
//! classifier and histogram.

use std::thread;

use p9_2016_inclined_tnos::known_objects::extended_high_i_objects;
use p9_2016_inclined_tnos::kozai_lidov::{Pathway, classify_pathway};
use p9_2016_inclined_tnos::simulation::{
    InclinedTnoConfig, SimulationResult, TnoSnapshot, inclination_histogram, particle_histories,
};
use p9_core::constants::{GM_SUN, RAD2DEG, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets::neptune_j2000;
use p9_core::initial_conditions::scattered_disk::{ScatteredDiskConfig, generate_scattered_disk};
use p9_core::integrator::hybrid::hybrid_step_with_forces;
use p9_core::types::{OrbitalElements, SimConfig, StateVector, cartesian_to_elements};
use rand::SeedableRng;
use serde_json::{Value, json};

/// Reduced-scale run.
const N_PARTICLES: usize = 48;
const T_MYR: f64 = 50.0;
const DT_DAYS: f64 = 3000.0;
const SNAPSHOT_MYR: f64 = 0.5;
const N_THREADS: usize = 48;
const SEED: u64 = 2016;

/// Inclination that counts as "highly inclined" (degrees), as in the crate's
/// pathway classifier.
const HIGH_I_DEG: f64 = 50.0;

/// Tracks exported for the time-series panel.
const N_TRACKS: usize = 8;

/// One chunk of particles integrated with its own copy of the planets; returns
/// each particle's elements at every snapshot until it is removed.
fn integrate(config: &InclinedTnoConfig, start: Vec<StateVector>) -> Vec<Vec<OrbitalElements>> {
    let mut particles = start;
    let mut active = vec![true; particles.len()];
    let mut bodies = vec![neptune_j2000(), config.p9.to_body()];
    let sim = SimConfig {
        dt: DT_DAYS,
        t_start: 0.0,
        t_end: T_MYR * 1e6 * YEAR_DAYS,
        removal_inner_au: 5.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: SNAPSHOT_MYR * 1e6 * YEAR_DAYS,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };
    let extra = [ExtraForce::J2Jsu];
    let n_steps = (sim.t_end / DT_DAYS).ceil() as usize;
    let snap_every = (sim.snapshot_interval_days / DT_DAYS).ceil() as usize;
    let mut out: Vec<Vec<OrbitalElements>> = vec![Vec::new(); particles.len()];
    for step in 0..=n_steps {
        if step > 0 {
            hybrid_step_with_forces(
                &mut bodies,
                &mut particles,
                &mut active,
                DT_DAYS,
                &sim,
                &extra,
            );
        }
        if step % snap_every != 0 {
            continue;
        }
        for (k, p) in particles.iter().enumerate() {
            if !active[k] {
                continue;
            }
            let el = cartesian_to_elements(p, GM_SUN);
            if el.e < 1.0 && el.a > 0.0 && out[k].len() == step / snap_every {
                out[k].push(el);
            }
        }
    }
    out
}

pub fn export() -> Value {
    let config = InclinedTnoConfig::nominal();
    let mut rng = rand::rngs::StdRng::seed_from_u64(SEED);
    let disk = ScatteredDiskConfig {
        a_min: config.a_min,
        a_max: config.a_max,
        q_min: config.q_min,
        q_max: config.q_max,
        sigma_i: config.sigma_i,
        n_particles: N_PARTICLES,
    };
    let start = generate_scattered_disk(&disk, &mut rng);

    let chunk = N_PARTICLES.div_ceil(N_THREADS);
    let series: Vec<Vec<OrbitalElements>> = thread::scope(|s| {
        let handles: Vec<_> = start
            .chunks(chunk)
            .map(|c| {
                let config = &config;
                let c = c.to_vec();
                s.spawn(move || integrate(config, c))
            })
            .collect();
        handles
            .into_iter()
            .flat_map(|h| h.join().unwrap())
            .collect()
    });

    // Rebuild the crate's snapshot structure so its own analysis applies.
    let n_snap = series.iter().map(|s| s.len()).max().unwrap_or(0);
    let snapshots: Vec<TnoSnapshot> = (0..n_snap)
        .map(|k| {
            let mut elements = Vec::new();
            let mut ids = Vec::new();
            for (id, s) in series.iter().enumerate() {
                if k < s.len() {
                    elements.push(s[k]);
                    ids.push(id);
                }
            }
            TnoSnapshot {
                t: k as f64 * SNAPSHOT_MYR * 1e6 * YEAR_DAYS,
                active_count: elements.len(),
                elements,
                ids,
                total_count: N_PARTICLES,
            }
        })
        .collect();
    let result = SimulationResult {
        snapshots,
        config_summary: String::new(),
        seed: SEED,
        bs_failures: 0,
    };

    let histories = particle_histories(&result);
    let max_i = |s: &[OrbitalElements]| s.iter().map(|e| e.i * RAD2DEG).fold(0.0, f64::max);
    let n_high = series.iter().filter(|s| max_i(s) > HIGH_I_DEG).count();
    let n_retro = series.iter().filter(|s| max_i(s) > 90.0).count();
    let n_kozai = histories
        .iter()
        .filter(|h| classify_pathway(h) == Pathway::KozaiLidov)
        .count();
    let n_delivered = histories
        .iter()
        .filter(|h| classify_pathway(h) == Pathway::HighInclination)
        .count();
    let n_survive = series.iter().filter(|s| s.len() == n_snap).count();

    // The most strongly excited particles, as time series.
    let mut order: Vec<usize> = (0..series.len()).collect();
    order.sort_by(|&a, &b| max_i(&series[b]).partial_cmp(&max_i(&series[a])).unwrap());
    let tracks: Vec<Value> = order
        .iter()
        .take(N_TRACKS)
        .map(|&id| {
            let s = &series[id];
            json!({
                "a0_au": s[0].a,
                "t_myr": (0..s.len()).map(|k| k as f64 * SNAPSHOT_MYR).collect::<Vec<_>>(),
                "i_deg": s.iter().map(|e| e.i * RAD2DEG).collect::<Vec<_>>(),
                "q_au": s.iter().map(|e| e.a * (1.0 - e.e)).collect::<Vec<_>>(),
                "a_au": s.iter().map(|e| e.a).collect::<Vec<_>>(),
                "max_i_deg": max_i(s),
            })
        })
        .collect();

    // Where each particle started and the most inclined state it reached.
    let cloud: Vec<Value> = series
        .iter()
        .map(|s| {
            let peak = s
                .iter()
                .max_by(|x, y| x.i.partial_cmp(&y.i).unwrap())
                .unwrap();
            json!({
                "a0_au": s[0].a,
                "i0_deg": s[0].i * RAD2DEG,
                "a_peak_au": peak.a,
                "q_peak_au": peak.a * (1.0 - peak.e),
                "i_peak_deg": peak.i * RAD2DEG,
            })
        })
        .collect();

    let first = &result.snapshots[0];
    let last = result.snapshots.last().unwrap();
    let (centres, start_counts) = inclination_histogram(first, 0.0, 10_000.0, 18);
    let (_, end_counts) = inclination_histogram(last, 0.0, 10_000.0, 18);

    let known: Vec<Value> = extended_high_i_objects()
        .iter()
        .map(|o| {
            json!({
                "name": o.name,
                "a_au": o.elements.a,
                "q_au": o.elements.a * (1.0 - o.elements.e),
                "i_deg": o.elements.i * RAD2DEG,
            })
        })
        .collect();

    json!({
        "n_particles": N_PARTICLES,
        "t_myr": T_MYR,
        "planet": {
            "mass_earth": config.p9.mass_earth,
            "a_au": config.p9.a,
            "e": config.p9.e,
            "i_deg": config.p9.i * RAD2DEG,
        },
        "high_i_deg": HIGH_I_DEG,
        "n_high": n_high,
        "high_fraction": n_high as f64 / N_PARTICLES as f64,
        "n_retrograde": n_retro,
        "n_kozai": n_kozai,
        "n_delivered": n_delivered,
        "n_survive": n_survive,
        "max_i_deg": series.iter().map(|s| max_i(s)).fold(0.0, f64::max),
        "tracks": tracks,
        "cloud": cloud,
        "histogram": {"i_deg": centres, "start": start_counts, "end": end_counts},
        "known": known,
    })
}
