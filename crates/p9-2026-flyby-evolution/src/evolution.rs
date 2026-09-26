//! Long-term evolution of the post-flyby population under the giant planets
//! (paper Section 2, second paragraph).
//!
//! Paper set-up: GENGA hybrid symplectic integrator, the four giant planets,
//! 23-day steps, 4.5 Gyr. Here: p9-core's Chambers hybrid integrator with
//! Neptune as the direct massive body (Bulirsch–Stoer through close
//! encounters) and Jupiter–Uranus as the orbit-averaged J2/J4 ring, 2500-day
//! steps. Tracers reaching 12 AU — where the ring expansion stops being
//! meaningful — are dropped; the paper finds 99% of tracers injected into
//! the planetary region are ejected anyway, so this truncates the "inner"
//! group's survival rather than the outer-population statistics.
//!
//! [`EvolutionConfig::reduced_scale`] follows 512 tracers for 20 Myr, the
//! phase the paper singles out (13% of the 30 < q < 40 AU tracers gone
//! within 10 Myr, the loss rate 100× the present one); the
//! [`EvolutionConfig::paper_scale`] 4.5 Gyr run is behind `#[ignore]`.

use p9_core::constants::{GM_SUN, GYR_DAYS, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets::neptune_j2000;
use p9_core::integrator::hybrid::hybrid_step_with_forces;
use p9_core::types::{
    cartesian_to_elements, elements_to_cartesian, OrbitalElements, SimConfig, StateVector,
};
use rand::seq::SliceRandom;
use rand::SeedableRng;
use rayon::prelude::*;

use crate::groups::{classify, surface_density_weight, Census, DynamicalGroup};

/// Long-term integration set-up.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct EvolutionConfig {
    /// Integration span (days).
    pub t_days: f64,
    /// Step (days); must resolve Neptune (≤ 3000 d).
    pub dt_days: f64,
    /// Number of post-flyby tracers followed (random subset), or all.
    pub n_particles: Option<usize>,
    /// Snapshot cadence (days).
    pub snapshot_every_days: f64,
    /// Tracers inside this heliocentric distance are dropped (AU).
    pub removal_inner_au: f64,
    pub seed: u64,
}

impl EvolutionConfig {
    /// 512 tracers, 20 Myr, snapshots every 2 Myr.
    pub fn reduced_scale() -> Self {
        Self {
            t_days: 20.0e6 * YEAR_DAYS,
            dt_days: 2500.0,
            n_particles: Some(512),
            snapshot_every_days: 2.0e6 * YEAR_DAYS,
            removal_inner_au: 12.0,
            seed: 20260903,
        }
    }

    /// Every bound tracer for 4.56 Gyr, snapshots every 100 Myr (hours).
    pub fn paper_scale() -> Self {
        Self {
            t_days: 4.56 * GYR_DAYS,
            dt_days: 2500.0,
            n_particles: None,
            snapshot_every_days: 100.0e6 * YEAR_DAYS,
            removal_inner_au: 12.0,
            seed: 20260903,
        }
    }
}

/// Elements of every followed tracer at one epoch (`None` once lost).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EvolutionSnapshot {
    pub t_days: f64,
    pub elements: Vec<Option<OrbitalElements>>,
}

impl EvolutionSnapshot {
    pub fn t_myr(&self) -> f64 {
        self.t_days / (1.0e6 * YEAR_DAYS)
    }

    pub fn bound(&self) -> impl Iterator<Item = &OrbitalElements> {
        self.elements.iter().flatten()
    }
}

/// Output of [`evolve`].
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EvolutionResult {
    pub config: EvolutionConfig,
    /// Pre-flyby radius of each followed tracer (AU).
    pub initial_radii: Vec<f64>,
    pub snapshots: Vec<EvolutionSnapshot>,
}

impl EvolutionResult {
    pub fn initial(&self) -> &EvolutionSnapshot {
        &self.snapshots[0]
    }

    pub fn final_snapshot(&self) -> &EvolutionSnapshot {
        self.snapshots.last().unwrap()
    }

    /// Snapshot nearest to `t_myr`.
    pub fn snapshot_at(&self, t_myr: f64) -> &EvolutionSnapshot {
        self.snapshots
            .iter()
            .min_by(|a, b| {
                (a.t_myr() - t_myr)
                    .abs()
                    .partial_cmp(&(b.t_myr() - t_myr).abs())
                    .unwrap()
            })
            .unwrap()
    }

    /// Weighted Table 1 census of the survivors at `snap`.
    pub fn census(&self, snap: &EvolutionSnapshot) -> Census {
        Census::of(
            self.initial_radii
                .iter()
                .zip(&snap.elements)
                .filter_map(|(&r, e)| e.as_ref().map(|e| (r, e))),
        )
    }

    /// Table 1 "N_x / N_init": weight of `group` at `snap` over its weight
    /// at the start (inflow included).
    pub fn retention(&self, snap: &EvolutionSnapshot, group: DynamicalGroup) -> f64 {
        let w0 = self.census(self.initial()).weight(group);
        if w0 <= 0.0 {
            f64::NAN
        } else {
            self.census(snap).weight(group) / w0
        }
    }

    /// Fraction of the tracers initially satisfying `region` that are, at
    /// `snap`, lost or no longer inside it (the paper's "loss from a
    /// region").
    pub fn region_loss(
        &self,
        snap: &EvolutionSnapshot,
        region: impl Fn(&OrbitalElements) -> bool,
    ) -> f64 {
        let mut w0 = 0.0;
        let mut still = 0.0;
        for ((e0, e1), &r0) in self
            .initial()
            .elements
            .iter()
            .zip(&snap.elements)
            .zip(&self.initial_radii)
        {
            if let Some(e0) = e0 {
                if region(e0) {
                    let w = surface_density_weight(r0);
                    w0 += w;
                    if e1.map(|e| region(&e)).unwrap_or(false) {
                        still += w;
                    }
                }
            }
        }
        if w0 <= 0.0 {
            f64::NAN
        } else {
            1.0 - still / w0
        }
    }

    /// Fraction of the initially-`region` tracers that are gone (unbound or
    /// inside the removal radius) at `snap`.
    pub fn ejected_fraction(
        &self,
        snap: &EvolutionSnapshot,
        region: impl Fn(&OrbitalElements) -> bool,
    ) -> f64 {
        let mut w0 = 0.0;
        let mut gone = 0.0;
        for ((e0, e1), &r0) in self
            .initial()
            .elements
            .iter()
            .zip(&snap.elements)
            .zip(&self.initial_radii)
        {
            if let Some(e0) = e0 {
                if region(e0) {
                    let w = surface_density_weight(r0);
                    w0 += w;
                    if e1.is_none() {
                        gone += w;
                    }
                }
            }
        }
        if w0 <= 0.0 {
            f64::NAN
        } else {
            gone / w0
        }
    }

    /// Survivors at `snap` paired with their pre-flyby radii.
    pub fn survivors(&self, snap: &EvolutionSnapshot) -> Vec<(f64, OrbitalElements)> {
        self.initial_radii
            .iter()
            .zip(&snap.elements)
            .filter_map(|(&r, e)| e.map(|e| (r, e)))
            .collect()
    }
}

/// Evolve a post-flyby population. `bound` pairs pre-flyby radius with
/// post-flyby heliocentric elements (see `FlybyOutcome::bound`).
pub fn evolve(bound: &[(f64, OrbitalElements)], config: &EvolutionConfig) -> EvolutionResult {
    let mut rng = rand::rngs::StdRng::seed_from_u64(config.seed);
    let mut chosen: Vec<(f64, OrbitalElements)> = bound.to_vec();
    if let Some(n) = config.n_particles {
        chosen.shuffle(&mut rng);
        chosen.truncate(n);
    }
    let initial_radii: Vec<f64> = chosen.iter().map(|c| c.0).collect();
    let states: Vec<StateVector> = chosen
        .iter()
        .map(|c| elements_to_cartesian(&c.1, GM_SUN))
        .collect();

    let sim = SimConfig {
        dt: config.dt_days,
        t_start: 0.0,
        t_end: config.t_days,
        removal_inner_au: config.removal_inner_au,
        removal_outer_au: 100_000.0,
        snapshot_interval_days: config.snapshot_every_days,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };
    let n_steps = (config.t_days / config.dt_days).ceil() as usize;
    let snap_every = (config.snapshot_every_days / config.dt_days)
        .round()
        .max(1.0) as usize;
    let n_snaps = n_steps / snap_every + 1;
    let chunk = (states.len() / (2 * rayon::current_num_threads())).max(4);

    let per_chunk: Vec<Vec<Vec<Option<OrbitalElements>>>> = states
        .par_chunks(chunk)
        .map(|part| integrate_chunk(part, &sim, n_steps, snap_every, n_snaps))
        .collect();

    let mut snapshots: Vec<EvolutionSnapshot> = (0..n_snaps)
        .map(|k| EvolutionSnapshot {
            t_days: (k * snap_every) as f64 * config.dt_days,
            elements: Vec::with_capacity(states.len()),
        })
        .collect();
    for chunk_snaps in per_chunk {
        for (k, els) in chunk_snaps.into_iter().enumerate() {
            snapshots[k].elements.extend(els);
        }
    }
    EvolutionResult {
        config: *config,
        initial_radii,
        snapshots,
    }
}

fn integrate_chunk(
    initial: &[StateVector],
    sim: &SimConfig,
    n_steps: usize,
    snap_every: usize,
    n_snaps: usize,
) -> Vec<Vec<Option<OrbitalElements>>> {
    let mut bodies = vec![neptune_j2000()];
    let extra = [ExtraForce::J2Jsu];
    let mut particles = initial.to_vec();
    let mut active = vec![true; particles.len()];
    let record = |particles: &[StateVector], active: &mut [bool]| -> Vec<Option<OrbitalElements>> {
        particles
            .iter()
            .zip(active.iter_mut())
            .map(|(s, a)| {
                if !*a {
                    return None;
                }
                let e = cartesian_to_elements(s, GM_SUN);
                if e.e >= 1.0 || e.a <= 0.0 {
                    // Unbound: drop it for good.
                    *a = false;
                    return None;
                }
                Some(e)
            })
            .collect()
    };
    let mut out = Vec::with_capacity(n_snaps);
    out.push(record(&particles, &mut active));
    for step in 1..=n_steps {
        hybrid_step_with_forces(
            &mut bodies,
            &mut particles,
            &mut active,
            sim.dt,
            sim,
            &extra,
        );
        if step % snap_every == 0 {
            out.push(record(&particles, &mut active));
        }
    }
    out
}

/// The 30 < q < 40 AU region the paper tracks most closely.
pub fn near_neptune(el: &OrbitalElements) -> bool {
    let q = el.a * (1.0 - el.e);
    (30.0..40.0).contains(&q)
}

/// Group membership predicate.
pub fn in_group(group: DynamicalGroup) -> impl Fn(&OrbitalElements) -> bool {
    move |el| classify(el) == Some(group)
}
