//! Long-term integration of a distant test-particle population under the
//! giant planets, Planet Nine, the Galactic tide and passing stars (paper
//! Section 2.2).
//!
//! Paper set-up: MERCURY6 Bulirsch–Stoer, Neptune and Planet Nine as active
//! bodies, Jupiter–Uranus as a solar J2, Galactic tide, impulse-approximation
//! passing stars, 4 Gyr, statistics over the final Gyr sampled every 10 Myr.
//!
//! Here: p9-core's Wisdom–Holman integrator with Planet Nine as the massive
//! body, all four giants absorbed into the orbit-averaged J2/J4 ring
//! (`j2_jsun_force`; the analysed particles have q > 40 AU, well outside
//! Neptune's scattering region, so the paper's direct Neptune is not needed
//! for the width statistic), the Galactic tide, and the paper's passing
//! stars. The default [`RunConfig::reduced_scale`] integrates for 100 Myr —
//! long enough to show that the width is an invariant of the secular
//! dynamics — while [`RunConfig::paper_scale`] is the 4 Gyr configuration
//! (hours; used behind `#[ignore]`).
//!
//! Particles are massless, so the run is split into independent chunks that
//! evolve in parallel, each with its own copy of Planet Nine (whose motion
//! is identical in every chunk).

use p9_core::constants::{DEG2RAD, GM_SUN, GYR_DAYS, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::scattered_disk_sim::j2_jsun_force;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{
    cartesian_to_elements, MassiveBody, OrbitalElements, P9Params, SimConfig, StateVector,
};
use rand::SeedableRng;
use rayon::prelude::*;

use crate::clustering::perihelion_concentration;
use crate::population::{generate, BirthEnvironment};
use crate::stars::{Encounter, PassingStarModel};
use crate::width::{population_width, SelectionCuts};

/// A Planet Nine configuration of the paper's Table 1.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct P9Config {
    pub mass_earth: f64,
    /// Semi-major axis (AU)
    pub a: f64,
    pub e: f64,
    /// Inclination (deg)
    pub i_deg: f64,
}

impl P9Config {
    /// Perihelion distance (AU).
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// p9-core parameter set. The paper does not state ω, Ω or the phase;
    /// the workspace's usual Batygin–Brown-style values are used.
    pub fn params(&self) -> P9Params {
        P9Params {
            mass_earth: self.mass_earth,
            a: self.a,
            e: self.e,
            i: self.i_deg * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }
}

const fn p9(mass_earth: f64, a: f64, e: f64) -> P9Config {
    P9Config {
        mass_earth,
        a,
        e,
        i_deg: 20.0,
    }
}

/// Table 1, cluster-influenced rows: m₉ ∈ {5, 7.07, 10} M⊕, e₉ ∈ {0.2, 0.35,
/// 0.5}, a₉ chosen to keep q₉ ≈ 250–300 AU.
pub const CLUSTER_INFLUENCED_P9: [P9Config; 9] = [
    p9(5.0, 367.0, 0.2),
    p9(5.0, 420.0, 0.35),
    p9(5.0, 480.0, 0.5),
    p9(7.07, 356.0, 0.2),
    p9(7.07, 433.0, 0.35),
    p9(7.07, 497.0, 0.5),
    p9(10.0, 356.0, 0.2),
    p9(10.0, 433.0, 0.35),
    p9(10.0, 540.0, 0.5),
];

/// Table 1, cluster-free rows: 5 M⊕, a₉ ∈ {400, 500} AU, e₉ ∈ {0.25 … 0.55}.
pub const CLUSTER_FREE_P9: [P9Config; 8] = [
    p9(5.0, 400.0, 0.25),
    p9(5.0, 400.0, 0.35),
    p9(5.0, 400.0, 0.45),
    p9(5.0, 400.0, 0.55),
    p9(5.0, 500.0, 0.25),
    p9(5.0, 500.0, 0.35),
    p9(5.0, 500.0, 0.45),
    p9(5.0, 500.0, 0.55),
];

/// Full description of one integration.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RunConfig {
    pub environment: BirthEnvironment,
    pub p9: P9Config,
    pub n_particles: usize,
    /// Perihelion ceiling of the initial population (AU).
    pub q_max: f64,
    /// Total integration time (days).
    pub t_days: f64,
    /// Integrator step (days).
    pub dt_days: f64,
    /// Orbital-element snapshot cadence (days).
    pub snapshot_every_days: f64,
    /// Snapshots at or after this time enter the pooled statistics (days).
    pub analysis_from_days: f64,
    pub passing_stars: Option<PassingStarModel>,
    pub galactic_tide: bool,
    pub seed: u64,
}

impl RunConfig {
    /// 400 Myr, 10 000-day steps (27 yr, ~37 steps per orbit at a = 100 AU),
    /// snapshots every 10 Myr, statistics pooled over the second half; 192
    /// cluster-influenced particles (the cuts select a smaller fraction of
    /// that population) or 128 cluster-free ones. About 6 CPU-seconds per
    /// particle.
    pub fn reduced_scale(environment: BirthEnvironment, p9: P9Config) -> Self {
        let n_particles = match environment {
            BirthEnvironment::ClusterInfluenced => 192,
            BirthEnvironment::ClusterFree => 128,
        };
        Self {
            environment,
            p9,
            n_particles,
            q_max: environment.default_q_max(),
            t_days: 400.0e6 * YEAR_DAYS,
            dt_days: 10_000.0,
            snapshot_every_days: 10.0e6 * YEAR_DAYS,
            analysis_from_days: 200.0e6 * YEAR_DAYS,
            passing_stars: Some(PassingStarModel::paper()),
            galactic_tide: true,
            seed: 20260717, // arXiv:2607.15646 posting date
        }
    }

    /// The paper's scale: 10⁴ (cluster-influenced) or 10³ (cluster-free)
    /// particles, 4 Gyr, 10 Myr snapshots pooled over the final Gyr.
    pub fn paper_scale(environment: BirthEnvironment, p9: P9Config) -> Self {
        let n_particles = match environment {
            BirthEnvironment::ClusterInfluenced => 10_000,
            BirthEnvironment::ClusterFree => 1_000,
        };
        Self {
            environment,
            p9,
            n_particles,
            q_max: environment.default_q_max(),
            t_days: 4.0 * GYR_DAYS,
            dt_days: 10_000.0,
            snapshot_every_days: 10.0e6 * YEAR_DAYS,
            analysis_from_days: 3.0 * GYR_DAYS,
            passing_stars: Some(PassingStarModel::paper()),
            galactic_tide: true,
            seed: 20260717,
        }
    }
}

/// Orbital elements of the surviving particles at one epoch.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Snapshot {
    pub t_days: f64,
    pub elements: Vec<OrbitalElements>,
}

/// Output of [`run`].
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RunResult {
    pub config: RunConfig,
    pub snapshots: Vec<Snapshot>,
    pub n_initial: usize,
    pub n_final: usize,
}

impl RunResult {
    /// Snapshots inside the analysis window.
    pub fn analysis_snapshots(&self) -> impl Iterator<Item = &Snapshot> {
        self.snapshots
            .iter()
            .filter(|s| s.t_days >= self.config.analysis_from_days)
    }

    /// Pooled orbital elements over the analysis window.
    pub fn pooled_elements(&self) -> Vec<OrbitalElements> {
        self.analysis_snapshots()
            .flat_map(|s| s.elements.iter().copied())
            .collect()
    }

    /// Width `w` (degrees) of the initial population under `cuts`.
    pub fn initial_width_with(&self, cuts: &SelectionCuts) -> f64 {
        population_width(&self.snapshots[0].elements, cuts).0
    }

    /// Width `w` (degrees) pooled over the analysis window under `cuts`.
    pub fn pooled_width(&self, cuts: &SelectionCuts) -> f64 {
        population_width(&self.pooled_elements(), cuts).0
    }

    /// Width `w` (degrees) of the initial population, paper selection.
    pub fn initial_width(&self) -> f64 {
        self.initial_width_with(&SelectionCuts::width_sample())
    }

    /// Width `w` (degrees) pooled over the analysis window (the paper's
    /// Table 1 quantity).
    pub fn final_width(&self) -> f64 {
        self.pooled_width(&SelectionCuts::width_sample())
    }

    /// Width `w` (degrees) at every snapshot under `cuts`, as
    /// `(t_myr, w_deg, n_selected)`, for diagnostics.
    pub fn width_history(&self, cuts: &SelectionCuts) -> Vec<(f64, f64, usize)> {
        self.snapshots
            .iter()
            .map(|s| {
                let (w, n) = population_width(&s.elements, cuts);
                (s.t_days / (1.0e6 * YEAR_DAYS), w, n)
            })
            .collect()
    }

    /// von Mises κ of the longitudes of perihelion pooled over the analysis
    /// window (Table 1).
    pub fn final_kappa(&self) -> f64 {
        perihelion_concentration(&self.pooled_elements())
    }

    /// Surviving fraction.
    pub fn survival(&self) -> f64 {
        self.n_final as f64 / self.n_initial as f64
    }
}

/// Integrate one population.
pub fn run(config: &RunConfig) -> RunResult {
    let mut rng = rand::rngs::StdRng::seed_from_u64(config.seed);
    let particles = generate(
        config.environment,
        config.n_particles,
        config.q_max,
        &mut rng,
    );
    let encounters: Vec<(f64, Encounter)> = config
        .passing_stars
        .map(|m| m.schedule(config.t_days, &mut rng))
        .unwrap_or_default();

    let mut extra = vec![j2_jsun_force()];
    if config.galactic_tide {
        extra.push(ExtraForce::GalacticTide);
    }
    let sim = SimConfig {
        dt: config.dt_days,
        t_start: 0.0,
        t_end: config.t_days,
        // Below ~12 AU the averaged giant-planet ring is meaningless; the
        // paper likewise drops strongly Neptune-crossing orbits.
        removal_inner_au: 12.0,
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

    let chunk = (config.n_particles / (2 * rayon::current_num_threads())).max(4);
    let p9_body = config.p9.params().to_body();
    let chunks: Vec<(Vec<Vec<OrbitalElements>>, usize)> = particles
        .par_chunks(chunk)
        .map(|part| {
            integrate_chunk(
                part,
                &p9_body,
                &extra,
                &sim,
                &encounters,
                n_steps,
                snap_every,
                n_snaps,
            )
        })
        .collect();

    let mut snapshots: Vec<Snapshot> = (0..n_snaps)
        .map(|k| Snapshot {
            t_days: (k * snap_every) as f64 * config.dt_days,
            elements: Vec::new(),
        })
        .collect();
    let mut n_final = 0;
    for (per_snap, survivors) in chunks {
        n_final += survivors;
        for (k, els) in per_snap.into_iter().enumerate() {
            snapshots[k].elements.extend(els);
        }
    }
    RunResult {
        config: config.clone(),
        snapshots,
        n_initial: config.n_particles,
        n_final,
    }
}

#[allow(clippy::too_many_arguments)]
fn integrate_chunk(
    initial: &[StateVector],
    p9_body: &MassiveBody,
    extra: &[ExtraForce],
    sim: &SimConfig,
    encounters: &[(f64, Encounter)],
    n_steps: usize,
    snap_every: usize,
    n_snaps: usize,
) -> (Vec<Vec<OrbitalElements>>, usize) {
    let mut integrator = WhmIntegrator::with_extra_forces(extra.to_vec());
    integrator.parallel = false;
    let mut bodies = vec![p9_body.clone()];
    let mut particles = initial.to_vec();
    let mut active = vec![true; particles.len()];
    let mut per_snap: Vec<Vec<OrbitalElements>> = Vec::with_capacity(n_snaps);
    let mut next_encounter = 0;

    let record = |particles: &[StateVector], active: &[bool]| -> Vec<OrbitalElements> {
        particles
            .iter()
            .zip(active)
            .filter(|(_, &a)| a)
            .map(|(s, _)| cartesian_to_elements(s, GM_SUN))
            .filter(|e| e.e < 1.0 && e.a > 0.0)
            .collect()
    };
    per_snap.push(record(&particles, &active));

    for step in 1..=n_steps {
        integrator.step(&mut bodies, &mut particles, &mut active, sim.dt, sim);
        let t = step as f64 * sim.dt;
        while next_encounter < encounters.len() && encounters[next_encounter].0 <= t {
            let enc = &encounters[next_encounter].1;
            for body in bodies.iter_mut() {
                body.state.vel += enc.impulse(&body.state.pos);
            }
            for (p, &a) in particles.iter_mut().zip(&active) {
                if a {
                    p.vel += enc.impulse(&p.pos);
                }
            }
            next_encounter += 1;
        }
        if step % snap_every == 0 {
            per_snap.push(record(&particles, &active));
        }
    }
    let survivors = active.iter().filter(|&&a| a).count();
    (per_snap, survivors)
}
