//! N-body simulation for narrow vs broad perihelion distribution comparison.
//!
//! Runs two simulations with identical Planet Nine but different initial
//! perihelion distance ranges:
//! - Narrow: q ∈ (30, 36) AU — produces only anti-aligned objects
//! - Broad: q ∈ (30, 300) AU — produces both aligned and anti-aligned
//!
//! Physics configuration (Khain+ 2018 reproduce):
//! - Jupiter+Saturn+Uranus as the orbit-averaged J2 field with monopole
//!   boost (`ExtraForce::J2Jsu`), Neptune and Planet Nine direct;
//! - dt = 3000 d ≈ P_Neptune/20 (the WHM resolution criterion for the
//!   innermost direct body);
//! - the q < 36 AU population grazes Neptune, so steps use the Chambers
//!   hybrid (`p9_core::integrator::hybrid::hybrid_step`) which hands
//!   encounter particles to Bulirsch-Stoer; the J2 extra force is applied
//!   as a symmetric (Strang) half-kick around each hybrid step.

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::forces::{total_extra_acceleration, ExtraForce};
use p9_core::initial_conditions::planets;
use p9_core::initial_conditions::scattered_disk::{generate_scattered_disk, ScatteredDiskConfig};
use p9_core::integrator::hybrid::hybrid_step;
use p9_core::types::*;
use p9_core::units::{days, Time};

/// Configuration for the Kuiper Belt generation simulation.
#[derive(Debug, Clone)]
pub struct KuiperBeltConfig {
    /// Planet Nine parameters
    pub p9: P9Params,
    /// Number of test particles
    pub n_particles: usize,
    /// Semi-major axis range (AU)
    pub a_min: f64,
    pub a_max: f64,
    /// Perihelion distance range (AU)
    pub q_min: f64,
    pub q_max: f64,
    /// Inclination dispersion (rad)
    pub sigma_i: f64,
    /// Total integration time (days)
    pub t_total: f64,
    /// Integration timestep (days)
    pub dt: f64,
    /// Snapshot interval (days)
    pub snapshot_interval: f64,
}

impl KuiperBeltConfig {
    /// Planet Nine for Khain+ 2018: coplanar, a=700, e=0.6, 10 M_Earth.
    fn p9_khain_2018() -> P9Params {
        P9Params {
            mass_earth: 10.0,
            a: 700.0,
            e: 0.6,
            i: 0.0,
            omega: 0.0,
            omega_big: PI, // Anti-aligned with test particles
            mean_anomaly: 0.0,
        }
    }

    /// Narrow perihelion distribution: q ∈ (30, 36) AU.
    pub fn narrow() -> Self {
        Self {
            p9: Self::p9_khain_2018(),
            n_particles: 400,
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 36.0,
            sigma_i: 15.0 * DEG2RAD,
            t_total: 4.0 * GYR_DAYS,
            // J+S+U are J2-averaged; Neptune is the innermost direct body,
            // so dt <= P_N/20 ~ 3000 d.
            dt: 3000.0,
            snapshot_interval: 100e6 * YEAR_DAYS,
        }
    }

    /// Broad perihelion distribution: q ∈ (30, 300) AU.
    pub fn broad() -> Self {
        Self {
            q_max: 300.0,
            ..Self::narrow()
        }
    }

    /// Quick test with narrow distribution.
    pub fn quick_narrow() -> Self {
        Self {
            n_particles: 20,
            t_total: 1e5 * YEAR_DAYS,
            snapshot_interval: 1e4 * YEAR_DAYS,
            ..Self::narrow()
        }
    }

    /// Quick test with broad distribution.
    pub fn quick_broad() -> Self {
        Self {
            n_particles: 20,
            t_total: 1e5 * YEAR_DAYS,
            snapshot_interval: 1e4 * YEAR_DAYS,
            ..Self::broad()
        }
    }
}

/// Snapshot of particle state at a given time. Element slots are indexed by
/// particle identity (None = removed), so per-particle time series — the
/// basis of the Δϖ(t) alignment classification — can be assembled across
/// snapshots.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct KbSnapshot {
    pub t: f64,
    pub elements: Vec<Option<OrbitalElements>>,
    pub active_count: usize,
    pub total_count: usize,
}

impl KbSnapshot {
    /// Snapshot epoch as a typed [`Time`] (see [`Self::t`], in days).
    pub fn time(&self) -> Time {
        days(self.t)
    }
}

/// Result of the full simulation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimulationResult {
    pub snapshots: Vec<KbSnapshot>,
    pub config_summary: String,
    /// Planet Nine's longitude of perihelion for alignment classification
    pub varpi_p9: f64,
}

/// Run the Kuiper Belt generation simulation.
pub fn run_simulation(config: &KuiperBeltConfig, seed: u64) -> SimulationResult {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let disk_config = ScatteredDiskConfig {
        a_min: config.a_min,
        a_max: config.a_max,
        q_min: config.q_min,
        q_max: config.q_max,
        sigma_i: config.sigma_i,
        n_particles: config.n_particles,
    };
    let mut particles = generate_scattered_disk(&disk_config, &mut rng);
    let n_actual = particles.len();
    let mut active: Vec<bool> = vec![true; n_actual];

    // Neptune direct (innermost body resolved by dt = 3000 d) + Planet Nine;
    // Jupiter/Saturn/Uranus enter through the J2-averaged extra force.
    let mut bodies = vec![planets::neptune_j2000()];
    bodies.push(config.p9.to_body());
    let extra_forces = vec![ExtraForce::J2Jsu];

    let sim_config = SimConfig {
        dt: config.dt,
        t_start: 0.0,
        t_end: config.t_total,
        removal_inner_au: 5.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: config.snapshot_interval,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    let varpi_p9 = config.p9.omega_big + config.p9.omega;

    let mut snapshots = Vec::new();
    let mut t = 0.0;
    // First periodic snapshot due at one interval (the t = 0 snapshot is
    // taken below); the previous `t >= next + interval` trigger skipped it.
    let mut next_snapshot = config.snapshot_interval;

    snapshots.push(take_snapshot(&particles, &active, t, n_actual));

    let half_dt = 0.5 * config.dt;
    let n_steps = (config.t_total / config.dt).ceil() as usize;
    for _ in 0..n_steps {
        // Strang composition: J2-averaged giants half-kick, hybrid
        // (WHM + Bulirsch-Stoer encounters) full step, half-kick.
        apply_extra_kick(&extra_forces, &mut bodies, &mut particles, &active, half_dt);
        hybrid_step(
            &mut bodies,
            &mut particles,
            &mut active,
            config.dt,
            &sim_config,
        );
        apply_extra_kick(&extra_forces, &mut bodies, &mut particles, &active, half_dt);
        t += config.dt;

        if t >= next_snapshot {
            snapshots.push(take_snapshot(&particles, &active, t, n_actual));
            next_snapshot += config.snapshot_interval;
        }
    }

    snapshots.push(take_snapshot(&particles, &active, t, n_actual));

    SimulationResult {
        snapshots,
        config_summary: format!(
            "n={}, a=({},{}), q=({},{}), P9: a={}AU e={}",
            config.n_particles,
            config.a_min,
            config.a_max,
            config.q_min,
            config.q_max,
            config.p9.a,
            config.p9.e,
        ),
        varpi_p9,
    }
}

/// Velocity kick from the extra (J2-averaged giant) forces, applied to
/// bodies and active particles for time `dt`.
fn apply_extra_kick(
    forces: &[ExtraForce],
    bodies: &mut [MassiveBody],
    particles: &mut [StateVector],
    active: &[bool],
    dt: f64,
) {
    for body in bodies.iter_mut() {
        body.state.vel += dt * total_extra_acceleration(forces, &body.state.pos);
    }
    for (i, particle) in particles.iter_mut().enumerate() {
        if active[i] {
            particle.vel += dt * total_extra_acceleration(forces, &particle.pos);
        }
    }
}

fn take_snapshot(particles: &[StateVector], active: &[bool], t: f64, total: usize) -> KbSnapshot {
    let elements: Vec<Option<OrbitalElements>> = particles
        .iter()
        .zip(active)
        .map(|(p, &alive)| alive.then(|| cartesian_to_elements(p, GM_SUN)))
        .collect();
    let active_count = elements.iter().flatten().count();
    KbSnapshot {
        t,
        elements,
        active_count,
        total_count: total,
    }
}

/// Compute perihelion distances for surviving particles in a snapshot,
/// indexed by particle (None = removed).
pub fn perihelion_distances(snapshot: &KbSnapshot) -> Vec<Option<f64>> {
    snapshot
        .elements
        .iter()
        .map(|e| e.as_ref().map(|e| e.a * (1.0 - e.e)))
        .collect()
}

/// Wrap an angle to [-π, π].
fn wrap_pi(mut x: f64) -> f64 {
    while x > PI {
        x -= TWO_PI;
    }
    while x < -PI {
        x += TWO_PI;
    }
    x
}

/// Compute delta-varpi (Δϖ = ϖ_particle - ϖ_P9) for all particles in a
/// snapshot, indexed by particle (None = removed).
///
/// Returns values in [-π, π], where:
///   ~0 → aligned with P9
///   ~±π → anti-aligned with P9
pub fn delta_varpi_values(snapshot: &KbSnapshot, varpi_p9: f64) -> Vec<Option<f64>> {
    snapshot
        .elements
        .iter()
        .map(|e| {
            e.as_ref()
                .map(|e| wrap_pi(e.omega_big + e.omega - varpi_p9))
        })
        .collect()
}

/// Assemble the per-particle Δϖ(t) time series across all snapshots:
/// series[i] holds the Δϖ values of particle i at every snapshot where it
/// was still active.
pub fn delta_varpi_series(result: &SimulationResult) -> Vec<Vec<f64>> {
    let n = result
        .snapshots
        .first()
        .map(|s| s.elements.len())
        .unwrap_or(0);
    let mut series = vec![Vec::new(); n];
    for snapshot in &result.snapshots {
        for (i, dv) in delta_varpi_values(snapshot, result.varpi_p9)
            .into_iter()
            .enumerate()
        {
            if let Some(dv) = dv {
                series[i].push(dv);
            }
        }
    }
    series
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_narrow_config() {
        let config = KuiperBeltConfig::narrow();
        assert_eq!(config.n_particles, 400);
        assert!((config.q_max - 36.0).abs() < 1e-10);
        assert!((config.p9.a - 700.0).abs() < 1e-10);
        assert!((config.p9.i).abs() < 1e-10); // Coplanar
                                              // Timestep resolves Neptune, the innermost direct body.
        assert!(config.dt <= NEPTUNE_PERIOD_DAYS / 20.0);
    }

    #[test]
    fn test_broad_config() {
        let config = KuiperBeltConfig::broad();
        assert_eq!(config.n_particles, 400);
        assert!((config.q_max - 300.0).abs() < 1e-10);
    }

    #[test]
    fn test_quick_narrow_runs() {
        let config = KuiperBeltConfig::quick_narrow();
        let result = run_simulation(&config, 42);

        assert!(
            result.snapshots.len() >= 2,
            "Should have at least initial and final snapshots"
        );

        let initial = &result.snapshots[0];
        assert!(initial.active_count > 0);

        // Snapshot trigger fires at every interval (the off-by-one
        // previously skipped the first): 1e5 yr / 1e4 yr = ~10 periodic
        // snapshots + initial + final.
        assert!(
            result.snapshots.len() >= 10,
            "expected ~12 snapshots, got {}",
            result.snapshots.len()
        );

        // Element slots keep particle identity across snapshots.
        for s in &result.snapshots {
            assert_eq!(s.elements.len(), initial.elements.len());
        }
    }

    #[test]
    fn test_quick_broad_runs() {
        let config = KuiperBeltConfig::quick_broad();
        let result = run_simulation(&config, 42);

        assert!(
            result.snapshots.len() >= 2,
            "Should have at least initial and final snapshots"
        );
    }

    #[test]
    fn test_delta_varpi_series_assembled() {
        let config = KuiperBeltConfig::quick_narrow();
        let result = run_simulation(&config, 42);
        let series = delta_varpi_series(&result);
        assert_eq!(series.len(), result.snapshots[0].elements.len());
        // Every surviving particle has one Δϖ entry per snapshot.
        let n_snaps = result.snapshots.len();
        assert!(series.iter().any(|s| s.len() == n_snaps));
        for s in &series {
            for &dv in s {
                assert!((-PI..=PI).contains(&dv));
            }
        }
    }

    #[test]
    fn typed_snapshot_time_matches_days() {
        use uom::si::time::day;
        let snapshot = KbSnapshot {
            t: 1234.5,
            elements: vec![],
            active_count: 0,
            total_count: 0,
        };
        assert!((snapshot.time().get::<day>() - snapshot.t).abs() < 1e-9);
    }

    #[test]
    fn test_perihelion_distances() {
        let snapshot = KbSnapshot {
            t: 0.0,
            elements: vec![
                Some(OrbitalElements {
                    a: 300.0,
                    e: 0.8,
                    i: 0.0,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                }),
                None,
                Some(OrbitalElements {
                    a: 200.0,
                    e: 0.5,
                    i: 0.0,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                }),
            ],
            active_count: 2,
            total_count: 3,
        };

        let qs = perihelion_distances(&snapshot);
        assert!((qs[0].unwrap() - 60.0).abs() < 1e-10); // 300 * 0.2
        assert!(qs[1].is_none());
        assert!((qs[2].unwrap() - 100.0).abs() < 1e-10); // 200 * 0.5
    }

    #[test]
    fn test_delta_varpi_values() {
        let varpi_p9 = PI; // P9 at 180°
        let snapshot = KbSnapshot {
            t: 0.0,
            elements: vec![
                // Aligned with P9 (varpi = π)
                Some(OrbitalElements {
                    a: 300.0,
                    e: 0.8,
                    i: 0.0,
                    omega: PI / 2.0,
                    omega_big: PI / 2.0,
                    mean_anomaly: 0.0,
                }),
                // Anti-aligned (varpi = 0)
                Some(OrbitalElements {
                    a: 200.0,
                    e: 0.5,
                    i: 0.0,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                }),
            ],
            active_count: 2,
            total_count: 2,
        };

        let dvs = delta_varpi_values(&snapshot, varpi_p9);
        assert!(dvs[0].unwrap().abs() < 0.01, "Aligned should have Δϖ ≈ 0");
        assert!(
            (dvs[1].unwrap().abs() - PI).abs() < 0.01,
            "Anti-aligned should have |Δϖ| ≈ π"
        );
    }

    /// Paper-scale regression (Khain+ 2018): the narrow-q run produces an
    /// anti-aligned-dominated surviving population. Full 4 Gyr x 400
    /// particles — run only deliberately.
    #[test]
    #[ignore]
    fn test_paper_scale_narrow_anti_aligned() {
        use crate::population_analysis::{population_statistics_series, Alignment};
        let result = run_simulation(&KuiperBeltConfig::narrow(), 42);
        let stats = population_statistics_series(&result);
        assert!(
            stats.anti_aligned.count >= stats.aligned.count,
            "narrow run should be anti-aligned dominated: {:?}",
            stats
        );
        let _ = Alignment::AntiAligned;
    }
}
