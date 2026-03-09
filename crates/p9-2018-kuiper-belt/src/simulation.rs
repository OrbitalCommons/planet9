//! N-body simulation for narrow vs broad perihelion distribution comparison.
//!
//! Runs two simulations with identical Planet Nine but different initial
//! perihelion distance ranges:
//! - Narrow: q ∈ (30, 36) AU — produces only anti-aligned objects
//! - Broad: q ∈ (30, 300) AU — produces both aligned and anti-aligned

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::initial_conditions::planets;
use p9_core::initial_conditions::scattered_disk::{generate_scattered_disk, ScatteredDiskConfig};
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::*;

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
            dt: 300.0,
            snapshot_interval: 500e6 * YEAR_DAYS,
        }
    }

    /// Broad perihelion distribution: q ∈ (30, 300) AU.
    pub fn broad() -> Self {
        Self {
            p9: Self::p9_khain_2018(),
            n_particles: 400,
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 300.0,
            sigma_i: 15.0 * DEG2RAD,
            t_total: 4.0 * GYR_DAYS,
            dt: 300.0,
            snapshot_interval: 500e6 * YEAR_DAYS,
        }
    }

    /// Quick test with narrow distribution.
    pub fn quick_narrow() -> Self {
        Self {
            n_particles: 20,
            t_total: 1e4 * YEAR_DAYS,
            snapshot_interval: 5e3 * YEAR_DAYS,
            ..Self::narrow()
        }
    }

    /// Quick test with broad distribution.
    pub fn quick_broad() -> Self {
        Self {
            n_particles: 20,
            t_total: 1e4 * YEAR_DAYS,
            snapshot_interval: 5e3 * YEAR_DAYS,
            ..Self::broad()
        }
    }
}

/// Snapshot of particle state at a given time.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct KbSnapshot {
    pub t: f64,
    pub elements: Vec<OrbitalElements>,
    pub active_count: usize,
    pub total_count: usize,
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

    let mut bodies = planets::giant_planets_j2000();
    bodies.push(config.p9.to_body());

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

    let integrator = WhmIntegrator::new();
    let varpi_p9 = config.p9.omega_big + config.p9.omega;

    let mut snapshots = Vec::new();
    let mut t = 0.0;
    let mut next_snapshot = 0.0;

    snapshots.push(take_snapshot(&particles, &active, t, n_actual));

    let n_steps = (config.t_total / config.dt).ceil() as usize;
    for _ in 0..n_steps {
        integrator.step(
            &mut bodies,
            &mut particles,
            &mut active,
            config.dt,
            &sim_config,
        );
        t += config.dt;

        if t >= next_snapshot + config.snapshot_interval {
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

fn take_snapshot(particles: &[StateVector], active: &[bool], t: f64, total: usize) -> KbSnapshot {
    let mut elements = Vec::new();
    for (i, p) in particles.iter().enumerate() {
        if active[i] {
            elements.push(cartesian_to_elements(p, GM_SUN));
        }
    }
    let active_count = elements.len();
    KbSnapshot {
        t,
        elements,
        active_count,
        total_count: total,
    }
}

/// Compute perihelion distances for all particles in a snapshot.
pub fn perihelion_distances(snapshot: &KbSnapshot) -> Vec<f64> {
    snapshot
        .elements
        .iter()
        .map(|e| e.a * (1.0 - e.e))
        .collect()
}

/// Compute delta-varpi (Δϖ = ϖ_particle - ϖ_P9) for all particles.
///
/// Returns values in [-π, π], where:
///   ~0 → aligned with P9
///   ~±π → anti-aligned with P9
pub fn delta_varpi_values(snapshot: &KbSnapshot, varpi_p9: f64) -> Vec<f64> {
    snapshot
        .elements
        .iter()
        .map(|e| {
            let varpi = e.omega_big + e.omega;
            let mut dv = varpi - varpi_p9;
            while dv > PI {
                dv -= TWO_PI;
            }
            while dv < -PI {
                dv += TWO_PI;
            }
            dv
        })
        .collect()
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
    fn test_perihelion_distances() {
        let snapshot = KbSnapshot {
            t: 0.0,
            elements: vec![
                OrbitalElements {
                    a: 300.0,
                    e: 0.8,
                    i: 0.0,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 200.0,
                    e: 0.5,
                    i: 0.0,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 2,
            total_count: 2,
        };

        let qs = perihelion_distances(&snapshot);
        assert!((qs[0] - 60.0).abs() < 1e-10); // 300 * 0.2
        assert!((qs[1] - 100.0).abs() < 1e-10); // 200 * 0.5
    }

    #[test]
    fn test_delta_varpi_values() {
        let varpi_p9 = PI; // P9 at 180°
        let snapshot = KbSnapshot {
            t: 0.0,
            elements: vec![
                // Aligned with P9 (varpi = π)
                OrbitalElements {
                    a: 300.0,
                    e: 0.8,
                    i: 0.0,
                    omega: PI / 2.0,
                    omega_big: PI / 2.0,
                    mean_anomaly: 0.0,
                },
                // Anti-aligned (varpi = 0)
                OrbitalElements {
                    a: 200.0,
                    e: 0.5,
                    i: 0.0,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 2,
            total_count: 2,
        };

        let dvs = delta_varpi_values(&snapshot, varpi_p9);
        assert!(dvs[0].abs() < 0.01, "Aligned should have Δϖ ≈ 0");
        assert!(
            (dvs[1].abs() - PI).abs() < 0.01,
            "Anti-aligned should have |Δϖ| ≈ π"
        );
    }
}
