//! Synthetic scattered disk simulations.
//!
//! Reproduces the planar (Section 5.1) and 3D (Section 5.2) scattered disk
//! experiments from Batygin & Brown (2016).
//!
//! Planar: 400 particles, a ∈ (50,550), q ∈ (30,50), coplanar P9
//! 3D: 320 particles, a ∈ (150,550), half-normal i (σ=15°), inclined P9

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::initial_conditions::planets;
use p9_core::initial_conditions::scattered_disk::{
    generate_planar_disk, generate_scattered_disk, ScatteredDiskConfig,
};
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::*;

/// Result snapshot from a scattered disk simulation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DiskSnapshot {
    pub t: f64,
    pub elements: Vec<OrbitalElements>,
    pub active_count: usize,
    pub total_count: usize,
}

/// Configuration for a scattered disk simulation run.
#[derive(Debug, Clone)]
pub struct DiskSimConfig {
    pub p9: P9Params,
    pub n_particles: usize,
    pub a_min: f64,
    pub a_max: f64,
    pub q_min: f64,
    pub q_max: f64,
    pub sigma_i: f64,
    pub t_total: f64,
    pub dt: f64,
    pub snapshot_interval: f64,
    /// Include Neptune as a direct integrator body
    pub include_neptune: bool,
}

impl DiskSimConfig {
    /// Standard planar configuration from Section 5.1.
    pub fn planar_nominal() -> Self {
        Self {
            p9: P9Params::nominal_2016(),
            n_particles: 400,
            a_min: 50.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 0.0,
            t_total: 4.0 * GYR_DAYS,
            dt: JUPITER_PERIOD_DAYS / 20.0, // 1/20 Jupiter's period ≈ 216.6 days
            snapshot_interval: 1e9 * YEAR_DAYS, // every 1 Gyr
            include_neptune: false,
        }
    }

    /// Standard 3D configuration from Section 5.2.
    pub fn inclined_nominal() -> Self {
        let mut p9 = P9Params::nominal_2016();
        p9.i = 30.0 * DEG2RAD;
        p9.omega = 150.0 * DEG2RAD;

        Self {
            p9,
            n_particles: 320,
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 15.0 * DEG2RAD,
            t_total: 4.0 * GYR_DAYS,
            dt: JUPITER_PERIOD_DAYS / 20.0,
            snapshot_interval: 1e9 * YEAR_DAYS,
            include_neptune: true,
        }
    }

    /// Quick test config (10 kyr, few particles).
    pub fn quick_test() -> Self {
        Self {
            p9: P9Params::nominal_2016(),
            n_particles: 10,
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 15.0 * DEG2RAD,
            t_total: 1e4 * YEAR_DAYS,
            dt: 300.0,
            snapshot_interval: 5e3 * YEAR_DAYS,
            include_neptune: false,
        }
    }
}

/// Run a scattered disk simulation with the given configuration.
///
/// Returns a series of snapshots at the specified interval.
pub fn run_scattered_disk(config: &DiskSimConfig, seed: u64) -> Vec<DiskSnapshot> {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    // Generate initial particle population
    let particles_init = if config.sigma_i > 0.0 {
        let disk_config = ScatteredDiskConfig {
            a_min: config.a_min,
            a_max: config.a_max,
            q_min: config.q_min,
            q_max: config.q_max,
            sigma_i: config.sigma_i,
            n_particles: config.n_particles,
        };
        generate_scattered_disk(&disk_config, &mut rng)
    } else {
        generate_planar_disk(
            config.a_min,
            config.a_max,
            config.q_min,
            config.q_max,
            config.n_particles,
            &mut rng,
        )
    };

    let n_actual = particles_init.len();
    let mut particles = particles_init;
    let mut active: Vec<bool> = vec![true; n_actual];

    // Set up massive bodies
    let mut bodies = Vec::new();

    if config.include_neptune {
        let neptune_bodies = planets::neptune_only_j2000();
        bodies.extend(neptune_bodies);
    }

    // Add Planet Nine
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

    let mut snapshots = Vec::new();
    let mut t = 0.0;
    let mut next_snapshot = 0.0;

    // Record initial state
    snapshots.push(take_snapshot(&particles, &active, t, n_actual));

    let n_steps = (config.t_total / config.dt).ceil() as usize;
    for step in 0..n_steps {
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

            // Progress reporting for long runs
            if step % 100_000 == 0 {
                let active_count = active.iter().filter(|&&a| a).count();
                eprintln!(
                    "  t = {:.2} Gyr, {}/{} particles active",
                    t / GYR_DAYS,
                    active_count,
                    n_actual
                );
            }
        }
    }

    // Final snapshot
    snapshots.push(take_snapshot(&particles, &active, t, n_actual));

    snapshots
}

fn take_snapshot(particles: &[StateVector], active: &[bool], t: f64, total: usize) -> DiskSnapshot {
    let mut elements = Vec::new();
    for (i, p) in particles.iter().enumerate() {
        if active[i] {
            elements.push(cartesian_to_elements(p, GM_SUN));
        }
    }
    let active_count = elements.len();
    DiskSnapshot {
        t,
        elements,
        active_count,
        total_count: total,
    }
}

/// Compute clustering statistics for a snapshot.
///
/// Returns (mean_varpi, std_varpi, mean_omega, std_omega, mean_omega_big, std_omega_big)
/// for particles with a > a_min_cluster.
pub fn clustering_statistics(
    snapshot: &DiskSnapshot,
    a_min_cluster: f64,
    varpi_p9: f64,
) -> Option<(f64, f64, f64, f64)> {
    let mut dvarpis = Vec::new();
    let mut omegas = Vec::new();

    for elem in &snapshot.elements {
        if elem.a < a_min_cluster {
            continue;
        }
        let varpi = elem.omega + elem.omega_big;
        let mut dv = varpi - varpi_p9;
        while dv > PI {
            dv -= 2.0 * PI;
        }
        while dv < -PI {
            dv += 2.0 * PI;
        }
        dvarpis.push(dv);
        omegas.push(elem.omega);
    }

    if dvarpis.is_empty() {
        return None;
    }

    let n = dvarpis.len() as f64;
    let mean_dv = dvarpis.iter().sum::<f64>() / n;
    let std_dv = (dvarpis.iter().map(|&x| (x - mean_dv).powi(2)).sum::<f64>() / n).sqrt();

    // Circular mean for ω
    let sin_sum: f64 = omegas.iter().map(|&w| w.sin()).sum();
    let cos_sum: f64 = omegas.iter().map(|&w| w.cos()).sum();
    let mean_omega = sin_sum.atan2(cos_sum);
    let r_bar = (sin_sum * sin_sum + cos_sum * cos_sum).sqrt() / n;
    let std_omega = (-2.0 * r_bar.ln()).sqrt(); // circular standard deviation

    Some((mean_dv, std_dv, mean_omega, std_omega))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_scattered_disk_runs() {
        let config = DiskSimConfig::quick_test();
        let snapshots = run_scattered_disk(&config, 42);

        assert!(
            snapshots.len() >= 2,
            "Should have at least initial and final snapshots"
        );

        let initial = &snapshots[0];
        assert!(
            initial.active_count > 0,
            "Should start with active particles"
        );
        assert_eq!(initial.total_count, config.n_particles);

        // Check final snapshot has valid elements
        let final_snap = snapshots.last().unwrap();
        for elem in &final_snap.elements {
            assert!(elem.a > 0.0, "Semi-major axis should be positive");
            assert!(elem.e >= 0.0, "Eccentricity should be non-negative");
        }
    }

    #[test]
    fn test_clustering_statistics_computation() {
        let snapshot = DiskSnapshot {
            t: 0.0,
            elements: vec![
                OrbitalElements {
                    a: 300.0,
                    e: 0.7,
                    i: 0.0,
                    omega: 5.5,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 400.0,
                    e: 0.8,
                    i: 0.0,
                    omega: 5.6,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 100.0,
                    e: 0.5,
                    i: 0.0,
                    omega: 1.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 3,
            total_count: 3,
        };

        let stats = clustering_statistics(&snapshot, 250.0, 2.5);
        assert!(stats.is_some());

        let (mean_dv, std_dv, _mean_omega, _std_omega) = stats.unwrap();
        // 2 particles with a > 250 should be analyzed
        assert!(mean_dv.is_finite());
        assert!(std_dv.is_finite());
    }
}
