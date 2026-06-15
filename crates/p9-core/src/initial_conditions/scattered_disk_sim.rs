//! Synthetic scattered disk simulations.
//!
//! Reproduces the planar (Section 5.1) and 3D (Section 5.2) scattered disk
//! experiments from Batygin & Brown (2016).
//!
//! Planar: 400 particles, a ∈ (50,550), q ∈ (30,50), coplanar P9
//! 3D: 320 particles, a ∈ (150,550), half-normal i (σ=15°), inclined P9
//!
//! The giant planets are included as the paper did: orbit-averaged into an
//! enhanced solar quadrupole (J+S+U+N for the planar runs; J+S+U with
//! Neptune direct for the 3D runs). This term sets the KBO apsidal
//! precession — omitting it (as a previous version did, integrating bare
//! Sun + P9) shifts precession rates, libration centers and survival
//! fractions.
//!
//! This is the simulation orchestrator: it builds the massive bodies,
//! integrator and averaged giant-planet field, then drives the test-particle
//! population produced by the generators in [`super::scattered_disk`]
//! (`generate_scattered_disk` / `generate_planar_disk`). It does not
//! duplicate those generators — it wraps them.

use std::sync::Arc;

use crate::analysis::circular::{circular_mean, circular_std, wrap_to_pi};
use crate::constants::*;
use crate::forces::{j2_secular, ExtraForce};
use crate::initial_conditions::planets;
use crate::initial_conditions::scattered_disk::{
    generate_planar_disk, generate_scattered_disk, ScatteredDiskConfig,
};
use crate::integrator::whm::WhmIntegrator;
use crate::types::*;

/// Orbit-averaged J2/J4 field absorbing all four giants (J+S+U from p9-core
/// plus Neptune folded in the same way), with the matching monopole boost.
/// Used when no planet is integrated directly (planar runs). Valid for
/// perihelia q ≳ 30 AU.
pub fn j2_jsun_force() -> ExtraForce {
    let (mut j2, mut j4, mut gm_boost) = j2_secular::combined_j2_jsu();
    j2 += j2_secular::effective_j2(MASS_NEPTUNE_SOLAR, A_NEPTUNE_AU, 1.0);
    j4 += -0.375 * MASS_NEPTUNE_SOLAR * A_NEPTUNE_AU.powi(4);
    gm_boost += GM_NEPTUNE;
    ExtraForce::Custom(Arc::new(move |pos| {
        j2_secular::secular_j2_acceleration(pos, GM_SUN, j2, j4, gm_boost, 1.0)
    }))
}

/// Result snapshot from a scattered disk simulation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DiskSnapshot {
    pub t: f64,
    pub elements: Vec<OrbitalElements>,
    pub active_count: usize,
    pub total_count: usize,
    /// Planet Nine's longitude of perihelion at this snapshot (it precesses
    /// under the averaged giant-planet quadrupole), the reference direction
    /// for Δϖ clustering statistics.
    pub varpi_p9: f64,
    /// RNG seed of the run (recorded for reproducibility).
    pub seed: u64,
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
            // All four giants are absorbed into the J2-averaged field, so
            // the timestep only needs to resolve the q ≈ 30 AU perihelion
            // passages of the particles (~P(30 AU)/60).
            dt: 1000.0,
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
            // Neptune direct (innermost massive body): dt ≲ P_N/20 ≈ 3000 d;
            // 1000 d also resolves the particles' perihelion passages.
            dt: 1000.0,
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

    // Set up massive bodies and the orbit-averaged giant-planet field:
    // Neptune direct + J2-averaged J+S+U, or all four giants averaged when
    // Neptune is not integrated directly (planar runs).
    let mut bodies = Vec::new();

    let extra_force = if config.include_neptune {
        let neptune_bodies = planets::neptune_only_j2000();
        bodies.extend(neptune_bodies);
        ExtraForce::J2Jsu
    } else {
        j2_jsun_force()
    };

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

    let integrator = WhmIntegrator::with_extra_forces(vec![extra_force]);

    let mut snapshots = Vec::new();
    let mut t = 0.0;
    // Snapshot trigger convention (shared with phase_portrait and the
    // inclined-TNO sim): record at t = 0, then whenever t crosses the next
    // multiple of the snapshot interval.
    let mut next_snapshot = config.snapshot_interval;

    // Record initial state
    snapshots.push(take_snapshot(
        &particles, &active, &bodies, t, n_actual, seed,
    ));

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

        if t >= next_snapshot {
            snapshots.push(take_snapshot(
                &particles, &active, &bodies, t, n_actual, seed,
            ));
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

    // Final snapshot (unless the loop already recorded this epoch)
    if snapshots.last().map(|s| s.t) != Some(t) {
        snapshots.push(take_snapshot(
            &particles, &active, &bodies, t, n_actual, seed,
        ));
    }

    snapshots
}

fn take_snapshot(
    particles: &[StateVector],
    active: &[bool],
    bodies: &[MassiveBody],
    t: f64,
    total: usize,
    seed: u64,
) -> DiskSnapshot {
    let mut elements = Vec::new();
    for (i, p) in particles.iter().enumerate() {
        if active[i] {
            elements.push(cartesian_to_elements(p, GM_SUN));
        }
    }
    let active_count = elements.len();
    // P9 is always the last massive body.
    let p9_elem = cartesian_to_elements(&bodies.last().unwrap().state, GM_SUN);
    DiskSnapshot {
        t,
        elements,
        active_count,
        total_count: total,
        varpi_p9: (p9_elem.omega + p9_elem.omega_big).rem_euclid(TWO_PI),
        seed,
    }
}

/// Compute clustering statistics for a snapshot, using circular statistics
/// throughout (Δϖ relative to P9 is an angle, frequently bimodal, so linear
/// means/stds are not meaningful for it).
///
/// Returns (mean_dvarpi, circ_std_dvarpi, mean_omega, circ_std_omega)
/// for particles with a > a_min_cluster. The mean Δϖ is wrapped to (−π, π].
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
        dvarpis.push(wrap_to_pi(varpi - varpi_p9));
        omegas.push(elem.omega);
    }

    if dvarpis.is_empty() {
        return None;
    }

    let mean_dv = wrap_to_pi(circular_mean(&dvarpis)?);
    let std_dv = circular_std(&dvarpis);
    let mean_omega = circular_mean(&omegas)?;
    let std_omega = circular_std(&omegas);

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
            varpi_p9: 2.5,
            seed: 0,
        };

        let stats = clustering_statistics(&snapshot, 250.0, 2.5);
        assert!(stats.is_some());

        let (mean_dv, std_dv, _mean_omega, _std_omega) = stats.unwrap();
        // 2 particles with a > 250 should be analyzed
        assert!(mean_dv.is_finite());
        assert!(std_dv.is_finite());
    }
}
