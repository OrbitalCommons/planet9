//! Planar N-body simulation for resonance characterization.
//!
//! Sweeps P9 eccentricity from 0 to 0.7 with a₉ = 600 AU, using the
//! J2 approximation for the inner solar system (giant planets replaced
//! by a quadrupole moment on the central body).

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::initial_conditions::scattered_disk::{generate_scattered_disk, ScatteredDiskConfig};
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::*;

use crate::resonance_catalog::{extended_catalog, identify_resonance, resonance_census, Resonance};

/// Configuration for the resonance characterization simulation.
#[derive(Debug, Clone)]
pub struct ResonanceSimConfig {
    /// Planet Nine semi-major axis (AU) — fixed at 600 for this paper
    pub a_p9: f64,
    /// Planet Nine eccentricity (swept from 0 to 0.7)
    pub e_p9: f64,
    /// Planet Nine mass (Earth masses)
    pub mass_earth: f64,
    /// Number of test particles
    pub n_particles: usize,
    /// Semi-major axis range for test particles (AU)
    pub a_min: f64,
    pub a_max: f64,
    /// Perihelion distance range (AU)
    pub q_min: f64,
    pub q_max: f64,
    /// Total integration time (days)
    pub t_total: f64,
    /// Integration timestep (days)
    pub dt: f64,
    /// Snapshot interval (days)
    pub snapshot_interval: f64,
}

impl ResonanceSimConfig {
    /// Nominal config for a given e₉ value.
    pub fn for_eccentricity(e_p9: f64) -> Self {
        Self {
            a_p9: 600.0,
            e_p9,
            mass_earth: 10.0,
            n_particles: 400,
            a_min: 50.0,
            a_max: 1000.0,
            q_min: 30.0,
            q_max: 50.0,
            t_total: 4.0 * GYR_DAYS,
            dt: NEPTUNE_PERIOD_DAYS / 8.0,
            snapshot_interval: 500e6 * YEAR_DAYS,
        }
    }

    /// Quick test config.
    pub fn quick_test(e_p9: f64) -> Self {
        Self {
            n_particles: 20,
            a_min: 150.0,
            a_max: 550.0,
            t_total: 1e4 * YEAR_DAYS,
            snapshot_interval: 5e3 * YEAR_DAYS,
            ..Self::for_eccentricity(e_p9)
        }
    }

    fn p9_params(&self) -> P9Params {
        P9Params {
            mass_earth: self.mass_earth,
            a: self.a_p9,
            e: self.e_p9,
            i: 0.0,
            omega: 0.0,
            omega_big: PI,
            mean_anomaly: 0.0,
        }
    }
}

/// Snapshot with resonance identification.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ResonanceSnapshot {
    pub t: f64,
    pub elements: Vec<OrbitalElements>,
    pub active_count: usize,
    pub total_count: usize,
}

/// Result of the eccentricity sweep.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EccentricitySweepResult {
    /// Eccentricity values tested
    pub eccentricities: Vec<f64>,
    /// For each e₉, the final snapshot
    pub final_snapshots: Vec<ResonanceSnapshot>,
    /// For each e₉, the resonance census
    pub census_results: Vec<Vec<(Resonance, usize)>>,
    /// Planet Nine semi-major axis used
    pub a_p9: f64,
}

/// Run a single planar simulation for a given e₉.
pub fn run_planar_simulation(config: &ResonanceSimConfig, seed: u64) -> ResonanceSnapshot {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let disk_config = ScatteredDiskConfig {
        a_min: config.a_min,
        a_max: config.a_max,
        q_min: config.q_min,
        q_max: config.q_max,
        sigma_i: 0.001, // Near-planar
        n_particles: config.n_particles,
    };
    let mut particles = generate_scattered_disk(&disk_config, &mut rng);
    let n_actual = particles.len();
    let mut active: Vec<bool> = vec![true; n_actual];

    // Sun with J2 (representing giant planets) + Planet Nine
    let p9 = config.p9_params();
    let mut bodies = vec![p9.to_body()];

    let sim_config = SimConfig {
        dt: config.dt,
        t_start: 0.0,
        t_end: config.t_total,
        removal_inner_au: 30.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: config.snapshot_interval,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    let integrator = WhmIntegrator::new();

    let n_steps = (config.t_total / config.dt).ceil() as usize;
    for _ in 0..n_steps {
        integrator.step(
            &mut bodies,
            &mut particles,
            &mut active,
            config.dt,
            &sim_config,
        );
    }

    let mut elements = Vec::new();
    for (i, p) in particles.iter().enumerate() {
        if active[i] {
            elements.push(cartesian_to_elements(p, GM_SUN));
        }
    }
    let active_count = elements.len();

    ResonanceSnapshot {
        t: config.t_total,
        elements,
        active_count,
        total_count: n_actual,
    }
}

/// Run the full eccentricity sweep.
///
/// Sweeps e₉ from 0 to 0.7 in steps of 0.1 (8 simulations).
pub fn eccentricity_sweep(quick: bool) -> EccentricitySweepResult {
    let eccentricities: Vec<f64> = (0..=7).map(|i| i as f64 * 0.1).collect();
    let catalog = extended_catalog();

    let mut final_snapshots = Vec::new();
    let mut census_results = Vec::new();

    for (idx, &e) in eccentricities.iter().enumerate() {
        let config = if quick {
            ResonanceSimConfig::quick_test(e)
        } else {
            ResonanceSimConfig::for_eccentricity(e)
        };

        let snapshot = run_planar_simulation(&config, 42 + idx as u64);
        let census = resonance_census(&snapshot.elements, config.a_p9, &catalog, 0.02);

        final_snapshots.push(snapshot);
        census_results.push(census);
    }

    EccentricitySweepResult {
        eccentricities,
        final_snapshots,
        census_results,
        a_p9: 600.0,
    }
}

/// Classify surviving particles by their resonance type.
pub fn classify_by_resonance_type(elements: &[OrbitalElements], a_p9: f64) -> ResonanceTypeStats {
    let catalog = extended_catalog();
    let mut n_over_1 = 0usize;
    let mut n_over_2 = 0usize;
    let mut high_order = 0usize;
    let mut non_resonant = 0usize;

    for elem in elements {
        match identify_resonance(elem.a, a_p9, &catalog, 0.02) {
            Some((res, _)) => {
                if res.is_n_over_1() {
                    n_over_1 += 1;
                } else if res.is_n_over_2() {
                    n_over_2 += 1;
                } else {
                    high_order += 1;
                }
            }
            None => non_resonant += 1,
        }
    }

    let total_resonant = n_over_1 + n_over_2 + high_order;
    ResonanceTypeStats {
        n_over_1,
        n_over_2,
        high_order,
        non_resonant,
        total: elements.len(),
        p_simple: if total_resonant > 0 {
            (n_over_1 + n_over_2) as f64 / total_resonant as f64
        } else {
            0.0
        },
    }
}

/// Statistics on resonance type distribution.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ResonanceTypeStats {
    pub n_over_1: usize,
    pub n_over_2: usize,
    pub high_order: usize,
    pub non_resonant: usize,
    pub total: usize,
    /// Probability that a resonant particle is in N/1 or N/2.
    pub p_simple: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_config_for_eccentricity() {
        let config = ResonanceSimConfig::for_eccentricity(0.3);
        assert!((config.a_p9 - 600.0).abs() < 1e-10);
        assert!((config.e_p9 - 0.3).abs() < 1e-10);
    }

    #[test]
    fn test_quick_simulation_runs() {
        let config = ResonanceSimConfig::quick_test(0.3);
        let snapshot = run_planar_simulation(&config, 42);
        assert!(snapshot.active_count <= snapshot.total_count);
    }

    #[test]
    fn test_classify_by_resonance_type() {
        let a_p9 = 600.0;
        // 1:2 exterior resonance (p=1 → N/1, period ratio 2)
        let a_12 = Resonance::new(1, 2).semimajor_axis(a_p9);
        // 5:3 resonance (p=5 → high-order)
        let a_53 = Resonance::new(5, 3).semimajor_axis(a_p9);

        let elements = vec![
            OrbitalElements {
                a: a_12,
                e: 0.5,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: a_53,
                e: 0.5,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            // Well away from any resonance
            OrbitalElements {
                a: 50.0,
                e: 0.5,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ];

        let stats = classify_by_resonance_type(&elements, a_p9);
        // Verify classification counts
        assert_eq!(stats.n_over_1, 1, "Expected 1 N/1 (the 1:2 resonance)");
        assert_eq!(
            stats.high_order, 1,
            "Expected 1 high-order (the 5:3 resonance)"
        );
        assert_eq!(stats.non_resonant, 1, "Expected 1 non-resonant (a=100)");
        assert_eq!(stats.total, 3);
    }
}
