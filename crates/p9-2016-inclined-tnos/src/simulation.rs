//! Full N-body simulation for inclined TNO generation.
//!
//! Integrates 3,200 test particles with all four giant planets + Planet Nine
//! for 4 Gyr using the hybrid Wisdom-Holman / Bulirsch-Stoer integrator.
//!
//! Key difference from the evidence paper: this uses direct N-body for all
//! giant planets (not J2 averaging), which is critical for capturing
//! Neptune scattering events that decouple particles from Planet Nine.

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::initial_conditions::planets;
use p9_core::initial_conditions::scattered_disk::{generate_scattered_disk, ScatteredDiskConfig};
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::*;

/// Configuration for the inclined TNO simulation.
#[derive(Debug, Clone)]
pub struct InclinedTnoConfig {
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

impl InclinedTnoConfig {
    /// Paper's nominal configuration: 3,200 particles, 4 Gyr.
    pub fn nominal() -> Self {
        let mut p9 = P9Params::nominal_2016();
        p9.a = 600.0;
        p9.e = 0.5;
        p9.i = 30.0 * DEG2RAD;
        p9.omega = 150.0 * DEG2RAD;
        p9.omega_big = PI; // Anti-aligned

        Self {
            p9,
            n_particles: 3200,
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 15.0 * DEG2RAD,
            t_total: 4.0 * GYR_DAYS,
            dt: 300.0,
            snapshot_interval: 100e6 * YEAR_DAYS, // 100 Myr
        }
    }

    /// Quick test configuration.
    pub fn quick_test() -> Self {
        let mut p9 = P9Params::nominal_2016();
        p9.a = 600.0;
        p9.e = 0.5;
        p9.i = 30.0 * DEG2RAD;
        p9.omega = 150.0 * DEG2RAD;
        p9.omega_big = PI;

        Self {
            p9,
            n_particles: 20,
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 15.0 * DEG2RAD,
            t_total: 1e4 * YEAR_DAYS,
            dt: 300.0,
            snapshot_interval: 5e3 * YEAR_DAYS,
        }
    }
}

/// Snapshot of particle state at a given time.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct TnoSnapshot {
    pub t: f64,
    pub elements: Vec<OrbitalElements>,
    pub active_count: usize,
    pub total_count: usize,
}

/// Result of the full simulation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SimulationResult {
    pub snapshots: Vec<TnoSnapshot>,
    pub config_summary: String,
}

/// Run the inclined TNO simulation.
///
/// This is the full simulation with all four giant planets in direct N-body.
pub fn run_simulation(config: &InclinedTnoConfig, seed: u64) -> SimulationResult {
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

    // All four giant planets + Planet Nine
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

    let mut snapshots = Vec::new();
    let mut t = 0.0;
    let mut next_snapshot = 0.0;

    // Initial snapshot
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

    SimulationResult {
        snapshots,
        config_summary: format!(
            "n={}, a=({},{}), q=({},{}), σi={:.0}°, P9: m={}Me a={}AU e={} i={:.0}°",
            config.n_particles,
            config.a_min,
            config.a_max,
            config.q_min,
            config.q_max,
            config.sigma_i * RAD2DEG,
            config.p9.mass_earth,
            config.p9.a,
            config.p9.e,
            config.p9.i * RAD2DEG,
        ),
    }
}

fn take_snapshot(particles: &[StateVector], active: &[bool], t: f64, total: usize) -> TnoSnapshot {
    let mut elements = Vec::new();
    for (i, p) in particles.iter().enumerate() {
        if active[i] {
            elements.push(cartesian_to_elements(p, GM_SUN));
        }
    }
    let active_count = elements.len();
    TnoSnapshot {
        t,
        elements,
        active_count,
        total_count: total,
    }
}

/// Identify highly inclined particles from a snapshot.
///
/// Returns elements with i > i_threshold (default 50°) and a < a_max (default 100 AU).
pub fn find_high_inclination(
    snapshot: &TnoSnapshot,
    i_threshold_deg: f64,
    a_max: f64,
) -> Vec<&OrbitalElements> {
    let i_thresh = i_threshold_deg * DEG2RAD;
    snapshot
        .elements
        .iter()
        .filter(|e| e.i > i_thresh && e.a < a_max)
        .collect()
}

/// Compute the inclination distribution histogram.
///
/// Returns (bin_centers_deg, counts) for binning inclinations in the given range.
pub fn inclination_histogram(
    snapshot: &TnoSnapshot,
    a_min: f64,
    a_max: f64,
    n_bins: usize,
) -> (Vec<f64>, Vec<usize>) {
    let bin_width = 180.0 / n_bins as f64;
    let mut counts = vec![0usize; n_bins];
    let centers: Vec<f64> = (0..n_bins).map(|i| (i as f64 + 0.5) * bin_width).collect();

    for elem in &snapshot.elements {
        if elem.a >= a_min && elem.a <= a_max {
            let i_deg = elem.i * RAD2DEG;
            let bin = (i_deg / bin_width).floor() as usize;
            if bin < n_bins {
                counts[bin] += 1;
            }
        }
    }

    (centers, counts)
}

/// Compute a 2D density map in (a, i) space.
///
/// Returns (a_bins, i_bins, density) where density[i_bin][a_bin] is the count.
pub fn density_map_ai(
    snapshot: &TnoSnapshot,
    a_range: (f64, f64),
    i_range: (f64, f64),
    n_a: usize,
    n_i: usize,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<usize>>) {
    let da = (a_range.1 - a_range.0) / n_a as f64;
    let di = (i_range.1 - i_range.0) / n_i as f64;

    let a_bins: Vec<f64> = (0..n_a)
        .map(|j| a_range.0 + (j as f64 + 0.5) * da)
        .collect();
    let i_bins: Vec<f64> = (0..n_i)
        .map(|j| i_range.0 + (j as f64 + 0.5) * di)
        .collect();

    let mut density = vec![vec![0usize; n_a]; n_i];

    for elem in &snapshot.elements {
        let i_deg = elem.i * RAD2DEG;
        let a_idx = ((elem.a - a_range.0) / da).floor() as isize;
        let i_idx = ((i_deg - i_range.0) / di).floor() as isize;

        if a_idx >= 0 && a_idx < n_a as isize && i_idx >= 0 && i_idx < n_i as isize {
            density[i_idx as usize][a_idx as usize] += 1;
        }
    }

    (a_bins, i_bins, density)
}

/// Compute a 2D density map in (q, i) space.
pub fn density_map_qi(
    snapshot: &TnoSnapshot,
    q_range: (f64, f64),
    i_range: (f64, f64),
    n_q: usize,
    n_i: usize,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<usize>>) {
    let dq = (q_range.1 - q_range.0) / n_q as f64;
    let di = (i_range.1 - i_range.0) / n_i as f64;

    let q_bins: Vec<f64> = (0..n_q)
        .map(|j| q_range.0 + (j as f64 + 0.5) * dq)
        .collect();
    let i_bins: Vec<f64> = (0..n_i)
        .map(|j| i_range.0 + (j as f64 + 0.5) * di)
        .collect();

    let mut density = vec![vec![0usize; n_q]; n_i];

    for elem in &snapshot.elements {
        let q = elem.a * (1.0 - elem.e);
        let i_deg = elem.i * RAD2DEG;
        let q_idx = ((q - q_range.0) / dq).floor() as isize;
        let i_idx = ((i_deg - i_range.0) / di).floor() as isize;

        if q_idx >= 0 && q_idx < n_q as isize && i_idx >= 0 && i_idx < n_i as isize {
            density[i_idx as usize][q_idx as usize] += 1;
        }
    }

    (q_bins, i_bins, density)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quick_simulation_runs() {
        let config = InclinedTnoConfig::quick_test();
        let result = run_simulation(&config, 42);

        assert!(
            result.snapshots.len() >= 2,
            "Should have at least initial and final snapshots"
        );

        let initial = &result.snapshots[0];
        assert!(initial.active_count > 0);
        assert_eq!(initial.total_count, config.n_particles);
    }

    #[test]
    fn test_find_high_inclination() {
        let snapshot = TnoSnapshot {
            t: 0.0,
            elements: vec![
                OrbitalElements {
                    a: 50.0,
                    e: 0.5,
                    i: 110.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 300.0,
                    e: 0.7,
                    i: 120.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 40.0,
                    e: 0.3,
                    i: 10.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 3,
            total_count: 3,
        };

        // i > 50°, a < 100 AU
        let high_i = find_high_inclination(&snapshot, 50.0, 100.0);
        assert_eq!(high_i.len(), 1, "Only one object with i>50° and a<100 AU");
        assert!((high_i[0].a - 50.0).abs() < 0.1);
    }

    #[test]
    fn test_inclination_histogram() {
        let snapshot = TnoSnapshot {
            t: 0.0,
            elements: vec![
                OrbitalElements {
                    a: 50.0,
                    e: 0.3,
                    i: 15.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 60.0,
                    e: 0.4,
                    i: 105.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
                OrbitalElements {
                    a: 70.0,
                    e: 0.5,
                    i: 150.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                },
            ],
            active_count: 3,
            total_count: 3,
        };

        let (centers, counts) = inclination_histogram(&snapshot, 30.0, 100.0, 18);
        assert_eq!(centers.len(), 18);
        assert_eq!(counts.len(), 18);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 3);
    }

    #[test]
    fn test_density_map_construction() {
        let snapshot = TnoSnapshot {
            t: 0.0,
            elements: vec![OrbitalElements {
                a: 50.0,
                e: 0.3,
                i: 100.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            }],
            active_count: 1,
            total_count: 1,
        };

        let (a_bins, i_bins, density) =
            density_map_ai(&snapshot, (0.0, 200.0), (0.0, 180.0), 10, 9);
        assert_eq!(a_bins.len(), 10);
        assert_eq!(i_bins.len(), 9);
        let total: usize = density.iter().flat_map(|row| row.iter()).sum();
        assert_eq!(total, 1);
    }
}
