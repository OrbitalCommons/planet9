//! Parameter grid and survey runner for the planar survey.
//!
//! Brown & Batygin (2016) systematically explore:
//!   a₉ ∈ (200, 2000) AU in 100 AU steps → 19 values
//!   e₉ ∈ (0.1, 0.9) in 0.1 steps → 9 values
//!   m₉ ∈ {0.1, 1, 10, 20, 30} M_Earth → 5 values
//!
//! Total: 19 × 9 × 5 = 855 parameter combinations
//! (paper says ~320 per mass × 5 masses ≈ 1600; they may skip some unphysical)
//!
//! Each grid point is evaluated by running the planar scattered-disk
//! simulation from p9-2016-evidence (giant planets absorbed into the
//! J2-averaged solar quadrupole, P9 integrated directly) and comparing the
//! surviving population's Δϖ confinement against the observed samples.

use p9_2016_evidence::scattered_disk_sim::{run_scattered_disk, DiskSimConfig};
use p9_core::constants::*;
use p9_core::types::P9Params;
use rayon::prelude::*;

use crate::clustering_metric::{evaluate_clustering, observed_etno_r_bar};

/// A single point in the parameter grid.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct GridPoint {
    pub mass_earth: f64,
    pub a: f64,
    pub e: f64,
    pub perihelion: f64,
}

impl GridPoint {
    pub fn to_p9_params(&self) -> P9Params {
        P9Params {
            mass_earth: self.mass_earth,
            a: self.a,
            e: self.e,
            i: 0.0,
            omega: 0.0,
            omega_big: std::f64::consts::PI, // Anti-aligned with KBOs
            mean_anomaly: 0.0,
        }
    }
}

/// Generate the full parameter grid from the paper.
///
/// Filters out unphysical combinations where perihelion < 30 AU
/// (would have been detected) or aphelion < 200 AU (too close).
pub fn generate_grid() -> Vec<GridPoint> {
    let masses = [0.1, 1.0, 10.0, 20.0, 30.0];
    let a_values: Vec<f64> = (2..=20).map(|i| i as f64 * 100.0).collect();
    let e_values: Vec<f64> = (1..=9).map(|i| i as f64 * 0.1).collect();

    let mut grid = Vec::new();

    for &m in &masses {
        for &a in &a_values {
            for &e in &e_values {
                let q = a * (1.0 - e);
                let big_q = a * (1.0 + e);

                // Skip if perihelion too small (would be detected)
                // or if the orbit doesn't extend far enough
                if q < 30.0 || big_q < 200.0 {
                    continue;
                }

                grid.push(GridPoint {
                    mass_earth: m,
                    a,
                    e,
                    perihelion: q,
                });
            }
        }
    }

    grid
}

/// Generate a reduced grid for testing (fewer points).
pub fn generate_test_grid() -> Vec<GridPoint> {
    let masses = [10.0];
    let a_values = [400.0, 700.0, 1000.0];
    let e_values = [0.3, 0.6];

    let mut grid = Vec::new();

    for &m in &masses {
        for &a in &a_values {
            for &e in &e_values {
                let q = a * (1.0 - e);
                if q < 30.0 {
                    continue;
                }
                grid.push(GridPoint {
                    mass_earth: m,
                    a,
                    e,
                    perihelion: q,
                });
            }
        }
    }

    grid
}

/// Result of evaluating a single grid point.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GridResult {
    pub point: GridPoint,
    /// Mean resultant length R̄ of the survivors' Δϖ relative to P9.
    pub r_bar_dvarpi: f64,
    /// Fraction of surviving particles with anti-aligned Δϖ (|Δϖ| > π/2)
    pub clustering_fraction: f64,
    /// Fraction of survivors detached to high perihelion
    /// (q > `clustering_metric::HIGH_Q_THRESHOLD_AU` = 60 AU, Sedna-like)
    pub high_perihelion_fraction: f64,
    /// Number of qualifying particles surviving to end of integration
    pub n_survivors: usize,
    /// Total particles at start
    pub n_total: usize,
    /// Whether this parameter set is "accepted"
    pub accepted: bool,
    /// RNG seed used for the simulation (recorded for reproducibility).
    pub seed: u64,
}

/// Minimum survivor count for a grid point to be evaluable. **Assumption**:
/// the paper does not publish an explicit cut; 7 (one more than the observed
/// sample size, so a 6-draw is non-trivial) marks runs with enough survivors
/// to measure confinement at all.
pub const MIN_SURVIVORS: usize = 7;

/// Acceptance criterion for a grid point, recast as a comparison of the
/// simulated confinement statistics against the observed sample rather than
/// invented thresholds:
///
/// 1. ≥ `MIN_SURVIVORS` qualifying survivors (documented assumption above);
/// 2. the survivors' Δϖ mean resultant length reaches the *observed* ETNO
///    sample's R̄ (computed from `p9_core::data::etno`, not a fixed 0.5);
/// 3. the run produces at least one Sedna-like detached object
///    (q > 60 AU; the previous 5% floor was untraceable to the paper).
pub fn evaluate_acceptance(
    r_bar_dvarpi: f64,
    high_perihelion_fraction: f64,
    n_survivors: usize,
) -> bool {
    n_survivors >= MIN_SURVIVORS
        && r_bar_dvarpi >= observed_etno_r_bar()
        && high_perihelion_fraction > 0.0
}

/// Scale settings for the survey runner. The defaults are reduced-scale so
/// the grid is tractable in tests; `full_scale` reproduces the paper's
/// 400-particle, 4-Gyr planar runs (hours of CPU — gate behind `#[ignore]`
/// or an explicit binary).
#[derive(Debug, Clone, Copy)]
pub struct SurveyRunConfig {
    pub n_particles: usize,
    /// Integration horizon (days).
    pub t_total: f64,
    /// Timestep (days); must resolve the particles' q ≈ 30 AU perihelia.
    pub dt: f64,
    /// Clustering analysis window in semi-major axis (AU).
    pub a_min_cluster: f64,
    pub a_max_cluster: f64,
}

impl SurveyRunConfig {
    /// Reduced scale: 24 particles, 1 Myr. Captures the differential apsidal
    /// precession and early sculpting, not the 4-Gyr steady state.
    pub fn reduced() -> Self {
        Self {
            n_particles: 24,
            t_total: 1e6 * YEAR_DAYS,
            dt: 2000.0,
            a_min_cluster: 300.0,
            a_max_cluster: 700.0,
        }
    }

    /// Paper scale: 400 particles, 4 Gyr.
    pub fn full_scale() -> Self {
        Self {
            n_particles: 400,
            t_total: 4.0 * GYR_DAYS,
            dt: 1000.0,
            a_min_cluster: 300.0,
            a_max_cluster: 700.0,
        }
    }
}

/// Run the planar scattered-disk simulation for one grid point and evaluate
/// the clustering acceptance (this is the runner that fills `GridResult`;
/// previously the struct was scaffolding no code ever populated).
pub fn run_grid_point(point: &GridPoint, run: &SurveyRunConfig, seed: u64) -> GridResult {
    let mut config = DiskSimConfig::planar_nominal();
    config.p9 = point.to_p9_params();
    config.n_particles = run.n_particles;
    config.t_total = run.t_total;
    config.dt = run.dt;
    config.snapshot_interval = run.t_total;

    let snapshots = run_scattered_disk(&config, seed);
    let last = snapshots.last().expect("at least the initial snapshot");

    let eval = evaluate_clustering(
        &last.elements,
        last.varpi_p9,
        run.a_min_cluster,
        run.a_max_cluster,
    );

    GridResult {
        point: *point,
        r_bar_dvarpi: eval.r_bar_dvarpi,
        clustering_fraction: eval.anti_aligned_fraction,
        high_perihelion_fraction: eval.high_perihelion_fraction,
        n_survivors: eval.n_qualifying,
        n_total: last.total_count,
        accepted: evaluate_acceptance(
            eval.r_bar_dvarpi,
            eval.high_perihelion_fraction,
            eval.n_qualifying,
        ),
        seed,
    }
}

/// Run the survey over a grid (in parallel). Each point gets a deterministic
/// seed derived from `base_seed` and its index.
pub fn run_grid(grid: &[GridPoint], run: &SurveyRunConfig, base_seed: u64) -> Vec<GridResult> {
    grid.par_iter()
        .enumerate()
        .map(|(k, point)| run_grid_point(point, run, base_seed.wrapping_add(k as u64)))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::clustering_metric::HIGH_Q_THRESHOLD_AU;

    #[test]
    fn test_grid_generation() {
        let grid = generate_grid();
        assert!(!grid.is_empty());

        // All points should have perihelion > 30 AU
        for point in &grid {
            assert!(
                point.perihelion >= 30.0,
                "q = {:.1} should be >= 30",
                point.perihelion
            );
        }

        // Should contain the nominal point (700, 0.6, 10)
        let has_nominal = grid.iter().any(|p| {
            (p.a - 700.0).abs() < 1.0
                && (p.e - 0.6).abs() < 0.01
                && (p.mass_earth - 10.0).abs() < 0.1
        });
        assert!(has_nominal, "Grid should include nominal P9 parameters");
    }

    #[test]
    fn test_grid_size() {
        let grid = generate_grid();
        // Should be a reasonable number of points
        assert!(
            grid.len() > 100 && grid.len() < 2000,
            "Grid has {} points",
            grid.len()
        );
    }

    #[test]
    fn test_acceptance_criteria() {
        let r_obs = observed_etno_r_bar();
        // Confined at the observed level, with detached objects
        assert!(evaluate_acceptance(r_obs + 0.05, 0.1, 10));
        // Not enough survivors
        assert!(!evaluate_acceptance(r_obs + 0.05, 0.1, MIN_SURVIVORS - 1));
        // Less confined than the observed sample
        assert!(!evaluate_acceptance(r_obs - 0.2, 0.1, 10));
        // No detached production
        assert!(!evaluate_acceptance(r_obs + 0.05, 0.0, 10));
    }

    #[test]
    fn test_high_q_threshold_single_source() {
        // L3 regression: the documented Sedna-like threshold and the metric
        // code share one constant (previously docs said 100 AU, code 60 AU).
        assert_eq!(HIGH_Q_THRESHOLD_AU, 60.0);
    }

    #[test]
    fn test_run_grid_point_reduced() {
        // Real (reduced-scale) runner on the nominal point: the GridResult
        // is filled from an actual simulation, deterministically.
        let point = GridPoint {
            mass_earth: 10.0,
            a: 700.0,
            e: 0.6,
            perihelion: 280.0,
        };
        let mut run = SurveyRunConfig::reduced();
        run.n_particles = 12;
        run.t_total = 2e5 * YEAR_DAYS;

        let result = run_grid_point(&point, &run, 17);
        assert_eq!(result.seed, 17);
        assert!(result.n_total >= 12); // planar generator may round up
        assert!(result.n_survivors <= result.n_total);
        assert!(result.r_bar_dvarpi >= 0.0 && result.r_bar_dvarpi <= 1.0);

        // Deterministic under the same seed
        let result2 = run_grid_point(&point, &run, 17);
        assert_eq!(result.n_survivors, result2.n_survivors);
        assert_eq!(result.r_bar_dvarpi, result2.r_bar_dvarpi);
    }

    #[test]
    #[ignore = "full-scale survey over the test grid (4 Gyr, 400 particles per point; hours)"]
    fn test_run_test_grid_full_scale() {
        let grid = generate_test_grid();
        let results = run_grid(&grid, &SurveyRunConfig::full_scale(), 99);
        assert_eq!(results.len(), grid.len());
        // The nominal-like point (a=700, e=0.6) should be among the better
        // confined points of this small grid.
        let nominal = results
            .iter()
            .find(|r| (r.point.a - 700.0).abs() < 1.0 && (r.point.e - 0.6).abs() < 0.01)
            .unwrap();
        let max_r = results.iter().map(|r| r.r_bar_dvarpi).fold(0.0, f64::max);
        assert!(nominal.r_bar_dvarpi >= 0.5 * max_r);
    }
}
