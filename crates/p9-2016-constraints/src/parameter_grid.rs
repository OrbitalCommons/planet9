//! Parameter grid for the planar survey.
//!
//! Brown & Batygin (2016) systematically explore:
//!   a₉ ∈ (200, 2000) AU in 100 AU steps → 19 values
//!   e₉ ∈ (0.1, 0.9) in 0.1 steps → 9 values
//!   m₉ ∈ {0.1, 1, 10, 20, 30} M_Earth → 5 values
//!
//! Total: 19 × 9 × 5 = 855 parameter combinations
//! (paper says ~320 per mass × 5 masses ≈ 1600; they may skip some unphysical)

use p9_core::types::P9Params;

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
    /// Fraction of surviving particles with anti-aligned Δϖ
    pub clustering_fraction: f64,
    /// Fraction of survivors with elevated perihelion (q > 100 AU)
    pub high_perihelion_fraction: f64,
    /// Number of particles surviving to end of integration
    pub n_survivors: usize,
    /// Total particles at start
    pub n_total: usize,
    /// Whether this parameter set is "accepted"
    pub accepted: bool,
}

/// The paper's acceptance criteria:
/// 1. ≥7 survivors with a ∈ [300,700], q < 80 AU, i < 50°
/// 2. Perihelion confinement: most (>50%) survivors have |Δϖ| > π/2 (anti-aligned)
/// 3. Some fraction of objects pushed to high perihelion (Sedna-like, q > 60 AU)
pub fn evaluate_acceptance(
    clustering_fraction: f64,
    high_perihelion_fraction: f64,
    n_survivors: usize,
) -> bool {
    n_survivors >= 7 && clustering_fraction > 0.5 && high_perihelion_fraction > 0.05
}

#[cfg(test)]
mod tests {
    use super::*;

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
        // Good clustering
        assert!(evaluate_acceptance(0.7, 0.1, 10));
        // Not enough survivors
        assert!(!evaluate_acceptance(0.7, 0.1, 5));
        // Poor clustering
        assert!(!evaluate_acceptance(0.3, 0.1, 10));
    }
}
