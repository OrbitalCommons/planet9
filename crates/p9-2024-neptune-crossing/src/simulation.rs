//! Simulation configurations for the P9-inclusive and P9-free models.
//!
//! The P9-inclusive model uses Planet Nine with m = 5 Earth masses, a = 500 AU,
//! e = 0.25, i = 20 deg. The P9-free null model includes only the known giant
//! planets and galactic tide.

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;
use serde::{Deserialize, Serialize};

/// Configuration for the P9-inclusive simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct P9InclusiveConfig {
    /// Planet Nine parameters
    pub p9: P9Params,
    /// Number of test particles
    pub n_particles: usize,
    /// Integration time in Gyr
    pub t_gyr: f64,
}

impl P9InclusiveConfig {
    /// Default configuration matching the paper: m=5 ME, a=500, e=0.25, i=20 deg.
    pub fn default_paper() -> Self {
        Self {
            p9: P9Params {
                mass_earth: 5.0,
                a: 500.0,
                e: 0.25,
                i: 20.0 * DEG2RAD,
                omega: 150.0 * DEG2RAD,
                omega_big: 100.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            n_particles: 10_000,
            t_gyr: 4.0,
        }
    }
}

/// Configuration for the P9-free null model simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct P9FreeConfig {
    /// Number of test particles
    pub n_particles: usize,
    /// Integration time in Gyr
    pub t_gyr: f64,
    /// Include galactic tide perturbation
    pub include_galactic_tide: bool,
}

impl P9FreeConfig {
    /// Default null model configuration.
    pub fn default_paper() -> Self {
        Self {
            n_particles: 10_000,
            t_gyr: 4.0,
            include_galactic_tide: true,
        }
    }
}

/// Result of a simulation run, storing orbital footprint statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationResult {
    /// Number of surviving particles
    pub n_surviving: usize,
    /// Number meeting selection criteria (a>100, i<40, q<30)
    pub n_selected: usize,
    /// Perihelion distances of selected particles (AU)
    pub selected_perihelia: Vec<f64>,
    /// Inclinations of selected particles (degrees)
    pub selected_inclinations: Vec<f64>,
    /// Semi-major axes of selected particles (AU)
    pub selected_semimajor: Vec<f64>,
    /// Model label
    pub model_name: String,
}

impl SimulationResult {
    /// Fraction of particles meeting the selection criteria.
    pub fn selection_fraction(&self) -> f64 {
        if self.n_surviving == 0 {
            return 0.0;
        }
        self.n_selected as f64 / self.n_surviving as f64
    }
}

/// Run a quick test simulation with a small number of particles.
///
/// This generates synthetic orbital elements from simplified distributions
/// to validate the analysis pipeline without full N-body integration.
pub fn quick_test_simulation(with_p9: bool) -> SimulationResult {
    use rand::SeedableRng;
    use rand_distr::{Distribution, Normal, Uniform};

    let mut rng = rand::rngs::StdRng::seed_from_u64(if with_p9 { 42 } else { 137 });
    let n = 500;

    let a_dist = Uniform::new(100.0, 800.0);
    let e_dist = Uniform::new(0.6, 0.99);
    let i_normal = Normal::new(15.0, 12.0).unwrap();

    let mut selected_perihelia = Vec::new();
    let mut selected_inclinations = Vec::new();
    let mut selected_semimajor = Vec::new();

    for _ in 0..n {
        let a = a_dist.sample(&mut rng);
        let e = e_dist.sample(&mut rng);
        let i: f64 = i_normal.sample(&mut rng);
        let i = i.abs();
        let q = a * (1.0 - e);

        // P9 model produces more Neptune-crossers at low inclination
        let i_adjusted = if with_p9 { i * 0.7 } else { i };
        let q_adjusted = if with_p9 { q * 0.85 } else { q };

        if a > 100.0 && i_adjusted < 40.0 && q_adjusted < 30.0 {
            selected_perihelia.push(q_adjusted);
            selected_inclinations.push(i_adjusted);
            selected_semimajor.push(a);
        }
    }

    let n_selected = selected_perihelia.len();
    let model_name = if with_p9 {
        "P9-inclusive".to_string()
    } else {
        "P9-free".to_string()
    };

    SimulationResult {
        n_surviving: n,
        n_selected,
        selected_perihelia,
        selected_inclinations,
        selected_semimajor,
        model_name,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_p9_inclusive_config() {
        let config = P9InclusiveConfig::default_paper();
        assert!((config.p9.mass_earth - 5.0).abs() < 0.01);
        assert!((config.p9.a - 500.0).abs() < 0.01);
        assert!((config.p9.e - 0.25).abs() < 0.01);
        let i_deg = config.p9.i / DEG2RAD;
        assert!((i_deg - 20.0).abs() < 0.1);
    }

    #[test]
    fn test_p9_free_config() {
        let config = P9FreeConfig::default_paper();
        assert_eq!(config.n_particles, 10_000);
        assert!(config.include_galactic_tide);
    }

    #[test]
    fn test_quick_simulation_p9() {
        let result = quick_test_simulation(true);
        assert_eq!(result.model_name, "P9-inclusive");
        assert!(result.n_selected > 0, "should select some particles");
        assert_eq!(result.selected_perihelia.len(), result.n_selected);
    }

    #[test]
    fn test_quick_simulation_null() {
        let result = quick_test_simulation(false);
        assert_eq!(result.model_name, "P9-free");
        assert!(result.n_selected > 0, "should select some particles");
    }

    #[test]
    fn test_selection_fraction() {
        let result = quick_test_simulation(true);
        let frac = result.selection_fraction();
        assert!(frac > 0.0 && frac <= 1.0, "fraction = {}", frac);
    }

    #[test]
    fn test_selected_within_criteria() {
        let result = quick_test_simulation(true);
        for &q in &result.selected_perihelia {
            assert!(q < 30.0, "q = {} should be < 30 AU", q);
        }
        for &i in &result.selected_inclinations {
            assert!(i < 40.0, "i = {} should be < 40 deg", i);
        }
        for &a in &result.selected_semimajor {
            assert!(a > 100.0, "a = {} should be > 100 AU", a);
        }
    }
}
