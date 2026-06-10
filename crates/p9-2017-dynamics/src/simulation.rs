//! Simulation configuration and test particle setup.
//!
//! Fiducial simulation: 6,000 test particles, a in (150-750) AU, q in (30-36) AU,
//! integrated for 4 Gyr. Removal at r < 20 AU or r > 10,000 AU.
//! Observability cuts: q <= 100 AU AND i <= 40 deg.

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;
use rand::Rng;
use rand_distr::{Distribution, Uniform};
use serde::{Deserialize, Serialize};

/// Simulation configuration matching the paper.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    /// Planet Nine parameters
    pub p9: P9Params,
    /// Number of test particles
    pub n_particles: usize,
    /// Integration time in Gyr
    pub t_gyr: f64,
    /// Minimum initial semi-major axis (AU)
    pub a_min: f64,
    /// Maximum initial semi-major axis (AU)
    pub a_max: f64,
    /// Minimum initial perihelion distance (AU)
    pub q_min: f64,
    /// Maximum initial perihelion distance (AU)
    pub q_max: f64,
    /// Inner removal boundary (AU)
    pub r_remove_inner: f64,
    /// Outer removal boundary (AU)
    pub r_remove_outer: f64,
}

impl SimConfig {
    /// Fiducial configuration: a_9=700, e_9=0.6, m_9=10 ME.
    pub fn default_paper() -> Self {
        Self {
            p9: P9Params {
                mass_earth: 10.0,
                a: 700.0,
                e: 0.6,
                i: 20.0 * DEG2RAD,
                omega: 150.0 * DEG2RAD,
                omega_big: 100.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            n_particles: 6_000,
            t_gyr: 4.0,
            a_min: 150.0,
            a_max: 750.0,
            q_min: 30.0,
            q_max: 36.0,
            r_remove_inner: 20.0,
            r_remove_outer: 10_000.0,
        }
    }

    /// Mass variation: m_9 = 5 Earth masses.
    pub fn mass_5() -> Self {
        let mut config = Self::default_paper();
        config.p9.mass_earth = 5.0;
        config
    }

    /// Mass variation: m_9 = 20 Earth masses.
    pub fn mass_20() -> Self {
        let mut config = Self::default_paper();
        config.p9.mass_earth = 20.0;
        config
    }
}

/// Observability selection criteria from the paper.
pub struct ObservabilityCriteria {
    pub q_max: f64,
    pub i_max_deg: f64,
}

impl ObservabilityCriteria {
    pub fn default_paper() -> Self {
        Self {
            q_max: 100.0,
            i_max_deg: 40.0,
        }
    }

    pub fn is_observable(&self, q: f64, i_deg: f64) -> bool {
        q <= self.q_max && i_deg <= self.i_max_deg
    }
}

/// Generate initial orbital elements for test particles.
///
/// Semi-major axes uniform in (a_min, a_max), perihelion distances
/// uniform in (q_min, q_max), inclinations drawn from Rayleigh distribution
/// with sigma = 5 deg.
pub fn generate_initial_conditions<R: Rng>(config: &SimConfig, rng: &mut R) -> Vec<TestParticle> {
    let a_dist = Uniform::new(config.a_min, config.a_max);
    let q_dist = Uniform::new(config.q_min, config.q_max);
    let angle_dist = Uniform::new(0.0_f64, std::f64::consts::TAU);

    let mut particles = Vec::with_capacity(config.n_particles);

    for _ in 0..config.n_particles {
        let a = a_dist.sample(rng);
        let q = q_dist.sample(rng);
        let e = 1.0 - q / a;

        if e < 0.0 || e >= 1.0 {
            continue;
        }

        // Rayleigh distribution for inclination with sigma = 5 deg
        let u: f64 = rng.gen();
        let i = 5.0 * DEG2RAD * (-2.0 * (1.0 - u).ln()).sqrt();

        let omega = angle_dist.sample(rng);
        let omega_big = angle_dist.sample(rng);
        let mean_anomaly = angle_dist.sample(rng);

        particles.push(TestParticle {
            a,
            e,
            i,
            omega,
            omega_big,
            mean_anomaly,
        });
    }

    particles
}

/// Test particle orbital elements.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestParticle {
    pub a: f64,
    pub e: f64,
    pub i: f64,
    pub omega: f64,
    pub omega_big: f64,
    pub mean_anomaly: f64,
}

impl TestParticle {
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    pub fn inclination_deg(&self) -> f64 {
        self.i / DEG2RAD
    }
}

/// Quick test simulation producing synthetic end-state particles.
pub fn quick_test_simulation(config: &SimConfig) -> Vec<TestParticle> {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    generate_initial_conditions(config, &mut rng)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = SimConfig::default_paper();
        assert!((config.p9.mass_earth - 10.0).abs() < 0.01);
        assert!((config.p9.a - 700.0).abs() < 0.01);
        assert_eq!(config.n_particles, 6_000);
    }

    #[test]
    fn test_generate_particles() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let config = SimConfig::default_paper();
        let particles = generate_initial_conditions(&config, &mut rng);
        assert!(!particles.is_empty());

        for p in &particles {
            assert!(p.a >= config.a_min && p.a <= config.a_max);
            assert!(p.e >= 0.0 && p.e < 1.0);
            let q = p.perihelion();
            assert!(q >= config.q_min - 1.0 && q <= config.q_max + 1.0);
        }
    }

    #[test]
    fn test_observability_criteria() {
        let criteria = ObservabilityCriteria::default_paper();
        assert!(criteria.is_observable(50.0, 20.0));
        assert!(!criteria.is_observable(150.0, 20.0));
        assert!(!criteria.is_observable(50.0, 50.0));
    }

    #[test]
    fn test_mass_variants() {
        let c5 = SimConfig::mass_5();
        let c20 = SimConfig::mass_20();
        assert!((c5.p9.mass_earth - 5.0).abs() < 0.01);
        assert!((c20.p9.mass_earth - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_high_inclination_fraction() {
        let particles = quick_test_simulation(&SimConfig::default_paper());
        let high_i = particles
            .iter()
            .filter(|p| p.inclination_deg() > 40.0)
            .count();
        let fraction = high_i as f64 / particles.len() as f64;
        // Paper reports 38% experience high-i excursions, but initial conditions
        // should have very few at i > 40 deg (Rayleigh with sigma=5)
        assert!(fraction < 0.1, "initial high-i fraction = {}", fraction);
    }
}
