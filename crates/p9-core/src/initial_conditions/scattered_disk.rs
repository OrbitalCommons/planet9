//! Synthetic scattered disk generators.
//!
//! Creates initial populations of test particles matching the
//! distributions used in the Planet Nine papers:
//! - Uniform semi-major axis in a given range
//! - Uniform perihelion distance in a given range
//! - Half-normal inclination distribution
//! - Uniform random angles (ω, Ω, M)

use rand::Rng;
use rand_distr::{Distribution, Normal, Uniform};

use crate::constants::*;
use crate::types::{elements_to_cartesian, OrbitalElements, StateVector};

/// Configuration for generating a synthetic scattered disk population.
#[derive(Debug, Clone)]
pub struct ScatteredDiskConfig {
    /// Semi-major axis range (AU)
    pub a_min: f64,
    pub a_max: f64,
    /// Perihelion distance range (AU)
    pub q_min: f64,
    pub q_max: f64,
    /// Standard deviation of the half-normal inclination distribution (radians)
    pub sigma_i: f64,
    /// Number of particles to generate
    pub n_particles: usize,
}

impl ScatteredDiskConfig {
    /// Standard config from Batygin & Brown (2016): evidence paper
    pub fn batygin_brown_2016() -> Self {
        Self {
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 15.0 * DEG2RAD,
            n_particles: 400,
        }
    }

    /// Config for the inclined TNOs paper (Batygin & Brown 2016b)
    pub fn inclined_tnos_2016() -> Self {
        Self {
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 50.0,
            sigma_i: 15.0 * DEG2RAD,
            n_particles: 3200,
        }
    }

    /// Broad perihelion distribution (Khain et al. 2018)
    pub fn broad_perihelion_2018() -> Self {
        Self {
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 300.0,
            sigma_i: 15.0 * DEG2RAD,
            n_particles: 400,
        }
    }

    /// Narrow perihelion distribution (Khain et al. 2018)
    pub fn narrow_perihelion_2018() -> Self {
        Self {
            a_min: 150.0,
            a_max: 550.0,
            q_min: 30.0,
            q_max: 36.0,
            sigma_i: 15.0 * DEG2RAD,
            n_particles: 400,
        }
    }
}

/// Generate a synthetic scattered disk population.
///
/// Returns state vectors in heliocentric ecliptic J2000 (AU, AU/day).
pub fn generate_scattered_disk<R: Rng>(
    config: &ScatteredDiskConfig,
    rng: &mut R,
) -> Vec<StateVector> {
    let a_dist = Uniform::new(config.a_min, config.a_max);
    let q_dist = Uniform::new(config.q_min, config.q_max);
    let angle_dist = Uniform::new(0.0, TWO_PI);
    let incl_normal = Normal::new(0.0, config.sigma_i).unwrap();

    let mut particles = Vec::with_capacity(config.n_particles);

    for _ in 0..config.n_particles {
        let a = a_dist.sample(rng);
        let q = q_dist.sample(rng);

        // Eccentricity from a and q: e = 1 - q/a
        let e = 1.0 - q / a;

        // Skip invalid configurations (q > a means e < 0, or q = 0)
        if e < 0.0 || e >= 1.0 {
            continue;
        }

        // Half-normal inclination: take absolute value of normal sample
        let i = incl_normal.sample(rng).abs();

        let omega = angle_dist.sample(rng);
        let omega_big = angle_dist.sample(rng);
        let mean_anomaly = angle_dist.sample(rng);

        let elements = OrbitalElements {
            a,
            e,
            i,
            omega_big,
            omega,
            mean_anomaly,
        };

        particles.push(elements_to_cartesian(&elements, GM_SUN));
    }

    particles
}

/// Generate a planar scattered disk (i = 0) for phase-space studies.
pub fn generate_planar_disk<R: Rng>(
    a_min: f64,
    a_max: f64,
    q_min: f64,
    q_max: f64,
    n_particles: usize,
    rng: &mut R,
) -> Vec<StateVector> {
    let config = ScatteredDiskConfig {
        a_min,
        a_max,
        q_min,
        q_max,
        sigma_i: 0.0,
        n_particles,
    };

    let a_dist = Uniform::new(config.a_min, config.a_max);
    let q_dist = Uniform::new(config.q_min, config.q_max);
    let angle_dist = Uniform::new(0.0, TWO_PI);

    let mut particles = Vec::with_capacity(n_particles);

    for _ in 0..n_particles {
        let a = a_dist.sample(rng);
        let q = q_dist.sample(rng);
        let e = 1.0 - q / a;

        if e < 0.0 || e >= 1.0 {
            continue;
        }

        let omega = angle_dist.sample(rng);
        let mean_anomaly = angle_dist.sample(rng);

        let elements = OrbitalElements {
            a,
            e,
            i: 0.0,
            omega_big: 0.0,
            omega,
            mean_anomaly,
        };

        particles.push(elements_to_cartesian(&elements, GM_SUN));
    }

    particles
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::types::cartesian_to_elements;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_scattered_disk_generation() {
        let mut rng = StdRng::seed_from_u64(42);
        let config = ScatteredDiskConfig::batygin_brown_2016();
        let particles = generate_scattered_disk(&config, &mut rng);

        assert!(!particles.is_empty());

        // Check that all particles have valid orbital elements within specified ranges
        for p in &particles {
            let elem = cartesian_to_elements(p, GM_SUN);
            assert!(
                elem.a >= config.a_min * 0.99 && elem.a <= config.a_max * 1.01,
                "a = {} out of range [{}, {}]",
                elem.a,
                config.a_min,
                config.a_max
            );

            let q = elem.perihelion();
            assert!(
                q >= config.q_min * 0.99 && q <= config.q_max * 1.01,
                "q = {} out of range [{}, {}]",
                q,
                config.q_min,
                config.q_max
            );

            assert!(
                elem.e >= 0.0 && elem.e < 1.0,
                "e = {} out of valid range",
                elem.e
            );

            assert!(elem.i >= 0.0, "i = {} should be non-negative", elem.i);
        }
    }

    #[test]
    fn test_inclination_distribution() {
        let mut rng = StdRng::seed_from_u64(123);
        let config = ScatteredDiskConfig {
            n_particles: 10000,
            ..ScatteredDiskConfig::batygin_brown_2016()
        };
        let particles = generate_scattered_disk(&config, &mut rng);

        // Compute mean inclination — for half-normal with sigma=15°, E[|X|] = sigma*sqrt(2/pi)
        let expected_mean = config.sigma_i * (2.0 / std::f64::consts::PI).sqrt();
        let actual_mean: f64 = particles
            .iter()
            .map(|p| cartesian_to_elements(p, GM_SUN).i)
            .sum::<f64>()
            / particles.len() as f64;

        let relative_error = (actual_mean - expected_mean).abs() / expected_mean;
        assert!(
            relative_error < 0.05,
            "Mean inclination {:.4} vs expected {:.4} (error {:.2}%)",
            actual_mean * RAD2DEG,
            expected_mean * RAD2DEG,
            relative_error * 100.0
        );
    }
}
