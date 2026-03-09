//! MCMC posterior parameters from Brown & Batygin (2021).
//!
//! The paper derives posterior distributions for Planet Nine's orbital
//! parameters using Markov Chain Monte Carlo analysis of 11 distant KBOs.

use rand::Rng;
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;

/// Posterior distribution for a single orbital parameter,
/// represented as an asymmetric Gaussian (median with upper/lower 1-sigma).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AsymmetricGaussian {
    pub median: f64,
    pub sigma_upper: f64,
    pub sigma_lower: f64,
}

impl AsymmetricGaussian {
    /// Sample from the asymmetric Gaussian using a split-normal distribution.
    /// Draws from the upper half-Gaussian if the uniform variate > 0.5,
    /// otherwise from the lower half-Gaussian.
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        let u: f64 = rng.gen();
        if u < 0.5 {
            let normal = Normal::new(0.0, self.sigma_lower).unwrap();
            self.median - normal.sample(rng).abs()
        } else {
            let normal = Normal::new(0.0, self.sigma_upper).unwrap();
            self.median + normal.sample(rng).abs()
        }
    }

    /// Sample with truncation at given min/max bounds.
    pub fn sample_truncated<R: Rng>(&self, rng: &mut R, min: f64, max: f64) -> f64 {
        loop {
            let val = self.sample(rng);
            if val >= min && val <= max {
                return val;
            }
        }
    }
}

/// Full MCMC posterior for Planet Nine orbital parameters.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct P9Posterior {
    /// Mass in Earth masses
    pub mass: AsymmetricGaussian,
    /// Semi-major axis in AU
    pub a: AsymmetricGaussian,
    /// Eccentricity (derived from a and q)
    pub e: AsymmetricGaussian,
    /// Inclination in degrees
    pub i: AsymmetricGaussian,
    /// Perihelion distance in AU
    pub perihelion: AsymmetricGaussian,
    /// Mean of the full posterior (for reference)
    pub mean: PosteriorSummary,
    /// Standard deviations of the full posterior
    pub std: PosteriorSummary,
}

/// Summary statistics (mean or std) for each parameter.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PosteriorSummary {
    pub mass_earth: f64,
    pub a_au: f64,
    pub e: f64,
    pub i_deg: f64,
    pub q_au: f64,
}

/// Returns the MCMC posterior from Brown & Batygin (2021).
///
/// Parameter constraints from the paper:
/// - Mass: 6.2 (+2.2/-1.3) Earth masses
/// - Semi-major axis: 380 (+140/-80) AU
/// - Inclination: 16 +/- 5 degrees
/// - Perihelion: 300 (+85/-60) AU
/// - Eccentricity: ~0.21 (derived from median a and q)
pub fn mcmc_2021_posterior() -> P9Posterior {
    let median_a = 380.0;
    let median_q = 300.0;
    let median_e = 1.0 - median_q / median_a;

    P9Posterior {
        mass: AsymmetricGaussian {
            median: 6.2,
            sigma_upper: 2.2,
            sigma_lower: 1.3,
        },
        a: AsymmetricGaussian {
            median: 380.0,
            sigma_upper: 140.0,
            sigma_lower: 80.0,
        },
        e: AsymmetricGaussian {
            median: median_e,
            sigma_upper: 0.15,
            sigma_lower: 0.10,
        },
        i: AsymmetricGaussian {
            median: 16.0,
            sigma_upper: 5.0,
            sigma_lower: 5.0,
        },
        perihelion: AsymmetricGaussian {
            median: 300.0,
            sigma_upper: 85.0,
            sigma_lower: 60.0,
        },
        mean: PosteriorSummary {
            mass_earth: 6.2,
            a_au: 380.0,
            e: median_e,
            i_deg: 16.0,
            q_au: 300.0,
        },
        std: PosteriorSummary {
            mass_earth: 1.75,
            a_au: 110.0,
            e: 0.12,
            i_deg: 5.0,
            q_au: 72.5,
        },
    }
}

/// Sample a set of P9 orbital parameters from the 2021 posterior.
///
/// Uses truncated asymmetric Gaussians to ensure physical validity:
/// - Mass > 0
/// - Semi-major axis > 0
/// - 0 < eccentricity < 1
/// - 0 < inclination < 90 degrees
/// - Perihelion > 0 and < semi-major axis
///
/// The argument of perihelion and longitude of ascending node are
/// drawn uniformly (the paper constrains longitude of perihelion
/// but not these individually).
pub fn sample_from_posterior<R: Rng>(posterior: &P9Posterior, rng: &mut R) -> P9Params {
    let mass = posterior.mass.sample_truncated(rng, 1.0, 20.0);
    let a = posterior.a.sample_truncated(rng, 200.0, 1000.0);
    let i_deg = posterior.i.sample_truncated(rng, 0.0, 90.0);
    let q = posterior.perihelion.sample_truncated(rng, 50.0, a);
    let e = 1.0 - q / a;

    let omega: f64 = rng.gen_range(0.0..std::f64::consts::TAU);
    let omega_big: f64 = rng.gen_range(0.0..std::f64::consts::TAU);
    let mean_anomaly: f64 = rng.gen_range(0.0..std::f64::consts::TAU);

    P9Params {
        mass_earth: mass,
        a,
        e,
        i: i_deg * DEG2RAD,
        omega,
        omega_big,
        mean_anomaly,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;

    #[test]
    fn test_posterior_median_values() {
        let post = mcmc_2021_posterior();
        assert_relative_eq!(post.mass.median, 6.2, epsilon = 0.01);
        assert_relative_eq!(post.a.median, 380.0, epsilon = 0.1);
        assert_relative_eq!(post.i.median, 16.0, epsilon = 0.1);
        assert_relative_eq!(post.perihelion.median, 300.0, epsilon = 0.1);
    }

    #[test]
    fn test_posterior_eccentricity_derived() {
        let post = mcmc_2021_posterior();
        let expected_e = 1.0 - 300.0 / 380.0;
        assert_relative_eq!(post.e.median, expected_e, epsilon = 0.01);
    }

    #[test]
    fn test_sample_physical_validity() {
        let post = mcmc_2021_posterior();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        for _ in 0..100 {
            let params = sample_from_posterior(&post, &mut rng);
            assert!(params.mass_earth > 0.0, "Mass must be positive");
            assert!(params.a > 0.0, "Semi-major axis must be positive");
            assert!(
                params.e > 0.0 && params.e < 1.0,
                "Eccentricity must be in (0,1)"
            );
            assert!(params.i > 0.0, "Inclination must be positive");
            assert!(params.perihelion() > 0.0, "Perihelion must be positive");
            assert!(
                params.perihelion() < params.a,
                "Perihelion must be less than a"
            );
        }
    }

    #[test]
    fn test_sample_distribution_reasonable() {
        let post = mcmc_2021_posterior();
        let mut rng = rand::rngs::StdRng::seed_from_u64(123);
        let n = 1000;

        let samples: Vec<P9Params> = (0..n)
            .map(|_| sample_from_posterior(&post, &mut rng))
            .collect();

        let mean_mass: f64 = samples.iter().map(|s| s.mass_earth).sum::<f64>() / n as f64;
        let mean_a: f64 = samples.iter().map(|s| s.a).sum::<f64>() / n as f64;

        assert!(
            (mean_mass - 6.2).abs() < 1.5,
            "Mean mass {:.2} should be near 6.2",
            mean_mass
        );
        assert!(
            (mean_a - 380.0).abs() < 80.0,
            "Mean a {:.1} should be near 380",
            mean_a
        );
    }
}
