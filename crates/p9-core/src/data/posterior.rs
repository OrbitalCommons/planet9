//! Summary of the *published* posterior of Brown & Batygin (2021).
//!
//! This module does not re-run the inference: the published posterior
//! summaries (median with asymmetric 1-sigma intervals) are emulated with
//! two-piece (split-normal) distributions, plus a documented correlation
//! assumption linking a and q (and hence e) — the published corner plot
//! shows a strong positive a–q degeneracy which independent marginal
//! sampling would erase. It is the stable representation of the paper's
//! *result*, consumed by the downstream survey crates via
//! `reference_population`.
//!
//! The from-scratch path (the paper's models 1-4: simulation grid -> KDE
//! likelihood -> GP emulator -> ensemble MCMC) lives in `sim_grid`/`kde`/
//! `gp`/`mcmc`, wired in `pipeline`. At the default reduced scale that
//! pipeline validates machinery, not the published numbers (the
//! `#[ignore]`d reduced end-to-end test checks the published medians fall
//! inside its loose 95% region; full fidelity needs the paper-scale grid —
//! see REPRODUCTION_NOTES.md).

use rand::Rng;
use rand_distr::{Distribution, Normal};
use serde::{Deserialize, Serialize};

use crate::constants::DEG2RAD;
use crate::types::P9Params;

/// Assumed correlation between the a and q marginals when sampling jointly.
///
/// Documented emulator assumption: Brown & Batygin (2021, Fig. 7) show a
/// strong positive a–q degeneracy (more distant solutions require larger
/// perihelia); 0.85 reproduces the visual elongation of that contour. This
/// is an assumption of the emulator, not a published number.
pub const A_Q_CORRELATION: f64 = 0.85;

/// Posterior distribution for a single orbital parameter, represented as a
/// two-piece normal (mode with upper/lower 1-sigma widths).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct AsymmetricGaussian {
    pub median: f64,
    pub sigma_upper: f64,
    pub sigma_lower: f64,
}

impl AsymmetricGaussian {
    /// Sample from the two-piece normal. The side is chosen with the proper
    /// weights — lower side with probability σ_l/(σ_l + σ_u) — so the
    /// density is continuous at the mode (a 50/50 split, the previous
    /// behavior, over-weights the narrow side and biases the mean).
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        let w_lower = self.sigma_lower / (self.sigma_lower + self.sigma_upper);
        if rng.gen::<f64>() < w_lower {
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

    /// Mean of the (untruncated) two-piece normal:
    /// mode + sqrt(2/π)(σ_u − σ_l).
    pub fn mean(&self) -> f64 {
        self.median + (2.0 / std::f64::consts::PI).sqrt() * (self.sigma_upper - self.sigma_lower)
    }
}

/// Full posterior summary for Planet Nine orbital parameters.
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

/// Returns the published posterior summaries from Brown & Batygin (2021).
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

/// Sample a set of P9 orbital parameters from the 2021 posterior emulator.
///
/// a and q (hence e) are sampled *jointly*: a is drawn from its marginal,
/// then q from its conditional given a, using a linear (Gaussian-style)
/// regression with correlation `A_Q_CORRELATION` applied to the two-piece
/// summaries:
///
///   q | a ~ TwoPiece(q_med + ρ (σ_q/σ_a)(a − a_med), σ's shrunk by √(1−ρ²))
///
/// Truncations keep the sample physical (q < a so 0 < e < 1, etc.). Mass and
/// inclination are drawn from their marginals; the angles the paper does not
/// constrain individually (ω, Ω, M) are uniform.
pub fn sample_from_posterior<R: Rng>(posterior: &P9Posterior, rng: &mut R) -> P9Params {
    let mass = posterior.mass.sample_truncated(rng, 1.0, 20.0);
    let a = posterior.a.sample_truncated(rng, 200.0, 1000.0);
    let i_deg = posterior.i.sample_truncated(rng, 0.0, 90.0);

    // Conditional q | a with the documented correlation assumption.
    let rho = A_Q_CORRELATION;
    let sigma_a = 0.5 * (posterior.a.sigma_upper + posterior.a.sigma_lower);
    let sigma_q = 0.5 * (posterior.perihelion.sigma_upper + posterior.perihelion.sigma_lower);
    let shrink = (1.0 - rho * rho).sqrt();
    let conditional = AsymmetricGaussian {
        median: posterior.perihelion.median + rho * (sigma_q / sigma_a) * (a - posterior.a.median),
        sigma_upper: posterior.perihelion.sigma_upper * shrink,
        sigma_lower: posterior.perihelion.sigma_lower * shrink,
    };
    let q = conditional.sample_truncated(rng, 50.0, 0.999 * a);
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
    use uom::si::length::astronomical_unit;

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

    /// H5 regression: the two-piece sampler must put mass
    /// σ_l/(σ_l + σ_u) below the mode (a 50/50 side split biases the mean
    /// low for σ_u > σ_l).
    #[test]
    fn test_two_piece_side_weights() {
        let g = AsymmetricGaussian {
            median: 0.0,
            sigma_upper: 3.0,
            sigma_lower: 1.0,
        };
        let mut rng = rand::rngs::StdRng::seed_from_u64(5);
        let n = 200_000;
        let below = (0..n).filter(|_| g.sample(&mut rng) < 0.0).count();
        let frac_below = below as f64 / n as f64;
        let expected = 1.0 / 4.0; // sigma_l / (sigma_l + sigma_u)
        assert!(
            (frac_below - expected).abs() < 0.005,
            "frac below mode = {frac_below:.4}, expected {expected}"
        );

        // And the sample mean matches the closed-form two-piece mean.
        let mean_mc: f64 = (0..n).map(|_| g.sample(&mut rng)).sum::<f64>() / n as f64;
        assert!(
            (mean_mc - g.mean()).abs() < 0.02,
            "MC mean {mean_mc:.3} vs analytic {:.3}",
            g.mean()
        );
    }

    /// H5: a and q must come out correlated (joint sampling), not
    /// independent.
    #[test]
    fn test_a_q_jointly_sampled() {
        let post = mcmc_2021_posterior();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let samples: Vec<P9Params> = (0..5000)
            .map(|_| sample_from_posterior(&post, &mut rng))
            .collect();

        let qs: Vec<f64> = samples
            .iter()
            .map(|s| s.perihelion_typed().get::<astronomical_unit>())
            .collect();
        let a_mean = samples.iter().map(|s| s.a).sum::<f64>() / samples.len() as f64;
        let q_mean = qs.iter().sum::<f64>() / qs.len() as f64;
        let mut sa = 0.0;
        let mut sq = 0.0;
        let mut saq = 0.0;
        for (s, &q) in samples.iter().zip(&qs) {
            let da = s.a - a_mean;
            let dq = q - q_mean;
            sa += da * da;
            sq += dq * dq;
            saq += da * dq;
        }
        let corr = saq / (sa.sqrt() * sq.sqrt());
        assert!(
            (corr - A_Q_CORRELATION).abs() < 0.1,
            "sample corr(a, q) = {corr:.3}, assumed {A_Q_CORRELATION}"
        );
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
            let q = params.perihelion_typed().get::<astronomical_unit>();
            assert!(q > 0.0, "Perihelion must be positive");
            assert!(q < params.a, "Perihelion must be less than a");
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
