//! Orbital-angle grouping significance of the ETNO sample.
//!
//! de la Fuente Marcos et al. (2016) argue the ETNOs are dynamically grouped:
//! their longitude of the ascending node Ω and longitude of perihelion ϖ are
//! not uniform on the circle. We measure this directly with the workspace's
//! single circular-statistics implementation (`p9_core::analysis::circular`):
//! the mean resultant length R̄ ∈ [0, 1] and the small-n-corrected Rayleigh
//! p-value. We then run a seeded uniform-angle control to confirm the observed
//! grouping is far stronger than a population whose angles are drawn at random.

use rand::SeedableRng;
use rand_distr::{Distribution, Uniform};

use p9_core::analysis::circular::{circular_mean, mean_resultant_length, rayleigh_p_value};
use p9_core::constants::{RAD2DEG, TWO_PI};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::units::{radians, Angle as UnitAngle};

/// Which orbital angle to test for grouping.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Angle {
    /// Longitude of the ascending node Ω.
    Node,
    /// Longitude of perihelion ϖ = ω + Ω.
    LongitudeOfPerihelion,
}

/// Grouping summary for one orbital angle of the observed sample.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct GroupingResult {
    /// Mean resultant length R̄ ∈ [0, 1]; 1 = perfectly aligned.
    pub r_bar: f64,
    /// Mean direction of the angle (radians, [0, 2π)).
    pub mean_direction: f64,
    /// Mean direction in degrees.
    pub mean_direction_deg: f64,
    /// Rayleigh test p-value against circular uniformity (small-n corrected).
    pub rayleigh_p: f64,
    /// Number of objects.
    pub n: usize,
}

impl GroupingResult {
    /// Mean direction of the grouped angle as a typed [`uom`] angle.
    pub fn mean_direction_typed(&self) -> UnitAngle {
        radians(self.mean_direction)
    }
}

/// Collect one orbital angle (radians) for every object in the vetted sample.
pub fn observed_angles(angle: Angle) -> Vec<f64> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|e| match angle {
            Angle::Node => (e.omega_big_deg * p9_core::constants::DEG2RAD).rem_euclid(TWO_PI),
            Angle::LongitudeOfPerihelion => e.longitude_of_perihelion(),
        })
        .collect()
}

/// Grouping statistics for one orbital angle of an explicit angle set.
pub fn grouping_of(angles: &[f64]) -> GroupingResult {
    let mean = circular_mean(angles).unwrap_or(0.0);
    GroupingResult {
        r_bar: mean_resultant_length(angles),
        mean_direction: mean,
        mean_direction_deg: (mean * RAD2DEG).rem_euclid(360.0),
        rayleigh_p: rayleigh_p_value(angles),
        n: angles.len(),
    }
}

/// Grouping statistics for one orbital angle of the observed ETNO sample.
pub fn observed_grouping(angle: Angle) -> GroupingResult {
    grouping_of(&observed_angles(angle))
}

/// Monte Carlo p-value: fraction of seeded uniform-angle controls (each of `n`
/// angles drawn uniformly on the circle) whose mean resultant length is at
/// least as large as `observed_r_bar`. This is the empirical, distribution-free
/// counterpart of the analytic Rayleigh p — it does not assume the asymptotic
/// series and so is the honest small-n significance of the observed grouping.
pub fn uniform_control_exceedance(
    observed_r_bar: f64,
    n: usize,
    n_iterations: usize,
    seed: u64,
) -> f64 {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let dist = Uniform::new(0.0, TWO_PI);
    let mut n_exceed = 0usize;
    for _ in 0..n_iterations {
        let angles: Vec<f64> = (0..n).map(|_| dist.sample(&mut rng)).collect();
        if mean_resultant_length(&angles) >= observed_r_bar {
            n_exceed += 1;
        }
    }
    n_exceed as f64 / n_iterations as f64
}

/// Mean resultant length of a single seeded uniform-angle control of `n`
/// angles — a baseline for "what a random population looks like".
pub fn uniform_control_r_bar(n: usize, seed: u64) -> f64 {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let dist = Uniform::new(0.0, TWO_PI);
    let angles: Vec<f64> = (0..n).map(|_| dist.sample(&mut rng)).collect();
    mean_resultant_length(&angles)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_observed_angles_match_sample_size() {
        assert_eq!(observed_angles(Angle::Node).len(), BROWN_2017_SAMPLE.len());
        assert_eq!(
            observed_angles(Angle::LongitudeOfPerihelion).len(),
            BROWN_2017_SAMPLE.len()
        );
    }

    /// The defining feature of the sample: both Ω and ϖ are grouped, with a
    /// Rayleigh p below a stated threshold and R̄ well above what a random
    /// population gives. Pinned threshold: p < 0.05 for ϖ, p < 0.2 for Ω
    /// (the node is the weaker of the two groupings in this 10-object set).
    #[test]
    fn test_observed_angles_cluster_significantly() {
        let varpi = observed_grouping(Angle::LongitudeOfPerihelion);
        assert!(
            varpi.rayleigh_p < 0.05,
            "ϖ Rayleigh p = {:.4} (R̄ = {:.3}, mean {:.1}°)",
            varpi.rayleigh_p,
            varpi.r_bar,
            varpi.mean_direction_deg
        );
        assert!(varpi.r_bar > 0.3, "ϖ R̄ = {:.3}", varpi.r_bar);

        let node = observed_grouping(Angle::Node);
        assert!(
            node.rayleigh_p < 0.2,
            "Ω Rayleigh p = {:.4} (R̄ = {:.3})",
            node.rayleigh_p,
            node.r_bar
        );
    }

    /// The observed grouping must beat a seeded uniform control: the MC
    /// exceedance probability (fraction of random samples at least as grouped)
    /// is small, and the observed R̄ exceeds a typical random R̄.
    #[test]
    fn test_observed_beats_uniform_control() {
        let n = BROWN_2017_SAMPLE.len();
        for angle in [Angle::Node, Angle::LongitudeOfPerihelion] {
            let obs = observed_grouping(angle);
            let p_mc = uniform_control_exceedance(obs.r_bar, n, 50_000, 1604);
            assert!(
                p_mc < 0.2,
                "{:?}: MC exceedance {:.4} (R̄ = {:.3})",
                angle,
                p_mc,
                obs.r_bar
            );
            // A single random control of the same size is, with high
            // probability, less grouped than the observed sample. (This
            // previously OR'd in `p_mc < 0.2`, already asserted above, so
            // the comparison could never fire.)
            let random_r = uniform_control_r_bar(n, 9999);
            assert!(
                obs.r_bar > random_r - 0.05,
                "{:?}: obs R̄ {:.3} vs random {:.3}",
                angle,
                obs.r_bar,
                random_r
            );
        }
    }

    /// The MC exceedance probability and the analytic Rayleigh p must agree to
    /// the same order of magnitude for ϖ (both are testing uniformity); this
    /// cross-validates the random control against the closed-form test.
    #[test]
    fn test_mc_matches_rayleigh_order_of_magnitude() {
        let obs = observed_grouping(Angle::LongitudeOfPerihelion);
        let p_mc = uniform_control_exceedance(obs.r_bar, obs.n, 200_000, 7);
        // Within a factor of ~3 of the analytic p.
        let ratio = (p_mc.max(1e-6)) / obs.rayleigh_p;
        assert!(
            (0.33..3.0).contains(&ratio),
            "MC p {:.5} vs Rayleigh p {:.5} (ratio {:.2})",
            p_mc,
            obs.rayleigh_p,
            ratio
        );
    }

    /// The typed mean-direction accessor matches its f64 (radian) source.
    #[test]
    fn test_mean_direction_typed_matches_f64() {
        use approx::assert_relative_eq;
        let g = observed_grouping(Angle::LongitudeOfPerihelion);
        assert_relative_eq!(
            (g.mean_direction_typed() / radians(1.0)).value,
            g.mean_direction,
            max_relative = 1e-12
        );
    }

    /// Seeded determinism.
    #[test]
    fn test_control_is_seeded() {
        let a = uniform_control_exceedance(0.5, 10, 5000, 42);
        let b = uniform_control_exceedance(0.5, 10, 5000, 42);
        assert_eq!(a, b);
    }

    /// A truly uniform population is not flagged as grouped: a random control's
    /// own R̄ is exceeded by random samples roughly half the time or more.
    #[test]
    fn test_uniform_population_not_grouped() {
        let r = uniform_control_r_bar(10, 123);
        let p = uniform_control_exceedance(r, 10, 50_000, 456);
        assert!(p > 0.1, "uniform control exceedance {p:.3} (R̄ {r:.3})");
    }
}
