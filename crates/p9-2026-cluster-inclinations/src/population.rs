//! Initial test-particle populations (paper Sections 2.1 and 2.3).
//!
//! * **Cluster-influenced**: the paper draws 10⁴ particles with
//!   a ∼ 100–5000 AU and q > 30 AU from the τ = 300 Myr snapshot of the
//!   `cluster_2` run of Nesvorný et al. (2023) — a population isotropically
//!   scattered by the strongest cluster perturbations the cold belt allows.
//!   That snapshot is not reproducible here, so the population is modelled
//!   as a mixture of an isotropic (uniform in cos i) component and a cold
//!   half-Gaussian core (σ = 15°). The isotropic fraction
//!   [`ISOTROPIC_FRACTION`] is the one free parameter; it is set so the
//!   *initial* fitted width lands in the paper's cluster-influenced band
//!   (26–27.5°, Table 1) — the paper's Figure 2 shows that band is produced
//!   by a distribution that is flat above ~17°, i.e. dominated by the
//!   isotropic part. Semi-major axes are log-uniform on the paper's range
//!   and perihelia uniform on (30, 100) AU.
//! * **Cluster-free**: the Batygin et al. (2019) population — 1000
//!   particles, half-Gaussian inclinations with σ = 15°, a ∈ (100, 800) AU,
//!   q ∈ (30, 100) AU — which is exactly p9-core's scattered-disk generator.
//!
//! Note on measures: a half-Gaussian of σ = 15° in `i` fits the paper's
//! `sin i · exp(−i²/2w²)` form with `w ≈ 10°`, not 15° (the model's `sin i`
//! factor expects more high-i orbits than a half-Gaussian has). The paper's
//! cluster-free runs therefore *start* at w ≈ 10° in its own measure, and
//! the 16–19.5° they end at is Planet Nine's forced-plane broadening of the
//! detached population (Anderson & Kaib 2021), not a preserved 15°.

use p9_core::constants::{DEG2RAD, GM_SUN, TWO_PI};
use p9_core::initial_conditions::scattered_disk::{generate_scattered_disk, ScatteredDiskConfig};
use p9_core::types::{elements_to_cartesian, OrbitalElements, StateVector};
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Which birth environment the initial conditions represent.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum BirthEnvironment {
    ClusterInfluenced,
    ClusterFree,
}

/// Cold-core inclination dispersion shared by both populations (degrees).
pub const COLD_SIGMA_DEG: f64 = 15.0;

/// Fraction of the cluster-influenced population drawn isotropically.
pub const ISOTROPIC_FRACTION: f64 = 0.83;

/// Inclination ceiling of the isotropic component (degrees). The width
/// statistic only sees i < 40°, so orbits generated far above the cut would
/// cost integration time without entering the statistic; the ceiling is
/// kept 10° above the cut so Planet Nine's inclination pumping across the
/// boundary is still represented in both directions.
pub const ISOTROPIC_I_MAX_DEG: f64 = 50.0;

/// Cluster-influenced ranges (paper Section 2.1).
pub const CLUSTER_A_RANGE: (f64, f64) = (100.0, 5000.0);
pub const CLUSTER_Q_RANGE: (f64, f64) = (30.0, 100.0);

/// Cluster-free ranges (paper Section 2.3, Batygin et al. 2019).
pub const FREE_A_RANGE: (f64, f64) = (100.0, 800.0);
pub const FREE_Q_RANGE: (f64, f64) = (30.0, 100.0);

impl BirthEnvironment {
    /// Upper perihelion bound of the environment's default population (AU).
    pub fn default_q_max(self) -> f64 {
        match self {
            BirthEnvironment::ClusterInfluenced => CLUSTER_Q_RANGE.1,
            BirthEnvironment::ClusterFree => FREE_Q_RANGE.1,
        }
    }
}

/// Generate `n` heliocentric state vectors for the given environment, with
/// perihelia drawn up to `q_max` (AU).
pub fn generate<R: Rng>(
    env: BirthEnvironment,
    n: usize,
    q_max: f64,
    rng: &mut R,
) -> Vec<StateVector> {
    match env {
        BirthEnvironment::ClusterFree => generate_scattered_disk(
            &ScatteredDiskConfig {
                a_min: FREE_A_RANGE.0,
                a_max: FREE_A_RANGE.1,
                q_min: FREE_Q_RANGE.0,
                q_max,
                sigma_i: COLD_SIGMA_DEG * DEG2RAD,
                n_particles: n,
            },
            rng,
        ),
        BirthEnvironment::ClusterInfluenced => (0..n)
            .map(|_| {
                elements_to_cartesian(
                    &cluster_influenced_orbit(rng, ISOTROPIC_FRACTION, q_max),
                    GM_SUN,
                )
            })
            .collect(),
    }
}

/// One cluster-influenced orbit with the given isotropic fraction and
/// perihelion ceiling.
pub fn cluster_influenced_orbit<R: Rng>(
    rng: &mut R,
    isotropic_fraction: f64,
    q_max: f64,
) -> OrbitalElements {
    let cold = Normal::new(0.0, COLD_SIGMA_DEG * DEG2RAD).unwrap();
    let (ln_lo, ln_hi) = (CLUSTER_A_RANGE.0.ln(), CLUSTER_A_RANGE.1.ln());
    let (a, q) = loop {
        let a = rng.gen_range(ln_lo..ln_hi).exp();
        let q = rng.gen_range(CLUSTER_Q_RANGE.0..q_max);
        if q < a {
            break (a, q);
        }
    };
    let i = if rng.gen::<f64>() < isotropic_fraction {
        let cos_max = (ISOTROPIC_I_MAX_DEG * DEG2RAD).cos();
        rng.gen_range(cos_max..1.0f64).acos()
    } else {
        cold.sample(rng).abs()
    };
    OrbitalElements {
        a,
        e: 1.0 - q / a,
        i,
        omega_big: rng.gen_range(0.0..TWO_PI),
        omega: rng.gen_range(0.0..TWO_PI),
        mean_anomaly: rng.gen_range(0.0..TWO_PI),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::width::{population_width, SelectionCuts};
    use p9_core::types::cartesian_to_elements;
    use rand::SeedableRng;

    fn widths(env: BirthEnvironment) -> f64 {
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        let states = generate(env, 20_000, env.default_q_max(), &mut rng);
        let elements: Vec<_> = states
            .iter()
            .map(|s| cartesian_to_elements(s, GM_SUN))
            .collect();
        population_width(&elements, &SelectionCuts::width_sample()).0
    }

    #[test]
    fn cluster_influenced_starts_in_the_papers_band() {
        let w = widths(BirthEnvironment::ClusterInfluenced);
        assert!(
            (25.5..=28.5).contains(&w),
            "initial cluster-influenced w = {w}"
        );
    }

    #[test]
    fn cluster_free_half_gaussian_fits_to_ten_degrees() {
        let w = widths(BirthEnvironment::ClusterFree);
        assert!((9.0..=11.5).contains(&w), "initial cluster-free w = {w}");
    }
}
