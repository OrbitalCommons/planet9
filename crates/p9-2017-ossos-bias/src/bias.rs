//! The OSSOS VI bias experiment.
//!
//! Generate a large intrinsically isotropic ETNO population (uniform ϖ), run
//! it through the wide-field survey detection model, and compare the
//! longitude-of-perihelion clustering statistics of the parent and the
//! detected subsample.
//!
//! Headline reproduction (Shankman et al. 2017): the parent population is
//! isotropic in ϖ (mean resultant length R̄ ≈ 0, non-significant Rayleigh
//! test), yet the *detected* subsample shows materially elevated R̄ and a
//! small Rayleigh p-value — apparent clustering manufactured entirely by
//! detection bias. The bias-induced R̄ of the detected sample is the
//! reproduced quantity.

use p9_core::analysis::circular::{circular_mean, mean_resultant_length, rayleigh_p_value};
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::detection::{detected_subsample, SurveyModel};
use crate::population::{
    generate_population, longitudes_of_perihelion, PopulationParams, SyntheticEtno,
};

/// Clustering statistics of a longitude-of-perihelion sample.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ClusteringStats {
    /// Number of objects.
    pub n: usize,
    /// Mean resultant length R̄ ∈ [0, 1].
    pub r_bar: f64,
    /// Rayleigh p-value for uniformity.
    pub rayleigh_p: f64,
    /// Mean longitude of perihelion (radians), if defined.
    pub mean_varpi: Option<f64>,
}

impl ClusteringStats {
    fn from_varpis(varpis: &[f64]) -> Self {
        Self {
            n: varpis.len(),
            r_bar: mean_resultant_length(varpis),
            rayleigh_p: rayleigh_p_value(varpis),
            mean_varpi: circular_mean(varpis),
        }
    }
}

/// Result of the bias experiment: parent vs detected ϖ clustering.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BiasExperiment {
    /// Statistics of the intrinsic (parent) ϖ distribution.
    pub intrinsic: ClusteringStats,
    /// Statistics of the detected ϖ distribution.
    pub detected: ClusteringStats,
}

impl BiasExperiment {
    /// Excess clustering R̄(detected) − R̄(intrinsic): the bias-induced mean
    /// resultant length, i.e. the apparent clustering manufactured from an
    /// isotropic input.
    pub fn bias_induced_r_bar(&self) -> f64 {
        self.detected.r_bar - self.intrinsic.r_bar
    }
}

/// Run the experiment with `n` parent objects and the given seed/models.
pub fn run_experiment(
    n: usize,
    seed: u64,
    pop_params: &PopulationParams,
    model: &SurveyModel,
) -> BiasExperiment {
    let mut rng = StdRng::seed_from_u64(seed);
    let pop = generate_population(n, pop_params, &mut rng);

    let intrinsic_varpis = longitudes_of_perihelion(&pop);
    let detected = detected_subsample(&pop, model);
    let detected_varpis = longitudes_of_perihelion(&detected);

    BiasExperiment {
        intrinsic: ClusteringStats::from_varpis(&intrinsic_varpis),
        detected: ClusteringStats::from_varpis(&detected_varpis),
    }
}

/// Run the experiment with default population and survey models.
pub fn run_default(n: usize, seed: u64) -> BiasExperiment {
    run_experiment(
        n,
        seed,
        &PopulationParams::default(),
        &SurveyModel::default(),
    )
}

/// Detected subsample for a default run (convenience for callers that want
/// the objects, not just the statistics).
pub fn detected_default(n: usize, seed: u64) -> Vec<SyntheticEtno> {
    let mut rng = StdRng::seed_from_u64(seed);
    let pop = generate_population(n, &PopulationParams::default(), &mut rng);
    detected_subsample(&pop, &SurveyModel::default())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// HEADLINE: an isotropic parent (R̄ ≈ 0) yields a detected subsample with
    /// materially elevated R̄ and a small Rayleigh p — detection bias
    /// manufactures apparent ϖ clustering. Robust across seeds.
    #[test]
    fn test_bias_manufactures_clustering() {
        for seed in [1u64, 7, 42, 123, 2024] {
            let exp = run_default(60_000, seed);

            // Parent is isotropic: R̄ ≈ 1/√n ≈ 0.004 at n = 60 000.
            assert!(
                exp.intrinsic.r_bar < 0.03,
                "seed {seed}: intrinsic R̄ = {:.4} not ~0",
                exp.intrinsic.r_bar
            );

            // Enough detections to make the statistic meaningful.
            assert!(
                exp.detected.n >= 50,
                "seed {seed}: only {} detections",
                exp.detected.n
            );

            // Detected subsample is materially clustered.
            assert!(
                exp.detected.r_bar > 0.25,
                "seed {seed}: detected R̄ = {:.4} should be materially elevated",
                exp.detected.r_bar
            );
            assert!(
                exp.detected.rayleigh_p < 1e-3,
                "seed {seed}: detected Rayleigh p = {:.3e} should be significant",
                exp.detected.rayleigh_p
            );

            // The bias-induced excess is large and positive.
            assert!(
                exp.bias_induced_r_bar() > 0.2,
                "seed {seed}: bias-induced R̄ = {:.4} too small",
                exp.bias_induced_r_bar()
            );
        }
    }

    /// Reproduced quantity: pin the magnitude of the bias-induced R̄ for the
    /// default model at a fixed seed. With the PUBLISHED OSSOS block layout
    /// (two clumps roughly opposite on the circle, Bannister et al. 2018)
    /// the apparent clustering strength produced from isotropy is R̄ ≈ 0.29 —
    /// weaker than the ~0.5–0.65 the previously invented single-side window
    /// layout manufactured (opposite clumps partially cancel the resultant),
    /// but still Rayleigh-significant at survey sample sizes: selection alone
    /// creates apparent ϖ structure.
    #[test]
    fn test_bias_induced_r_bar_magnitude() {
        let exp = run_default(80_000, 42);
        let r = exp.detected.r_bar;
        assert!(
            r > 0.15 && r < 0.5,
            "detected R̄ = {r:.4} outside expected bias-induced band [0.15, 0.5]"
        );
        // Document the value in the test output.
        eprintln!(
            "OSSOS VI bias: intrinsic R̄ = {:.4} (p = {:.3}), detected R̄ = {:.4} (p = {:.2e}), \
             n_det = {}, bias-induced ΔR̄ = {:.4}",
            exp.intrinsic.r_bar,
            exp.intrinsic.rayleigh_p,
            exp.detected.r_bar,
            exp.detected.rayleigh_p,
            exp.detected.n,
            exp.bias_induced_r_bar()
        );
    }

    /// The detected sample's mean ϖ should sit near one of the covered
    /// discovery-longitude blocks (the bias has a preferred direction), not be
    /// undefined.
    #[test]
    fn test_detected_has_defined_mean_direction() {
        let exp = run_default(60_000, 7);
        assert!(exp.detected.mean_varpi.is_some());
    }
}
