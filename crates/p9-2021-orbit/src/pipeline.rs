//! End-to-end inference pipeline: simulation grid (Model 1) -> conditional
//! KDE likelihood (Model 2) -> GP emulator (Model 3) -> ensemble MCMC
//! (Model 4), wired as in Brown & Batygin (2021).
//!
//! Prior: uniform over the grid extent (m9, a9, e9, i9) — the paper's
//! primary prior. The paper's alternative *cluster-scattering* prior (a
//! Frechet distribution over m9/a9 from a planet-scattering population
//! synthesis, their Section 5) is **not implemented**: it requires the
//! scattering-simulation suite of Batygin & Brown (2021b), which is out of
//! scope here.
//!
//! Fidelity: at `GridScale::Smoke`/`Reduced` the pipeline validates the
//! *machinery* (deterministic grid -> finite likelihood surface -> healthy
//! GP cross-validation -> converged MCMC). It does NOT reproduce the
//! published posterior (m9 = 6.2, a9 = 380 AU, ...): that requires the
//! `Paper`-scale grid (121 points x 16,800 particles x 4 Gyr, estimated
//! ~3,000 CPU-days — see REPRODUCTION_NOTES.md). The `#[ignore]`d
//! paper-scale entry point is `run_pipeline(PipelineConfig::paper())`.

use serde::{Deserialize, Serialize};

use p9_core::data::etno::BROWN_2017_SAMPLE;

use crate::gp::{GpFitSummary, Matern32Gp};
use crate::kde::{profiled_log_likelihood, BinSpec, GridPointKde};
use crate::mcmc::{init_walkers_in_box, quantile, EnsembleSampler, PaperMcmcSettings};
use crate::sim_grid::{
    run_grid, GridPoint, GridScale, SimGridConfig, A9_RANGE, E9_RANGE, I9_RANGE_DEG, M9_RANGE,
};

/// Pipeline configuration. The grid scale drives the MCMC budget.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    pub grid: SimGridConfig,
    /// Orientation-profiling scan resolution (Model 2; n x n over
    /// (varpi9, Omega9)).
    pub kde_scan: usize,
    pub n_walkers: usize,
    pub n_steps: usize,
    pub burn_in: usize,
    pub thin: usize,
    /// Base seed for the MCMC stage (grid seeding lives in `grid.seed`).
    pub mcmc_seed: u64,
}

impl PipelineConfig {
    /// Default-test budget: 4-point grid, 12 x 12 orientation scan,
    /// 16 walkers x 400 steps.
    pub fn smoke() -> Self {
        Self {
            grid: SimGridConfig::smoke(),
            kde_scan: 12,
            n_walkers: 16,
            n_steps: 400,
            burn_in: 100,
            thin: 5,
            mcmc_seed: 777,
        }
    }

    /// Reduced scale (release mode, behind `#[ignore]`): 12-point grid,
    /// 24 x 24 scan, 32 walkers x 2000 steps.
    pub fn reduced() -> Self {
        Self {
            grid: SimGridConfig::reduced(),
            kde_scan: 24,
            n_walkers: 32,
            n_steps: 2000,
            burn_in: 500,
            thin: 10,
            mcmc_seed: 777,
        }
    }

    /// Paper-scale entry point (121-point surrogate grid, 4 Gyr, the
    /// paper's MCMC settings). Estimated CPU-days; `#[ignore]`d only.
    pub fn paper() -> Self {
        let mcmc = PaperMcmcSettings::paper();
        Self {
            grid: SimGridConfig::paper(),
            kde_scan: 24,
            n_walkers: mcmc.n_walkers,
            n_steps: mcmc.n_steps,
            burn_in: mcmc.burn_in,
            thin: mcmc.thin,
            mcmc_seed: 777,
        }
    }

    /// The uniform-prior box over (m9, a9 [AU], e9, i9 [deg]): the grid
    /// extent, matching the paper's uniform prior over its grid.
    pub fn prior_bounds() -> [(f64, f64); 4] {
        [M9_RANGE, A9_RANGE, E9_RANGE, I9_RANGE_DEG]
    }
}

/// Posterior quantiles for one parameter.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ParamQuantiles {
    pub q025: f64,
    pub median: f64,
    pub q975: f64,
}

/// Full pipeline output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineResult {
    pub scale: GridScale,
    pub grid_seed: u64,
    pub mcmc_seed: u64,
    pub grid_points: Vec<GridPoint>,
    /// Profiled ETNO log-likelihood per grid point (the GP targets).
    pub log_likelihoods: Vec<f64>,
    pub gp: GpFitSummary,
    /// Flattened post-burn, thinned posterior samples of
    /// (m9, a9 [AU], e9, i9 [deg]).
    pub samples: Vec<Vec<f64>>,
    pub acceptance_rate: f64,
    pub r_hat: Vec<f64>,
    pub ess: Vec<f64>,
    /// Quantiles for (m9, a9, e9, i9).
    pub quantiles: [ParamQuantiles; 4],
}

impl PipelineResult {
    /// Fraction of posterior samples inside the prior box (must be 1.0 for
    /// a healthy run: the prior is -inf outside).
    pub fn fraction_in_prior(&self) -> f64 {
        let bounds = PipelineConfig::prior_bounds();
        let inside = self
            .samples
            .iter()
            .filter(|s| {
                s.iter()
                    .zip(bounds.iter())
                    .all(|(&v, &(lo, hi))| v >= lo && v <= hi)
            })
            .count();
        inside as f64 / self.samples.len().max(1) as f64
    }
}

/// Run the full grid -> KDE -> GP -> MCMC pipeline.
pub fn run_pipeline(config: &PipelineConfig) -> PipelineResult {
    // Model 1: the simulation grid.
    let grid_results = run_grid(&config.grid);
    let grid_points: Vec<GridPoint> = grid_results.iter().map(|r| r.point).collect();

    // Model 2: profiled ETNO log-likelihood per grid point.
    let log_likelihoods: Vec<f64> = grid_results
        .iter()
        .map(|r| {
            let kde = GridPointKde::from_samples(&r.samples, BinSpec::default());
            let (ll, _, _) = profiled_log_likelihood(&kde, &BROWN_2017_SAMPLE, config.kde_scan);
            ll
        })
        .collect();

    // Model 3: GP emulator over normalized grid coordinates.
    let inputs: Vec<[f64; 4]> = grid_points.iter().map(|p| p.normalized()).collect();
    let gp = Matern32Gp::fit(&inputs, &log_likelihoods);
    let gp_summary = gp.summary();

    // Model 4: ensemble MCMC over the uniform prior box x GP likelihood.
    let bounds = PipelineConfig::prior_bounds();
    let log_post = |theta: &[f64]| -> f64 {
        let inside = theta
            .iter()
            .zip(bounds.iter())
            .all(|(&v, &(lo, hi))| v >= lo && v <= hi);
        if !inside {
            return f64::NEG_INFINITY;
        }
        let x = [
            (theta[0] - M9_RANGE.0) / (M9_RANGE.1 - M9_RANGE.0),
            (theta[1] - A9_RANGE.0) / (A9_RANGE.1 - A9_RANGE.0),
            (theta[2] - E9_RANGE.0) / (E9_RANGE.1 - E9_RANGE.0),
            (theta[3] - I9_RANGE_DEG.0) / (I9_RANGE_DEG.1 - I9_RANGE_DEG.0),
        ];
        gp.predict(&x).0
    };

    let sampler = EnsembleSampler::new(config.n_walkers);
    let init = init_walkers_in_box(config.n_walkers, &bounds, config.mcmc_seed);
    let chain = sampler.run(log_post, &init, config.n_steps, config.mcmc_seed + 1);

    let r_hat = chain.gelman_rubin(config.burn_in);
    let ess = chain.effective_sample_size(config.burn_in);
    let samples = chain.flat_samples(config.burn_in, config.thin);

    let quantiles = std::array::from_fn(|dim| {
        let mut values: Vec<f64> = samples.iter().map(|s| s[dim]).collect();
        ParamQuantiles {
            q025: quantile(&mut values, 0.025),
            median: quantile(&mut values, 0.5),
            q975: quantile(&mut values, 0.975),
        }
    });

    PipelineResult {
        scale: config.grid.scale,
        grid_seed: config.grid.seed,
        mcmc_seed: config.mcmc_seed,
        grid_points,
        log_likelihoods,
        gp: gp_summary,
        samples,
        acceptance_rate: chain.acceptance_rate,
        r_hat,
        ess,
        quantiles,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn assert_machinery_health(result: &PipelineResult, r_hat_max: f64, ess_min: f64) {
        assert!(result.log_likelihoods.iter().all(|ll| ll.is_finite()));
        assert!(result.gp.loo_relative_error.is_finite());
        assert!(!result.samples.is_empty());
        assert!(
            (result.fraction_in_prior() - 1.0).abs() < 1e-12,
            "posterior mass escaped the prior box"
        );
        assert!(
            result.acceptance_rate > 0.05 && result.acceptance_rate < 0.95,
            "acceptance = {}",
            result.acceptance_rate
        );
        for (dim, &r) in result.r_hat.iter().enumerate() {
            assert!(r < r_hat_max, "R-hat[{dim}] = {r}");
        }
        for (dim, &e) in result.ess.iter().enumerate() {
            assert!(e > ess_min, "ESS[{dim}] = {e}");
        }
        for q in &result.quantiles {
            assert!(q.q025 <= q.median && q.median <= q.q975);
        }
    }

    /// Default-test smoke run of the full pipeline. Validates machinery
    /// only — at this scale the likelihood surface is noise (see module
    /// docs); no posterior values are asserted.
    #[test]
    fn test_pipeline_smoke() {
        let result = run_pipeline(&PipelineConfig::smoke());
        assert_eq!(result.grid_points.len(), 4);
        assert_eq!(result.log_likelihoods.len(), 4);
        // Short chains on a flat-ish surface: loose health bounds.
        assert_machinery_health(&result, 1.5, 50.0);
    }

    #[test]
    fn test_pipeline_smoke_is_deterministic() {
        let r1 = run_pipeline(&PipelineConfig::smoke());
        let r2 = run_pipeline(&PipelineConfig::smoke());
        assert_eq!(r1.log_likelihoods, r2.log_likelihoods);
        assert_eq!(r1.samples, r2.samples);
    }

    /// Reduced-scale end-to-end run (release mode):
    /// `cargo test --release -p p9-2021-orbit -- --ignored`.
    ///
    /// Asserts machinery health (R-hat < 1.2, ESS, posterior in prior, GP
    /// LOO sane) and the loose sanity check that the published medians lie
    /// inside the Reduced posterior's 95% box. It does NOT assert the
    /// published medians as point estimates — full-fidelity reproduction
    /// requires the Paper-scale grid (see REPRODUCTION_NOTES.md).
    #[test]
    #[ignore]
    fn test_pipeline_reduced_end_to_end() {
        let result = run_pipeline(&PipelineConfig::reduced());

        // Run record (REPRODUCTION_NOTES.md documents the measured values):
        eprintln!(
            "reduced pipeline: grid_seed={} mcmc_seed={} acceptance={:.3}",
            result.grid_seed, result.mcmc_seed, result.acceptance_rate
        );
        eprintln!(
            "  log-likelihoods = {:?}\n  GP: ell={} noise={} LOO_rel={:.4}\n  R-hat={:?}\n  ESS={:?}\n  quantiles={:?}",
            result.log_likelihoods,
            result.gp.length_scale,
            result.gp.noise_var,
            result.gp.loo_relative_error,
            result.r_hat,
            result.ess,
            result.quantiles
        );

        assert_eq!(result.grid_points.len(), 12);
        assert_machinery_health(&result, 1.2, 200.0);
        // Reduced 12-point grid: measured GP LOO relative error 0.239
        // (seeds 2021/777, run 2026-06-12; ~68 s wall / 12.9 CPU-min on 64
        // cores). The paper's < 5% applies to its 121-point grid — see
        // REPRODUCTION_NOTES.md. Pinned with headroom for cross-platform
        // FP jitter.
        assert!(
            result.gp.loo_relative_error < 0.35,
            "GP LOO relative error = {}",
            result.gp.loo_relative_error
        );

        // Loose sanity: published medians (Brown & Batygin 2021) inside the
        // Reduced posterior's 95% intervals. At reduced scale the posterior
        // is broad, so this is a weak containment check, not a reproduction.
        let published = [6.2, 380.0, 0.21, 16.0]; // m9, a9, e9 (1 - 300/380), i9
        for (dim, (&pub_val, q)) in published.iter().zip(result.quantiles.iter()).enumerate() {
            assert!(
                pub_val >= q.q025 && pub_val <= q.q975,
                "published value {pub_val} outside 95% interval [{}, {}] for dim {dim}",
                q.q025,
                q.q975
            );
        }
    }

    /// Paper-scale entry point: 121-point grid x 16,800 particles x 4 Gyr
    /// with the paper's MCMC settings. Estimated CPU-days — see
    /// REPRODUCTION_NOTES.md. Expected (not yet verified here): posterior
    /// medians within the published m9 = 6.2 (+2.2/-1.3), a9 = 380
    /// (+140/-80) AU, e9 = 0.3 (+0.1/-0.1), i9 = 16 +/- 5 deg intervals.
    #[test]
    #[ignore]
    fn test_pipeline_paper_scale() {
        let result = run_pipeline(&PipelineConfig::paper());
        assert_eq!(result.grid_points.len(), 121);
        assert_machinery_health(&result, 1.1, 1000.0);
        assert!(result.gp.loo_relative_error < 0.05);
    }
}
