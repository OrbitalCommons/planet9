//! End-to-end reproduction pipeline: clone generation → N-body
//! integration (p9-core WHM) → semi-major-axis diffusion → stability
//! classification → bias-aware clustering statistics → bimodal fit.
//!
//! This replaces the previous hard-coded "results": every number in
//! `PipelineResult` is computed.
//!
//! Dynamics: Neptune is integrated directly (the clones' perihelia reach
//! down to ~36 AU where the J2 ring expansion fails); Jupiter, Saturn and
//! Uranus enter through the orbit-averaged `ExtraForce::J2Jsu` field;
//! Planet Nine can be added as an optional massive body. The default
//! configuration is a reduced-scale smoke run (5 clones × 1 Myr); the
//! paper-scale 25 clones × 4 Gyr run uses the same code path and is
//! exercised by an `#[ignore]`d test.

use p9_core::analysis::circular::{mean_resultant_length, rayleigh_p_value};
use p9_core::constants::{GM_SUN, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets::neptune_j2000;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{cartesian_to_elements, P9Params, SimConfig, StateVector};
use p9_core::units::{days, julian_year, radians, Angle, Time};
use rand::SeedableRng;
use serde::{Deserialize, Serialize};

use crate::clone_generation::{distant_tno_sample, generate_clones, ObservedTno};
use crate::clustering::{
    bias_resampled_p_value, default_selection_weight, fit_bimodal_von_mises, BimodalVonMises,
};
use crate::diffusion::{aggregate_clone_diffusion, measure_diffusion};
use crate::stability::{classify, classify_batch, ClassificationSummary, StabilityClass};

/// Pipeline configuration.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PipelineConfig {
    /// Clones per TNO (paper: 25)
    pub n_clones: usize,
    /// Total integration time (yr; paper: 4e9)
    pub t_total_yr: f64,
    /// Timestep (days; must resolve Neptune, ≤ ~3000)
    pub dt_days: f64,
    /// Snapshot cadence for the a(t) series (yr)
    pub snapshot_interval_yr: f64,
    /// RNG seed (threaded through clone generation and the MC null)
    pub seed: u64,
    /// Monte Carlo samples for the bias-resampling null
    pub n_mc: usize,
    /// Optional Planet Nine
    pub p9: Option<P9Params>,
}

impl PipelineConfig {
    /// Reduced-scale smoke configuration: 5 clones × 1 Myr.
    pub fn smoke() -> Self {
        Self {
            n_clones: 5,
            t_total_yr: 1.0e6,
            dt_days: 2000.0,
            snapshot_interval_yr: 2000.0,
            seed: 20250610,
            n_mc: 200,
            p9: None,
        }
    }

    /// Paper-scale configuration: 25 clones × 4 Gyr (hours of CPU; used
    /// behind `#[ignore]`).
    pub fn paper() -> Self {
        Self {
            n_clones: 25,
            t_total_yr: 4.0e9,
            dt_days: 2000.0,
            snapshot_interval_yr: 1.0e6,
            seed: 20250610,
            n_mc: 10_000,
            p9: None,
        }
    }

    /// Total integration time as a dimension-checked [`Time`].
    pub fn total_time(&self) -> Time {
        julian_year() * self.t_total_yr
    }

    /// Integration timestep as a dimension-checked [`Time`].
    pub fn timestep(&self) -> Time {
        days(self.dt_days)
    }

    /// Snapshot cadence as a dimension-checked [`Time`].
    pub fn snapshot_interval(&self) -> Time {
        julian_year() * self.snapshot_interval_yr
    }
}

/// Per-TNO pipeline output.
#[derive(Debug, Clone, Serialize)]
pub struct TnoResult {
    pub name: &'static str,
    /// Observed longitude of perihelion (rad)
    pub varpi: f64,
    /// Clone-averaged measured diffusion coefficient (AU²/yr)
    pub d_mean: f64,
    /// Clone scatter of the diffusion coefficient
    pub d_std: f64,
    /// Analytical prediction at the observed perihelion
    pub d_analytical: f64,
    /// Stability class from the measured diffusion
    pub class: StabilityClass,
    /// Clones lost to the integration boundaries
    pub n_escaped: usize,
}

impl TnoResult {
    /// Observed longitude of perihelion ϖ as a dimension-checked [`Angle`].
    pub fn varpi_typed(&self) -> Angle {
        radians(self.varpi)
    }
}

/// Computed pipeline results (no hard-coded paper numbers).
#[derive(Debug, Clone, Serialize)]
pub struct PipelineResult {
    pub per_tno: Vec<TnoResult>,
    pub summary: ClassificationSummary,
    /// Mean resultant length of the stable+metastable ϖ sample
    pub r_bar: f64,
    /// Rayleigh test p-value (uniform null, small-n corrected)
    pub rayleigh_p: f64,
    /// Seeded bias-resampling MC p-value (survey-weighted null)
    pub bias_mc_p: f64,
    /// EM fit of a 2-component von Mises mixture to the ϖ of the
    /// unstable group (falls back to the full sample when the group has
    /// fewer than 4 members)
    pub bimodal_fit: BimodalVonMises,
    /// Seed used (recorded for reproducibility)
    pub seed: u64,
}

/// Integrate all clones of all TNOs as test particles in one WHM run and
/// return the per-clone a(t) series (years between snapshots =
/// `config.snapshot_interval_yr`). A series ends when its particle is
/// removed at the integration boundaries.
fn integrate_clone_series(
    tnos: &[ObservedTno],
    config: &PipelineConfig,
    rng: &mut rand::rngs::StdRng,
) -> (Vec<Vec<f64>>, Vec<usize>) {
    // Clone generation (seeded).
    let mut owners = Vec::new(); // particle index → TNO index
    let mut particles: Vec<StateVector> = Vec::new();
    for (t_idx, tno) in tnos.iter().enumerate() {
        for clone in generate_clones(tno, config.n_clones, rng) {
            owners.push(t_idx);
            particles.push(clone.elements().to_state_vector(GM_SUN));
        }
    }

    // Massive bodies: Neptune direct (+ optional P9); JSU via J2 field.
    let mut bodies = vec![neptune_j2000()];
    if let Some(p9) = &config.p9 {
        bodies.push(p9.to_body());
    }
    let integrator = WhmIntegrator::with_extra_forces(vec![ExtraForce::J2Jsu]);

    let t_total_days = config.t_total_yr * YEAR_DAYS;
    let sim = SimConfig {
        dt: config.dt_days,
        t_start: 0.0,
        t_end: t_total_days,
        removal_inner_au: 5.0,
        removal_outer_au: 20_000.0,
        snapshot_interval_days: config.snapshot_interval_yr * YEAR_DAYS,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    let n_steps = (t_total_days / config.dt_days).ceil() as usize;
    let snap_every =
        ((config.snapshot_interval_yr * YEAR_DAYS / config.dt_days).round() as usize).max(1);

    let mut active = vec![true; particles.len()];
    let mut series: Vec<Vec<f64>> = particles
        .iter()
        .map(|p| vec![cartesian_to_elements(p, GM_SUN).a])
        .collect();

    for step in 1..=n_steps {
        integrator.step(&mut bodies, &mut particles, &mut active, sim.dt, &sim);
        if step % snap_every == 0 {
            for (k, particle) in particles.iter().enumerate() {
                if active[k] {
                    series[k].push(cartesian_to_elements(particle, GM_SUN).a);
                }
            }
        }
    }

    (series, owners)
}

/// Run the full pipeline (see module docs).
pub fn run_pipeline(config: &PipelineConfig) -> PipelineResult {
    let tnos = distant_tno_sample();
    let mut rng = rand::rngs::StdRng::seed_from_u64(config.seed);

    let (series, owners) = integrate_clone_series(&tnos, config, &mut rng);

    // Per-clone diffusion → per-TNO aggregation → classification.
    let mut per_tno = Vec::with_capacity(tnos.len());
    let mut results_for_summary = Vec::with_capacity(tnos.len());
    for (t_idx, tno) in tnos.iter().enumerate() {
        let mut d_values = Vec::new();
        let mut n_escaped = 0usize;
        for (k, s) in series.iter().enumerate() {
            if owners[k] != t_idx {
                continue;
            }
            if s.len() < 4 {
                n_escaped += 1;
                continue;
            }
            d_values.push(measure_diffusion(s, config.snapshot_interval_yr, 5));
        }
        let agg = aggregate_clone_diffusion(tno.name, &d_values, tno.perihelion());
        let class = classify(&agg);
        per_tno.push(TnoResult {
            name: tno.name,
            varpi: tno.varpi(),
            d_mean: agg.d_mean,
            d_std: agg.d_std,
            d_analytical: agg.d_analytical,
            class,
            n_escaped,
        });
        results_for_summary.push(agg);
    }
    let summary = classify_batch(&results_for_summary);

    // Clustering of the observed ϖ, grouped by computed stability class.
    let group = |classes: &[StabilityClass]| -> Vec<f64> {
        per_tno
            .iter()
            .filter(|t| classes.contains(&t.class))
            .map(|t| t.varpi)
            .collect()
    };
    let mut clustered_varpis = group(&[StabilityClass::Stable, StabilityClass::Metastable]);
    if clustered_varpis.len() < 3 {
        clustered_varpis = per_tno.iter().map(|t| t.varpi).collect();
    }
    let r_bar = mean_resultant_length(&clustered_varpis);
    let rayleigh_p = rayleigh_p_value(&clustered_varpis);
    let bias_mc_p = bias_resampled_p_value(
        &clustered_varpis,
        &default_selection_weight,
        config.n_mc,
        &mut rng,
    );

    // Bimodal EM fit on the unstable group (or the full sample when the
    // group is too small to fit).
    let mut unstable_varpis = group(&[StabilityClass::Unstable]);
    if unstable_varpis.len() < 4 {
        unstable_varpis = per_tno.iter().map(|t| t.varpi).collect();
    }
    let bimodal_fit = fit_bimodal_von_mises(&unstable_varpis, 100);

    PipelineResult {
        per_tno,
        summary,
        r_bar,
        rayleigh_p,
        bias_mc_p,
        bimodal_fit,
        seed: config.seed,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::radian;
    use uom::si::time::day;

    #[test]
    fn config_typed_times_match_f64_sources() {
        let c = PipelineConfig::smoke();
        assert_relative_eq!(c.timestep().get::<day>(), c.dt_days, epsilon = 1e-9);
        // Years convert through the Julian-year length.
        let yr_d = julian_year().get::<day>();
        assert_relative_eq!(
            c.total_time().get::<day>(),
            c.t_total_yr * yr_d,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            c.snapshot_interval().get::<day>(),
            c.snapshot_interval_yr * yr_d,
            max_relative = 1e-12
        );
    }

    #[test]
    fn tno_result_varpi_typed_matches_f64() {
        let r = TnoResult {
            name: "x",
            varpi: 1.234,
            d_mean: 0.0,
            d_std: 0.0,
            d_analytical: 0.0,
            class: StabilityClass::Stable,
            n_escaped: 0,
        };
        assert_relative_eq!(r.varpi_typed().get::<radian>(), r.varpi, epsilon = 1e-12);
    }

    /// Reduced-scale smoke run: real integration, real statistics.
    #[test]
    fn smoke_pipeline_runs_and_is_coherent() {
        let mut config = PipelineConfig::smoke();
        config.n_clones = 3;
        config.t_total_yr = 2.0e5;
        config.snapshot_interval_yr = 1000.0;
        let result = run_pipeline(&config);

        // One result per vetted TNO.
        assert_eq!(result.per_tno.len(), 10);

        // Diffusion: finite, non-negative; analytical predictions present.
        for t in &result.per_tno {
            assert!(t.d_mean.is_finite() && t.d_mean >= 0.0, "{}", t.name);
            assert!(t.d_analytical > 0.0);
            assert!(t.n_escaped <= config.n_clones);
        }

        // The classification summary covers the sample.
        let total =
            result.summary.n_stable + result.summary.n_metastable + result.summary.n_unstable;
        assert_eq!(total, 10);

        // Detached ETNOs (q ≳ 36 AU) must not be wildly diffusive over
        // 0.2 Myr: most of the sample classifies as stable/metastable.
        assert!(
            result.summary.n_stable + result.summary.n_metastable >= 5,
            "summary = {:?}",
            result.summary
        );

        // Statistics are well-formed and computed (not hard-coded).
        assert!(result.r_bar > 0.0 && result.r_bar <= 1.0);
        assert!(result.rayleigh_p > 0.0 && result.rayleigh_p <= 1.0);
        assert!(result.bias_mc_p > 0.0 && result.bias_mc_p <= 1.0);
        assert!(result.bimodal_fit.w > 0.0 && result.bimodal_fit.w < 1.0);
        assert!(result.bimodal_fit.comp1.kappa >= 0.0);

        // The vetted ETNO sample is ϖ-clustered: R̄ should be sizable.
        assert!(result.r_bar > 0.3, "R̄ = {}", result.r_bar);
    }

    /// Determinism: the same seed reproduces the same statistics.
    #[test]
    fn pipeline_is_seeded_and_reproducible() {
        let mut config = PipelineConfig::smoke();
        config.n_clones = 2;
        config.t_total_yr = 5.0e4;
        config.snapshot_interval_yr = 1000.0;
        config.n_mc = 50;
        let r1 = run_pipeline(&config);
        let r2 = run_pipeline(&config);
        assert_eq!(r1.seed, r2.seed);
        assert!((r1.r_bar - r2.r_bar).abs() < 1e-15);
        assert!((r1.bias_mc_p - r2.bias_mc_p).abs() < 1e-15);
        for (a, b) in r1.per_tno.iter().zip(r2.per_tno.iter()) {
            assert!((a.d_mean - b.d_mean).abs() < 1e-15, "{}", a.name);
        }
    }

    /// Paper-scale configuration (25 clones × 4 Gyr) — same code path,
    /// hours of CPU.
    #[test]
    #[ignore = "paper scale: 25 clones x 4 Gyr integration (hours)"]
    fn paper_scale_pipeline() {
        let result = run_pipeline(&PipelineConfig::paper());
        assert_eq!(result.per_tno.len(), 10);
        assert!(result.rayleigh_p > 0.0 && result.rayleigh_p <= 1.0);
    }
}
