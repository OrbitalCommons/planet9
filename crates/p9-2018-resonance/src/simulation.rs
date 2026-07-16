//! Planar N-body simulation for resonance characterization.
//!
//! Sweeps P9 eccentricity from 0 to 0.7 with a₉ = 600 AU. The giant planets
//! enter as the orbit-averaged J2 quadrupole + monopole boost
//! (`p9_core::forces::ExtraForce::J2Jsu`), actually applied in the WHM kick
//! (a previous version documented the J2 field but integrated bare
//! Sun + P9).
//!
//! Timestep: dt = 2500 d. The resolution criterion is the *perihelion*
//! passage of the most eccentric particles: an orbit with q = 30 AU sweeps
//! perihelion on the timescale of a circular orbit at 30 AU
//! (P ≈ 60,000 d), requiring dt ≲ P/20 ≈ 3000 d. (The previous
//! dt = P_N/8 ≈ 7500 d produced pericenter step artifacts that read as
//! spurious chaos in the resonance census.)
//!
//! Resonance membership is decided from recorded resonant-angle time
//! series (libration vs circulation via `p9_core::analysis::resonance`),
//! not from single-snapshot semi-major-axis proximity.

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets;
use p9_core::initial_conditions::scattered_disk::{generate_scattered_disk, ScatteredDiskConfig};
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::*;
use p9_core::units::{au, radians, Angle, Length};

use crate::resonance_catalog::{
    extended_catalog, is_librating, resonant_angle_from_angles, Resonance,
};

/// Configuration for the resonance characterization simulation.
#[derive(Debug, Clone)]
pub struct ResonanceSimConfig {
    /// Planet Nine semi-major axis (AU) — fixed at 600 for this paper
    pub a_p9: f64,
    /// Planet Nine eccentricity (swept from 0 to 0.7)
    pub e_p9: f64,
    /// Planet Nine mass (Earth masses)
    pub mass_earth: f64,
    /// Number of test particles
    pub n_particles: usize,
    /// Semi-major axis range for test particles (AU)
    pub a_min: f64,
    pub a_max: f64,
    /// Perihelion distance range (AU)
    pub q_min: f64,
    pub q_max: f64,
    /// Total integration time (days)
    pub t_total: f64,
    /// Integration timestep (days)
    pub dt: f64,
    /// Interval between resonant-angle samples (days). Libration periods
    /// are ~10–100 Myr, so ~1 Myr sampling resolves them.
    pub angle_sample_interval: f64,
}

impl ResonanceSimConfig {
    /// Nominal config for a given e₉ value.
    pub fn for_eccentricity(e_p9: f64) -> Self {
        Self {
            a_p9: 600.0,
            e_p9,
            mass_earth: 10.0,
            n_particles: 400,
            a_min: 50.0,
            a_max: 1000.0,
            q_min: 30.0,
            q_max: 50.0,
            t_total: 4.0 * GYR_DAYS,
            // Perihelion-resolution criterion: q = 30 AU passages need
            // dt <= P(30 AU)/20 ~ 3000 d.
            dt: 2500.0,
            angle_sample_interval: 1e6 * YEAR_DAYS,
        }
    }

    /// Quick test config.
    pub fn quick_test(e_p9: f64) -> Self {
        Self {
            n_particles: 20,
            a_min: 150.0,
            a_max: 550.0,
            t_total: 1e5 * YEAR_DAYS,
            angle_sample_interval: 1e3 * YEAR_DAYS,
            ..Self::for_eccentricity(e_p9)
        }
    }

    fn p9_params(&self) -> P9Params {
        P9Params {
            mass_earth: self.mass_earth,
            a: self.a_p9,
            e: self.e_p9,
            i: 0.0,
            omega: 0.0,
            omega_big: PI,
            mean_anomaly: 0.0,
        }
    }
}

/// One sampled point of a particle's angle history.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct AngleSample {
    /// Mean longitude λ = M + ϖ (rad)
    pub lambda: f64,
    /// Longitude of perihelion ϖ (rad)
    pub varpi: f64,
    /// Osculating semi-major axis (AU)
    pub a: f64,
}

impl AngleSample {
    /// Mean longitude λ as a typed [`Angle`] (see [`Self::lambda`]).
    pub fn mean_longitude(&self) -> Angle {
        radians(self.lambda)
    }

    /// Longitude of perihelion ϖ as a typed [`Angle`] (see [`Self::varpi`]).
    pub fn longitude_of_perihelion(&self) -> Angle {
        radians(self.varpi)
    }

    /// Osculating semi-major axis as a typed [`Length`] (see [`Self::a`]).
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }
}

/// Result of one planar run: final elements plus the recorded angle series
/// needed for libration-based resonance classification.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ResonanceRunResult {
    pub t_final: f64,
    /// Final elements per particle slot (None = removed).
    pub final_elements: Vec<Option<OrbitalElements>>,
    /// Per-particle angle series (entries only while active), sampled at
    /// `angle_sample_interval`, index-aligned with `p9_samples`.
    pub particle_series: Vec<Vec<Option<AngleSample>>>,
    /// Planet Nine's angle samples at the same epochs.
    pub p9_samples: Vec<AngleSample>,
    pub active_count: usize,
    pub total_count: usize,
}

/// Run a single planar simulation for a given e₉ with the J2-averaged
/// giant-planet field applied in the kick stage.
pub fn run_planar_simulation(config: &ResonanceSimConfig, seed: u64) -> ResonanceRunResult {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let disk_config = ScatteredDiskConfig {
        a_min: config.a_min,
        a_max: config.a_max,
        q_min: config.q_min,
        q_max: config.q_max,
        sigma_i: 0.001, // Near-planar
        n_particles: config.n_particles,
    };
    let mut particles = generate_scattered_disk(&disk_config, &mut rng);
    let n_actual = particles.len();
    let mut active: Vec<bool> = vec![true; n_actual];

    let p9 = config.p9_params();
    // Neptune is integrated directly: it is the largest single ring term
    // (m·a² ≈ 15,500 M⊕AU² vs Saturn's 8,700) and the seeded q ∈ [30, 50] AU
    // scattered disk is Neptune-coupled — j2_secular's own docs say not to
    // average Neptune for Neptune-crossing perihelia. J+S+U stay in the
    // averaged J2 field below. dt = 2500 d resolves P_N/20 ≈ 3000 d.
    let mut bodies = vec![planets::neptune_j2000(), p9.to_body()];

    let sim_config = SimConfig {
        dt: config.dt,
        t_start: 0.0,
        t_end: config.t_total,
        removal_inner_au: 30.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: config.angle_sample_interval,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    // H1 fix: the giant-planet quadrupole is now genuinely applied.
    let integrator = WhmIntegrator::with_extra_forces(vec![ExtraForce::J2Jsu]);

    let mut particle_series: Vec<Vec<Option<AngleSample>>> = vec![Vec::new(); n_actual];
    let mut p9_samples: Vec<AngleSample> = Vec::new();

    let sample = |state: &StateVector| -> AngleSample {
        let elem = cartesian_to_elements(state, GM_SUN);
        let varpi = elem.omega + elem.omega_big;
        AngleSample {
            lambda: (elem.mean_anomaly + varpi).rem_euclid(TWO_PI),
            varpi: varpi.rem_euclid(TWO_PI),
            a: elem.a,
        }
    };

    let record = |particles: &[StateVector],
                  active: &[bool],
                  bodies: &[MassiveBody],
                  particle_series: &mut Vec<Vec<Option<AngleSample>>>,
                  p9_samples: &mut Vec<AngleSample>| {
        p9_samples.push(sample(&bodies[0].state));
        for (i, p) in particles.iter().enumerate() {
            particle_series[i].push(active[i].then(|| sample(p)));
        }
    };

    record(
        &particles,
        &active,
        &bodies,
        &mut particle_series,
        &mut p9_samples,
    );

    let n_steps = (config.t_total / config.dt).ceil() as usize;
    let sample_every = (config.angle_sample_interval / config.dt).ceil().max(1.0) as usize;
    for step in 1..=n_steps {
        integrator.step(
            &mut bodies,
            &mut particles,
            &mut active,
            config.dt,
            &sim_config,
        );
        if step % sample_every == 0 || step == n_steps {
            record(
                &particles,
                &active,
                &bodies,
                &mut particle_series,
                &mut p9_samples,
            );
        }
    }

    let final_elements: Vec<Option<OrbitalElements>> = particles
        .iter()
        .zip(&active)
        .map(|(p, &alive)| alive.then(|| cartesian_to_elements(p, GM_SUN)))
        .collect();
    let active_count = final_elements.iter().flatten().count();

    ResonanceRunResult {
        t_final: config.t_total,
        final_elements,
        particle_series,
        p9_samples,
        active_count,
        total_count: n_actual,
    }
}

/// Coarse preselection window (fractional in a) for candidate resonances;
/// membership itself is decided by libration of the angle series.
pub const CANDIDATE_A_TOLERANCE: f64 = 0.05;

/// Classify one particle: among catalog resonances within
/// `CANDIDATE_A_TOLERANCE` of its time-averaged semi-major axis, return the
/// one whose resonant-angle series librates with the smallest amplitude.
pub fn classify_particle(
    series: &[Option<AngleSample>],
    p9_samples: &[AngleSample],
    a_p9: f64,
    catalog: &[Resonance],
) -> Option<(Resonance, f64)> {
    let live: Vec<(usize, &AngleSample)> = series
        .iter()
        .enumerate()
        .filter_map(|(k, s)| s.as_ref().map(|s| (k, s)))
        .collect();
    if live.len() < 10 {
        return None;
    }
    let a_mean = live.iter().map(|(_, s)| s.a).sum::<f64>() / live.len() as f64;

    let mut best: Option<(Resonance, f64)> = None;
    for &res in catalog {
        let a_res = (res.semimajor_axis_typed(a_p9) / au(1.0)).value;
        if ((a_mean - a_res) / a_res).abs() > CANDIDATE_A_TOLERANCE {
            continue;
        }
        let angles: Vec<f64> = live
            .iter()
            .map(|(k, s)| {
                resonant_angle_from_angles(s.lambda, p9_samples[*k].lambda, s.varpi, &res)
            })
            .collect();
        let (librating, amplitude) = is_librating(&angles);
        if librating && best.map_or(true, |(_, amp)| amplitude < amp) {
            best = Some((res, amplitude));
        }
    }
    best
}

/// Libration-based resonance census: count librating particles per
/// resonance.
pub fn resonance_census_libration(
    result: &ResonanceRunResult,
    a_p9: f64,
    catalog: &[Resonance],
) -> Vec<(Resonance, usize)> {
    let mut counts: std::collections::HashMap<Resonance, usize> = std::collections::HashMap::new();
    for series in &result.particle_series {
        if let Some((res, _)) = classify_particle(series, &result.p9_samples, a_p9, catalog) {
            *counts.entry(res).or_insert(0) += 1;
        }
    }
    let mut census: Vec<(Resonance, usize)> = counts.into_iter().collect();
    census.sort_by(|a, b| b.1.cmp(&a.1));
    census
}

/// Result of the eccentricity sweep.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct EccentricitySweepResult {
    /// Eccentricity values tested
    pub eccentricities: Vec<f64>,
    /// For each e₉, the libration-based resonance census
    pub census_results: Vec<Vec<(Resonance, usize)>>,
    /// For each e₉, the surviving fraction
    pub survival_fractions: Vec<f64>,
    /// Planet Nine semi-major axis used
    pub a_p9: f64,
}

/// Run the full eccentricity sweep.
///
/// Sweeps e₉ from 0 to 0.7 in steps of 0.1 (8 simulations).
pub fn eccentricity_sweep(quick: bool) -> EccentricitySweepResult {
    let eccentricities: Vec<f64> = (0..=7).map(|i| i as f64 * 0.1).collect();
    let catalog = extended_catalog();

    let mut census_results = Vec::new();
    let mut survival_fractions = Vec::new();

    for (idx, &e) in eccentricities.iter().enumerate() {
        let config = if quick {
            ResonanceSimConfig::quick_test(e)
        } else {
            ResonanceSimConfig::for_eccentricity(e)
        };

        let result = run_planar_simulation(&config, 42 + idx as u64);
        let census = resonance_census_libration(&result, config.a_p9, &catalog);

        survival_fractions.push(result.active_count as f64 / result.total_count.max(1) as f64);
        census_results.push(census);
    }

    EccentricitySweepResult {
        eccentricities,
        census_results,
        survival_fractions,
        a_p9: 600.0,
    }
}

/// Classify surviving particles by their resonance type (libration-based).
pub fn classify_by_resonance_type(result: &ResonanceRunResult, a_p9: f64) -> ResonanceTypeStats {
    let catalog = extended_catalog();
    let mut n_over_1 = 0usize;
    let mut n_over_2 = 0usize;
    let mut high_order = 0usize;
    let mut non_resonant = 0usize;

    for series in &result.particle_series {
        match classify_particle(series, &result.p9_samples, a_p9, &catalog) {
            Some((res, _)) => {
                if res.is_n_over_1() {
                    n_over_1 += 1;
                } else if res.is_n_over_2() {
                    n_over_2 += 1;
                } else {
                    high_order += 1;
                }
            }
            None => non_resonant += 1,
        }
    }

    let total_resonant = n_over_1 + n_over_2 + high_order;
    ResonanceTypeStats {
        n_over_1,
        n_over_2,
        high_order,
        non_resonant,
        total: result.particle_series.len(),
        p_simple: if total_resonant > 0 {
            (n_over_1 + n_over_2) as f64 / total_resonant as f64
        } else {
            0.0
        },
    }
}

/// Statistics on resonance type distribution.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ResonanceTypeStats {
    pub n_over_1: usize,
    pub n_over_2: usize,
    pub high_order: usize,
    pub non_resonant: usize,
    pub total: usize,
    /// Probability that a resonant particle is in N/1 or N/2.
    pub p_simple: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_angle_sample_accessors_match_f64() {
        use uom::si::angle::radian;
        use uom::si::length::astronomical_unit;
        let s = AngleSample {
            lambda: 1.23,
            varpi: 2.34,
            a: 345.6,
        };
        assert!((s.mean_longitude().get::<radian>() - s.lambda).abs() < 1e-9);
        assert!((s.longitude_of_perihelion().get::<radian>() - s.varpi).abs() < 1e-9);
        assert!((s.semi_major_axis().get::<astronomical_unit>() - s.a).abs() < 1e-9);
    }

    #[test]
    fn test_config_for_eccentricity() {
        let config = ResonanceSimConfig::for_eccentricity(0.3);
        assert!((config.a_p9 - 600.0).abs() < 1e-10);
        assert!((config.e_p9 - 0.3).abs() < 1e-10);
        // M7: perihelion-resolution criterion dt <= P(q_min)/20.
        let p_at_qmin = TWO_PI * (config.q_min.powi(3) / GM_SUN).sqrt();
        assert!(
            config.dt <= p_at_qmin / 20.0,
            "dt = {} too coarse",
            config.dt
        );
    }

    #[test]
    fn test_quick_simulation_runs() {
        let config = ResonanceSimConfig::quick_test(0.3);
        let result = run_planar_simulation(&config, 42);
        assert!(result.active_count <= result.total_count);
        // Angle series recorded for every particle, aligned with P9 samples.
        assert_eq!(result.particle_series.len(), result.total_count);
        for s in &result.particle_series {
            assert_eq!(s.len(), result.p9_samples.len());
        }
        assert!(result.p9_samples.len() >= 10);
    }

    #[test]
    fn test_classify_particle_synthetic_librator() {
        // Build a synthetic 2:1 librator (period ratio 1/2): the angle
        // φ = qλ − pλ9 + (p−q)ϖ oscillates about π with amplitude 0.6.
        let res = Resonance::new(2, 1);
        let a_p9 = 600.0;
        let a_res = (res.semimajor_axis_typed(a_p9) / au(1.0)).value;
        let n = 200;
        let mut series = Vec::with_capacity(n);
        let mut p9_samples = Vec::with_capacity(n);
        for k in 0..n {
            let lambda9 = k as f64 * 0.37;
            let phi = PI + 0.6 * (k as f64 * 0.21).sin();
            let varpi = 1.1;
            // φ = qλ − pλ9 + (p−q)ϖ  ⇒  λ = (φ + pλ9 − (p−q)ϖ)/q
            let lambda = (phi + 2.0 * lambda9 - 1.0 * varpi) / 1.0;
            series.push(Some(AngleSample {
                lambda: lambda.rem_euclid(TWO_PI),
                varpi,
                a: a_res,
            }));
            p9_samples.push(AngleSample {
                lambda: lambda9.rem_euclid(TWO_PI),
                varpi: 0.0,
                a: a_p9,
            });
        }
        let catalog = extended_catalog();
        let (res_found, amplitude) =
            classify_particle(&series, &p9_samples, a_p9, &catalog).expect("should librate");
        assert_eq!(res_found, res);
        assert!((amplitude - 0.6).abs() < 0.1, "amplitude = {amplitude}");
    }

    #[test]
    fn test_classify_particle_circulator_rejected() {
        // Same semi-major axis as the 2:1, but λ drifts off-resonance: the
        // a-proximity preselect alone would have called this resonant.
        let res = Resonance::new(2, 1);
        let a_p9 = 600.0;
        let a_res = (res.semimajor_axis_typed(a_p9) / au(1.0)).value;
        let n = 200;
        let mut series = Vec::with_capacity(n);
        let mut p9_samples = Vec::with_capacity(n);
        for k in 0..n {
            let lambda9 = k as f64 * 0.37;
            let lambda = 2.13 * lambda9; // 6.5% off the resonant rate
            series.push(Some(AngleSample {
                lambda: lambda.rem_euclid(TWO_PI),
                varpi: 1.1,
                a: a_res,
            }));
            p9_samples.push(AngleSample {
                lambda: lambda9.rem_euclid(TWO_PI),
                varpi: 0.0,
                a: a_p9,
            });
        }
        let catalog = vec![res];
        assert!(classify_particle(&series, &p9_samples, a_p9, &catalog).is_none());
    }

    /// End-to-end (reduced scale, seeded): simulation → libration census →
    /// resonance-type classification → p_simple → P(all 6 simple). The
    /// headline statistic is now COMPUTED from the pipeline, not quoted —
    /// this is the wiring whose absence let the N/1-classifier inversion
    /// survive (issue #208).
    #[test]
    fn reduced_scale_pipeline_computes_p_simple() {
        use crate::probability_analysis::p_all_simple;
        let config = ResonanceSimConfig {
            n_particles: 40,
            t_total: 3e5 * YEAR_DAYS,
            angle_sample_interval: 2e3 * YEAR_DAYS,
            ..ResonanceSimConfig::for_eccentricity(0.5)
        };
        let result = run_planar_simulation(&config, 2018);
        let stats = classify_by_resonance_type(&result, config.a_p9);
        // Bookkeeping closes.
        assert_eq!(
            stats.n_over_1 + stats.n_over_2 + stats.high_order + stats.non_resonant,
            stats.total
        );
        assert!(stats.total > 0);
        // p_simple is a genuine census-derived probability, and the derived
        // headline P(all 6 in N/1 or N/2) follows from it.
        assert!(
            (0.0..=1.0).contains(&stats.p_simple),
            "p = {}",
            stats.p_simple
        );
        let p6 = p_all_simple(stats.p_simple, 6);
        assert!((0.0..=1.0).contains(&p6));
        assert!(
            (p6 - stats.p_simple.powi(6)).abs() < 1e-12,
            "P(all 6) must be p_simple^6"
        );
    }

    /// Paper-scale eccentricity sweep (8 x 400 particles x 4 Gyr). Run only
    /// deliberately. Asserts the crate's headline: with the full-length
    /// census, the probability that all 6 observed KBOs sit in N/1 or N/2
    /// resonances is below 5% (lib.rs's key finding, previously never
    /// computed end-to-end).
    #[test]
    #[ignore]
    fn test_paper_scale_sweep() {
        use crate::probability_analysis::p_all_simple;
        let sweep = eccentricity_sweep(false);
        assert_eq!(sweep.eccentricities.len(), 8);
        // Resonant occupation should exist for an eccentric P9.
        let last = sweep.census_results.last().unwrap();
        assert!(!last.is_empty());
        // Headline: recompute the classification at the nominal e9 = 0.5 and
        // pin P(all 6 simple) < 0.05.
        let config = ResonanceSimConfig::for_eccentricity(0.5);
        let result = run_planar_simulation(&config, 42);
        let stats = classify_by_resonance_type(&result, config.a_p9);
        let p6 = p_all_simple(stats.p_simple, 6);
        assert!(
            p6 < 0.05,
            "P(all 6 in N/1 or N/2) = {p6:.4} (p_simple = {:.3})",
            stats.p_simple
        );
    }
}
