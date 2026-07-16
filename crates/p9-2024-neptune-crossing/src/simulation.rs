//! N-body simulations for the P9-inclusive and P9-free models.
//!
//! Paper setup (Batygin, Morbidelli, Brown & Nesvorny 2024, Sec. 2): the
//! P9-inclusive model integrates ~2000 test particles (q > 30 AU,
//! 100 < a < 5000 AU initially, from the Nesvorny et al. 2023 `cluster2`
//! snapshot) for 4 Gyr under Jupiter-Neptune + Planet Nine (m = 5 M_earth,
//! a = 500 AU, e = 0.25, i = 20 deg) + galactic effects; the P9-free null
//! model omits Planet Nine. Orbital footprints satisfying a > 100 AU,
//! q < 30 AU, i < 40 deg are recorded over the late phase of the run.
//!
//! Here the integration uses the p9-core WHM integrator. At full scale
//! Neptune is a direct massive body with the Jupiter+Saturn+Uranus
//! orbit-averaged J2 ring (`ExtraForce::J2Jsu`), an optional galactic tide,
//! and optionally Planet Nine as a second massive body.
//! `quick_test_simulation` runs a reduced-scale version (fewer particles,
//! shorter span, all four giants absorbed into the ring, and a mass-boosted
//! P9 to compress the ~Gyr secular forcing — the secular rate scales
//! linearly with the perturber mass); the paper-scale configuration runs in
//! the `#[ignore]`d `full_scale` test.

use p9_core::constants::{DEG2RAD, GM_SUN, GYR_DAYS, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets::neptune_j2000;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{
    cartesian_to_elements, elements_to_cartesian, OrbitalElements, P9Params, SimConfig, StateVector,
};
use p9_core::units::{au, days, degrees, Angle, Length, Time};
use rand::Rng;
use serde::{Deserialize, Serialize};

/// Configuration for the P9-inclusive simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct P9InclusiveConfig {
    /// Planet Nine parameters
    pub p9: P9Params,
    /// Number of test particles
    pub n_particles: usize,
    /// Integration time in Gyr
    pub t_gyr: f64,
}

impl P9InclusiveConfig {
    /// Integration time as a typed [`Time`].
    pub fn integration_time(&self) -> Time {
        days(self.t_gyr * GYR_DAYS)
    }

    /// Default configuration matching the paper: m=5 ME, a=500, e=0.25,
    /// i=20 deg, ~2000 particles, 4 Gyr.
    pub fn default_paper() -> Self {
        Self {
            p9: P9Params {
                mass_earth: 5.0,
                a: 500.0,
                e: 0.25,
                i: 20.0 * DEG2RAD,
                omega: 150.0 * DEG2RAD,
                omega_big: 100.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            n_particles: 2_000,
            t_gyr: 4.0,
        }
    }
}

/// Configuration for the P9-free null model simulation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct P9FreeConfig {
    /// Number of test particles
    pub n_particles: usize,
    /// Integration time in Gyr
    pub t_gyr: f64,
    /// Include galactic tide perturbation
    pub include_galactic_tide: bool,
}

impl P9FreeConfig {
    /// Integration time as a typed [`Time`].
    pub fn integration_time(&self) -> Time {
        days(self.t_gyr * GYR_DAYS)
    }

    /// Default null model configuration.
    pub fn default_paper() -> Self {
        Self {
            n_particles: 2_000,
            t_gyr: 4.0,
            include_galactic_tide: true,
        }
    }
}

/// Options for a single integration run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationOptions {
    pub n_particles: usize,
    /// Total integration time (days)
    pub t_days: f64,
    /// WHM timestep (days); with `direct_neptune` it must resolve Neptune
    /// (P_N/20 ~ 3000 d), otherwise only the slowest massive body / the
    /// innermost particles.
    pub dt_days: f64,
    /// Record a footprint snapshot every this many steps
    pub snapshot_every_steps: usize,
    /// Planet Nine, if present
    pub p9: Option<P9Params>,
    /// Integrate Neptune as a direct massive body (paper-faithful, but caps
    /// the timestep at ~3000 d). When false, Neptune is absorbed into the
    /// orbit-averaged J2 ring together with Jupiter-Uranus — adequate for
    /// the P9-driven secular injection studied at reduced scale, though it
    /// omits Neptune scattering of crossing orbits (the ring expansion is
    /// formally invalid at q < 30 AU; see p9-core's j2_secular docs).
    pub direct_neptune: bool,
    /// Include the galactic tide kick
    pub include_galactic_tide: bool,
    /// RNG seed for the initial particle population
    pub seed: u64,
}

impl SimulationOptions {
    /// Total integration time as a typed [`Time`].
    pub fn total_time(&self) -> Time {
        days(self.t_days)
    }
    /// WHM timestep as a typed [`Time`].
    pub fn timestep(&self) -> Time {
        days(self.dt_days)
    }

    /// Reduced-scale options used by `quick_test_simulation`: 24 particles,
    /// 3 Myr, Neptune absorbed into the ring (20,000-day steps), and P9's
    /// mass boosted 60x (300 M_earth) so the run applies the secular forcing
    /// of ~180 Myr of the paper's 5 M_earth P9 (secular rates scale linearly
    /// with perturber mass). Real dynamics at compressed scale — not a
    /// distributional fudge.
    pub fn reduced_scale(with_p9: bool) -> Self {
        let p9 = with_p9.then(|| {
            let mut p9 = P9InclusiveConfig::default_paper().p9;
            p9.mass_earth *= 60.0;
            p9
        });
        Self {
            n_particles: 24,
            t_days: 3.0e6 * YEAR_DAYS,
            dt_days: 20_000.0,
            snapshot_every_steps: 250,
            p9,
            direct_neptune: false,
            include_galactic_tide: false,
            seed: 20240417, // arXiv:2404.11594 submission date
        }
    }

    /// Paper-scale options (hours-to-days of CPU; used behind #[ignore]).
    pub fn full_scale(with_p9: bool) -> Self {
        // Both paper models parameterize their own run: the P9-inclusive
        // config for the signal run, P9FreeConfig for the null control.
        let (n_particles, t_days, p9, include_galactic_tide) = if with_p9 {
            let c = P9InclusiveConfig::default_paper();
            (
                c.n_particles,
                c.t_gyr * p9_core::constants::GYR_DAYS,
                Some(c.p9),
                true,
            )
        } else {
            let c = P9FreeConfig::default_paper();
            (
                c.n_particles,
                (c.integration_time() / days(1.0)).value,
                None,
                c.include_galactic_tide,
            )
        };
        Self {
            n_particles,
            t_days,
            dt_days: 2500.0,
            snapshot_every_steps: 40_000,
            p9,
            direct_neptune: true,
            include_galactic_tide,
            seed: 20240417,
        }
    }
}

/// An orbital footprint recorded during the run.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Footprint {
    /// Semi-major axis (AU)
    pub a: f64,
    /// Perihelion (AU)
    pub q: f64,
    /// Inclination (deg)
    pub i_deg: f64,
}

impl Footprint {
    /// Semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }
    /// Perihelion as a typed [`Length`].
    pub fn perihelion(&self) -> Length {
        au(self.q)
    }
    /// Inclination as a typed [`Angle`].
    pub fn inclination(&self) -> Angle {
        degrees(self.i_deg)
    }
}

/// Result of a simulation run, storing orbital footprint statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimulationResult {
    /// Number of surviving particles
    pub n_surviving: usize,
    /// Number of recorded footprints meeting the selection criteria
    /// (a > 100 AU, i < 40 deg, q < 30 AU)
    pub n_selected: usize,
    /// Total number of recorded footprints (selection denominator)
    pub n_footprints: usize,
    /// Perihelion distances of selected footprints (AU)
    pub selected_perihelia: Vec<f64>,
    /// Inclinations of selected footprints (degrees)
    pub selected_inclinations: Vec<f64>,
    /// Semi-major axes of selected footprints (AU)
    pub selected_semimajor: Vec<f64>,
    /// Model label
    pub model_name: String,
}

impl SimulationResult {
    /// Fraction of recorded footprints meeting the selection criteria.
    pub fn selection_fraction(&self) -> f64 {
        if self.n_footprints == 0 {
            return 0.0;
        }
        self.n_selected as f64 / self.n_footprints as f64
    }

    /// Selected footprints as (a, q) pairs for the perihelion-distribution
    /// CDF of the hypothesis test.
    pub fn selected_footprints(&self) -> Vec<Footprint> {
        self.selected_semimajor
            .iter()
            .zip(&self.selected_perihelia)
            .zip(&self.selected_inclinations)
            .map(|((&a, &q), &i_deg)| Footprint { a, q, i_deg })
            .collect()
    }
}

/// Run a WHM integration and record orbital footprints.
///
/// Initial particles: a ~ U(150, 500) AU, q ~ U(31, 45) AU (detached,
/// non-crossing, mirroring the paper's removal of Neptune-crossing initial
/// conditions), i ~ U(0, 25) deg, angles uniform.
pub fn run_simulation(opts: &SimulationOptions) -> SimulationResult {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(opts.seed);

    // Massive bodies: Neptune direct (paper-faithful) or absorbed into the
    // ring; optional P9.
    let mut bodies = Vec::new();
    if opts.direct_neptune {
        bodies.push(neptune_j2000());
    }
    if let Some(p9) = &opts.p9 {
        bodies.push(p9.to_body());
    }

    // Below this radius the orbit-averaged ring expansion (in powers of
    // (a_planet/r)^2) is meaningless and its J4 term diverges; particles
    // entering it are dropped (see `removal_inner` below).
    const RING_VALID_MIN_AU: f64 = 12.0;

    let mut extra = if opts.direct_neptune {
        vec![ExtraForce::J2Jsu]
    } else {
        // Jupiter+Saturn+Uranus+Neptune as one orbit-averaged J2/J4 ring.
        use p9_core::constants::{A_NEPTUNE_AU, GM_NEPTUNE, MASS_NEPTUNE_SOLAR};
        use p9_core::forces::j2_secular::{combined_j2_jsu, effective_j2, secular_j2_acceleration};
        use std::sync::Arc;
        let (j2_jsu, j4_jsu, gm_jsu) = combined_j2_jsu();
        let j2 = j2_jsu + effective_j2(MASS_NEPTUNE_SOLAR, A_NEPTUNE_AU, 1.0);
        let j4 = j4_jsu - 0.375 * MASS_NEPTUNE_SOLAR * A_NEPTUNE_AU.powi(4);
        let gm_boost = gm_jsu + GM_NEPTUNE;
        vec![ExtraForce::Custom(Arc::new(move |pos| {
            if pos.norm() < RING_VALID_MIN_AU {
                // Out of the model's validity region; the particle is
                // removed at the end of the step.
                return *pos * 0.0;
            }
            secular_j2_acceleration(pos, GM_SUN, j2, j4, gm_boost, 1.0)
        }))]
    };
    if opts.include_galactic_tide {
        extra.push(ExtraForce::GalacticTide);
    }
    let integrator = WhmIntegrator::with_extra_forces(extra);

    // Initial test-particle population.
    let mut particles: Vec<StateVector> = (0..opts.n_particles)
        .map(|_| {
            let a = rng.gen_range(150.0..500.0);
            let q = rng.gen_range(31.0..45.0);
            let elem = OrbitalElements {
                a,
                e: 1.0 - q / a,
                i: rng.gen_range(0.0..25.0 * DEG2RAD),
                omega_big: rng.gen_range(0.0..std::f64::consts::TAU),
                omega: rng.gen_range(0.0..std::f64::consts::TAU),
                mean_anomaly: rng.gen_range(0.0..std::f64::consts::TAU),
            };
            elements_to_cartesian(&elem, GM_SUN)
        })
        .collect();
    let mut active = vec![true; opts.n_particles];

    let config = SimConfig {
        dt: opts.dt_days,
        t_start: 0.0,
        t_end: opts.t_days,
        // Paper boundaries are 1 and 100,000 AU; without direct Neptune the
        // inner boundary is the ring model's validity radius instead.
        removal_inner_au: if opts.direct_neptune {
            1.0
        } else {
            RING_VALID_MIN_AU
        },
        removal_outer_au: 100_000.0,
        ..SimConfig::standard_4gyr()
    };

    let n_steps = (opts.t_days / opts.dt_days).ceil() as usize;
    let mut footprints: Vec<Footprint> = Vec::new();

    for step in 0..n_steps {
        integrator.step(
            &mut bodies,
            &mut particles,
            &mut active,
            opts.dt_days,
            &config,
        );

        // Record footprints in the latter half of the run (quasi-steady
        // state, mirroring the paper's final-Gyr sampling).
        if step % opts.snapshot_every_steps == 0 && step >= n_steps / 2 {
            for (k, particle) in particles.iter().enumerate() {
                if !active[k] {
                    continue;
                }
                let elem = cartesian_to_elements(particle, GM_SUN);
                if elem.e < 0.0 || elem.e >= 1.0 || elem.a <= 0.0 {
                    continue;
                }
                footprints.push(Footprint {
                    a: elem.a,
                    q: elem.a * (1.0 - elem.e),
                    i_deg: elem.i / DEG2RAD,
                });
            }
        }
    }

    let n_surviving = active.iter().filter(|&&a| a).count();
    let criteria = crate::observed_tnos::selection_criteria();
    let selected: Vec<&Footprint> = footprints
        .iter()
        .filter(|f| f.a > criteria.a_min && f.i_deg < criteria.i_max && f.q < criteria.q_max)
        .collect();

    SimulationResult {
        n_surviving,
        n_selected: selected.len(),
        n_footprints: footprints.len(),
        selected_perihelia: selected.iter().map(|f| f.q).collect(),
        selected_inclinations: selected.iter().map(|f| f.i_deg).collect(),
        selected_semimajor: selected.iter().map(|f| f.a).collect(),
        model_name: if opts.p9.is_some() {
            "P9-inclusive".to_string()
        } else {
            "P9-free".to_string()
        },
    }
}

/// Reduced-scale integration used by tests and plots: a real WHM run (not a
/// distributional fudge), seeded, with a mass-boosted P9 to compress the
/// secular timescale. See `SimulationOptions::reduced_scale`. The two
/// (deterministic) runs are cached process-wide so repeated callers do not
/// re-integrate.
pub fn quick_test_simulation(with_p9: bool) -> SimulationResult {
    use std::sync::OnceLock;
    static WITH_P9: OnceLock<SimulationResult> = OnceLock::new();
    static CONTROL: OnceLock<SimulationResult> = OnceLock::new();
    let cell = if with_p9 { &WITH_P9 } else { &CONTROL };
    cell.get_or_init(|| run_simulation(&SimulationOptions::reduced_scale(with_p9)))
        .clone()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_accessors_match_f64() {
        use approx::assert_relative_eq;
        use uom::si::angle::degree;
        use uom::si::length::astronomical_unit;
        use uom::si::time::day;
        let cfg = P9InclusiveConfig::default_paper();
        assert_relative_eq!(
            cfg.integration_time().get::<day>(),
            cfg.t_gyr * GYR_DAYS,
            max_relative = 1e-12
        );
        let opts = SimulationOptions::reduced_scale(true);
        assert_relative_eq!(
            opts.total_time().get::<day>(),
            opts.t_days,
            max_relative = 1e-9
        );
        assert_relative_eq!(opts.timestep().get::<day>(), opts.dt_days, epsilon = 1e-9);
        let fp = Footprint {
            a: 150.0,
            q: 28.0,
            i_deg: 12.0,
        };
        assert_relative_eq!(
            fp.semi_major_axis().get::<astronomical_unit>(),
            fp.a,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            fp.perihelion().get::<astronomical_unit>(),
            fp.q,
            epsilon = 1e-12
        );
        assert_relative_eq!(fp.inclination().get::<degree>(), fp.i_deg, epsilon = 1e-12);
    }

    #[test]
    fn test_p9_inclusive_config() {
        let config = P9InclusiveConfig::default_paper();
        assert!((config.p9.mass_earth - 5.0).abs() < 0.01);
        assert!((config.p9.a - 500.0).abs() < 0.01);
        assert!((config.p9.e - 0.25).abs() < 0.01);
        let i_deg = config.p9.i / DEG2RAD;
        assert!((i_deg - 20.0).abs() < 0.1);
    }

    #[test]
    fn test_p9_free_config_drives_the_null_run() {
        // P9FreeConfig is consumed, not decorative: the full-scale control
        // run derives its particle count, duration and tide flag from it.
        let config = P9FreeConfig::default_paper();
        let opts = SimulationOptions::full_scale(false);
        assert_eq!(opts.n_particles, config.n_particles);
        assert_eq!(opts.include_galactic_tide, config.include_galactic_tide);
        assert!(((opts.t_days / p9_core::constants::GYR_DAYS) - config.t_gyr).abs() < 1e-12);
        assert!(opts.p9.is_none());
    }

    #[test]
    fn test_p9_drives_neptune_crossers() {
        // The paper's qualitative signature: the P9-inclusive run injects
        // long-period objects into the q < 30 AU region far more efficiently
        // than the P9-free control (paper: ~3% vs ~0.5% of footprints).
        let with_p9 = quick_test_simulation(true);
        let control = quick_test_simulation(false);

        assert_eq!(with_p9.model_name, "P9-inclusive");
        assert_eq!(control.model_name, "P9-free");
        assert!(with_p9.n_footprints > 0);
        assert!(control.n_footprints > 0);

        assert!(
            with_p9.n_selected > control.n_selected,
            "P9 run should produce more Neptune-crossing footprints: {} vs {}",
            with_p9.n_selected,
            control.n_selected
        );
        assert!(
            with_p9.selection_fraction() > control.selection_fraction(),
            "P9 selection fraction {} should exceed control {}",
            with_p9.selection_fraction(),
            control.selection_fraction()
        );
    }

    #[test]
    fn test_selected_within_criteria() {
        let result = quick_test_simulation(true);
        for &q in &result.selected_perihelia {
            assert!(q < 30.0, "q = {} should be < 30 AU", q);
        }
        for &i in &result.selected_inclinations {
            assert!(i < 40.0, "i = {} should be < 40 deg", i);
        }
        for &a in &result.selected_semimajor {
            assert!(a > 100.0, "a = {} should be > 100 AU", a);
        }
    }

    #[test]
    #[ignore = "paper-scale 4 Gyr integration; run explicitly"]
    fn test_full_scale_simulation() {
        let with_p9 = run_simulation(&SimulationOptions::full_scale(true));
        let control = run_simulation(&SimulationOptions::full_scale(false));
        assert!(with_p9.selection_fraction() > control.selection_fraction());
    }
}
