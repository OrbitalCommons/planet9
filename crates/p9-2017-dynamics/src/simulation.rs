//! N-body simulation for the paper's central experiment (Batygin &
//! Morbidelli 2017): scattered-disk test particles under the giant planets
//! and Planet Nine.
//!
//! Physics configuration:
//! - Jupiter+Saturn+Uranus absorbed into the orbit-averaged J2 field with
//!   monopole boost (`p9_core::forces::ExtraForce::J2Jsu`),
//! - Neptune and Planet Nine integrated directly,
//! - dt = 3000 d ≈ P_Neptune/20 (the WHM resolution criterion for the
//!   innermost direct body),
//! - removal at r < 20 AU or r > 10,000 AU as configured.
//!
//! Fiducial paper scale: 6,000 particles, a ∈ (150, 750) AU, q ∈ (30, 36)
//! AU, 4 Gyr (behind `#[ignore]`/`default_paper`); the reduced-scale default
//! (`reduced_scale`) integrates 100 particles for 50 Myr.

use p9_core::constants::{DEG2RAD, GM_SUN, GYR_DAYS, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{cartesian_to_elements, OrbitalElements, P9Params, SimConfig as CoreConfig};
use rand::Rng;
use rand_distr::{Distribution, Uniform};
use serde::{Deserialize, Serialize};

/// Simulation configuration matching the paper.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    /// Planet Nine parameters
    pub p9: P9Params,
    /// Number of test particles
    pub n_particles: usize,
    /// Integration time in Gyr
    pub t_gyr: f64,
    /// Timestep (days); 3000 d ≈ P_Neptune/20
    pub dt: f64,
    /// Interval between element snapshots (days)
    pub snapshot_interval: f64,
    /// Minimum initial semi-major axis (AU)
    pub a_min: f64,
    /// Maximum initial semi-major axis (AU)
    pub a_max: f64,
    /// Minimum initial perihelion distance (AU)
    pub q_min: f64,
    /// Maximum initial perihelion distance (AU)
    pub q_max: f64,
    /// Inner removal boundary (AU)
    pub r_remove_inner: f64,
    /// Outer removal boundary (AU)
    pub r_remove_outer: f64,
}

impl SimConfig {
    /// Fiducial paper-scale configuration: a_9=700, e_9=0.6, m_9=10 ME,
    /// 6000 particles, 4 Gyr. Run behind `#[ignore]`.
    pub fn default_paper() -> Self {
        Self {
            p9: P9Params {
                mass_earth: 10.0,
                a: 700.0,
                e: 0.6,
                i: 20.0 * DEG2RAD,
                omega: 150.0 * DEG2RAD,
                omega_big: 100.0 * DEG2RAD,
                mean_anomaly: 0.0,
            },
            n_particles: 6_000,
            t_gyr: 4.0,
            dt: 3000.0,
            snapshot_interval: 10e6 * YEAR_DAYS,
            a_min: 150.0,
            a_max: 750.0,
            q_min: 30.0,
            q_max: 36.0,
            r_remove_inner: 20.0,
            r_remove_outer: 10_000.0,
        }
    }

    /// Reduced-scale default: 100 particles, 50 Myr — runs in seconds in
    /// release mode while exercising the full physics stack.
    pub fn reduced_scale() -> Self {
        Self {
            n_particles: 100,
            t_gyr: 0.05,
            snapshot_interval: 1e6 * YEAR_DAYS,
            ..Self::default_paper()
        }
    }

    /// Mass variation: m_9 = 5 Earth masses.
    pub fn mass_5() -> Self {
        let mut config = Self::default_paper();
        config.p9.mass_earth = 5.0;
        config
    }

    /// Mass variation: m_9 = 20 Earth masses.
    pub fn mass_20() -> Self {
        let mut config = Self::default_paper();
        config.p9.mass_earth = 20.0;
        config
    }
}

/// Observability selection criteria from the paper.
pub struct ObservabilityCriteria {
    pub q_max: f64,
    pub i_max_deg: f64,
}

impl ObservabilityCriteria {
    pub fn default_paper() -> Self {
        Self {
            q_max: 100.0,
            i_max_deg: 40.0,
        }
    }

    pub fn is_observable(&self, q: f64, i_deg: f64) -> bool {
        q <= self.q_max && i_deg <= self.i_max_deg
    }
}

/// Generate initial orbital elements for test particles.
///
/// Semi-major axes uniform in (a_min, a_max), perihelion distances
/// uniform in (q_min, q_max), inclinations drawn from Rayleigh distribution
/// with sigma = 5 deg.
pub fn generate_initial_conditions<R: Rng>(config: &SimConfig, rng: &mut R) -> Vec<TestParticle> {
    let a_dist = Uniform::new(config.a_min, config.a_max);
    let q_dist = Uniform::new(config.q_min, config.q_max);
    let angle_dist = Uniform::new(0.0_f64, std::f64::consts::TAU);

    let mut particles = Vec::with_capacity(config.n_particles);

    while particles.len() < config.n_particles {
        let a = a_dist.sample(rng);
        let q = q_dist.sample(rng);
        let e = 1.0 - q / a;

        if !(0.0..1.0).contains(&e) {
            continue;
        }

        // Rayleigh distribution for inclination with sigma = 5 deg
        let u: f64 = rng.gen();
        let i = 5.0 * DEG2RAD * (-2.0 * (1.0 - u).ln()).sqrt();

        let omega = angle_dist.sample(rng);
        let omega_big = angle_dist.sample(rng);
        let mean_anomaly = angle_dist.sample(rng);

        particles.push(TestParticle {
            a,
            e,
            i,
            omega,
            omega_big,
            mean_anomaly,
        });
    }

    particles
}

/// Test particle orbital elements.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TestParticle {
    pub a: f64,
    pub e: f64,
    pub i: f64,
    pub omega: f64,
    pub omega_big: f64,
    pub mean_anomaly: f64,
}

impl TestParticle {
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    pub fn inclination_deg(&self) -> f64 {
        self.i / DEG2RAD
    }

    fn elements(&self) -> OrbitalElements {
        OrbitalElements {
            a: self.a,
            e: self.e,
            i: self.i,
            omega: self.omega,
            omega_big: self.omega_big,
            mean_anomaly: self.mean_anomaly,
        }
    }
}

/// Result of an integration run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RunResult {
    /// Final elements for each particle (None when removed).
    pub final_elements: Vec<Option<OrbitalElements>>,
    /// Per-particle maximum inclination (deg) reached over the run.
    pub max_inclination_deg: Vec<f64>,
    /// Fraction of the initial population surviving to the end.
    pub survival_fraction: f64,
    /// Fraction of all particles that reached i > 40 deg at some snapshot
    /// (the paper reports 38% with high-inclination excursions at 4 Gyr).
    pub high_i_fraction: f64,
    /// Survivors passing the paper's observability cuts.
    pub n_observable: usize,
}

/// Run the central experiment: WHM with Neptune + (optionally) Planet Nine
/// direct, J2-averaged Jupiter+Saturn+Uranus, configured removal radii, and
/// inclination tracking at every snapshot interval.
pub fn run_simulation(config: &SimConfig, seed: u64, include_p9: bool) -> RunResult {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let initial = generate_initial_conditions(config, &mut rng);
    let mut particles: Vec<_> = initial
        .iter()
        .map(|p| p.elements().to_state_vector(GM_SUN))
        .collect();
    let n = particles.len();
    let mut active = vec![true; n];

    let mut bodies = vec![planets::neptune_j2000()];
    if include_p9 {
        bodies.push(config.p9.to_body());
    }

    let core_config = CoreConfig {
        dt: config.dt,
        t_start: 0.0,
        t_end: config.t_gyr * GYR_DAYS,
        removal_inner_au: config.r_remove_inner,
        removal_outer_au: config.r_remove_outer,
        snapshot_interval_days: config.snapshot_interval,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    let integrator = WhmIntegrator::with_extra_forces(vec![ExtraForce::J2Jsu]);

    let mut max_i_deg: Vec<f64> = initial.iter().map(|p| p.inclination_deg()).collect();

    let t_total = config.t_gyr * GYR_DAYS;
    let n_steps = (t_total / config.dt).ceil() as usize;
    let snap_every = (config.snapshot_interval / config.dt).ceil().max(1.0) as usize;

    for step in 1..=n_steps {
        integrator.step(
            &mut bodies,
            &mut particles,
            &mut active,
            config.dt,
            &core_config,
        );

        if step % snap_every == 0 || step == n_steps {
            for (idx, state) in particles.iter().enumerate() {
                if !active[idx] {
                    continue;
                }
                let elem = cartesian_to_elements(state, GM_SUN);
                let i_deg = elem.i / DEG2RAD;
                if i_deg > max_i_deg[idx] {
                    max_i_deg[idx] = i_deg;
                }
            }
        }
    }

    let final_elements: Vec<Option<OrbitalElements>> = particles
        .iter()
        .zip(&active)
        .map(|(state, &alive)| alive.then(|| cartesian_to_elements(state, GM_SUN)))
        .collect();

    let n_survivors = active.iter().filter(|&&a| a).count();
    let high_i = max_i_deg.iter().filter(|&&i| i > 40.0).count();

    let criteria = ObservabilityCriteria::default_paper();
    let n_observable = final_elements
        .iter()
        .flatten()
        .filter(|e| criteria.is_observable(e.a * (1.0 - e.e), e.i / DEG2RAD))
        .count();

    RunResult {
        final_elements,
        max_inclination_deg: max_i_deg,
        survival_fraction: n_survivors as f64 / n as f64,
        high_i_fraction: high_i as f64 / n as f64,
        n_observable,
    }
}

/// Initial conditions only (kept for the analysis/plotting paths that need
/// the t = 0 population).
pub fn quick_test_simulation(config: &SimConfig) -> Vec<TestParticle> {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(42);
    generate_initial_conditions(config, &mut rng)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_config() {
        let config = SimConfig::default_paper();
        assert!((config.p9.mass_earth - 10.0).abs() < 0.01);
        assert!((config.p9.a - 700.0).abs() < 0.01);
        assert_eq!(config.n_particles, 6_000);
        // dt resolves Neptune: P_N / 20 ≈ 3009 d.
        assert!(config.dt <= p9_core::constants::NEPTUNE_PERIOD_DAYS / 20.0);
    }

    #[test]
    fn test_generate_particles() {
        use rand::SeedableRng;
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let config = SimConfig::default_paper();
        let particles = generate_initial_conditions(&config, &mut rng);
        assert_eq!(particles.len(), config.n_particles);

        for p in &particles {
            assert!(p.a >= config.a_min && p.a <= config.a_max);
            assert!(p.e >= 0.0 && p.e < 1.0);
            let q = p.perihelion();
            assert!(q >= config.q_min - 1.0 && q <= config.q_max + 1.0);
        }
    }

    #[test]
    fn test_observability_criteria() {
        let criteria = ObservabilityCriteria::default_paper();
        assert!(criteria.is_observable(50.0, 20.0));
        assert!(!criteria.is_observable(150.0, 20.0));
        assert!(!criteria.is_observable(50.0, 50.0));
    }

    #[test]
    fn test_mass_variants() {
        let c5 = SimConfig::mass_5();
        let c20 = SimConfig::mass_20();
        assert!((c5.p9.mass_earth - 5.0).abs() < 0.01);
        assert!((c20.p9.mass_earth - 20.0).abs() < 0.01);
    }

    #[test]
    fn test_initial_high_inclination_fraction() {
        let particles = quick_test_simulation(&SimConfig::default_paper());
        let high_i = particles
            .iter()
            .filter(|p| p.inclination_deg() > 40.0)
            .count();
        let fraction = high_i as f64 / particles.len() as f64;
        // Initial conditions (Rayleigh sigma = 5 deg) have essentially no
        // i > 40 deg members; excursions must come from the dynamics.
        assert!(fraction < 0.01, "initial high-i fraction = {}", fraction);
    }

    #[test]
    fn test_smoke_integration_runs() {
        // Tiny end-to-end run: integrator wired, removal applied, snapshots
        // tracked. 100 kyr -> ~12,000 steps x 20 particles.
        let config = SimConfig {
            n_particles: 20,
            t_gyr: 1e-4,
            snapshot_interval: 1e4 * YEAR_DAYS,
            ..SimConfig::default_paper()
        };
        let result = run_simulation(&config, 42, true);
        assert_eq!(result.final_elements.len(), 20);
        assert!(result.survival_fraction > 0.0);
        // Max inclination tracking initialized from the initial conditions.
        assert!(result.max_inclination_deg.iter().all(|&i| i >= 0.0));
    }

    /// Qualitative regression for the paper's central mechanism: with
    /// Planet Nine the population develops inclination excursions that the
    /// control run (giants only) does not. Reduced scale; run with
    /// `cargo test --release -- --ignored`.
    #[test]
    #[ignore]
    fn test_p9_drives_inclination_excursions() {
        let config = SimConfig::reduced_scale();
        let with_p9 = run_simulation(&config, 42, true);
        let control = run_simulation(&config, 42, false);

        let max_gain = |r: &RunResult| {
            r.max_inclination_deg
                .iter()
                .cloned()
                .fold(0.0_f64, f64::max)
        };
        let gain_p9 = max_gain(&with_p9);
        let gain_control = max_gain(&control);
        assert!(
            gain_p9 > gain_control + 1.0,
            "P9 run max i = {gain_p9:.1} deg should exceed control = {gain_control:.1} deg"
        );
    }

    /// Free physical-invariant test (the workspace previously never used
    /// `p9_core::integrator::whm::jacobi_constant`): for a *circular*,
    /// coplanar Planet Nine the rotating-frame Jacobi constant of a test
    /// particle is conserved along the WHM trajectory.
    #[test]
    fn test_jacobi_constant_conserved_circular_p9() {
        use nalgebra::Vector3;
        use p9_core::constants::TWO_PI;
        use p9_core::integrator::whm::jacobi_constant;
        use p9_core::types::StateVector;

        let a9 = 700.0;
        let m9 = 10.0 * p9_core::constants::EARTH_MASS_SOLAR;
        let mu = m9 / (1.0 + m9);
        let gm_total = GM_SUN * (1.0 + m9);
        let n9 = (gm_total / (a9 * a9 * a9)).sqrt(); // rad/day

        // Circular coplanar P9.
        let p9 = P9Params {
            mass_earth: 10.0,
            a: a9,
            e: 0.0,
            i: 0.0,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        let mut bodies = vec![p9.to_body()];

        // Coplanar eccentric test particle well interior to P9.
        let elem = OrbitalElements {
            a: 300.0,
            e: 0.5,
            i: 0.0,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 1.0,
        };
        let mut particles = vec![elem.to_state_vector(GM_SUN)];
        let mut active = vec![true];

        let core_config = CoreConfig {
            dt: 2000.0,
            t_start: 0.0,
            t_end: 1e9,
            removal_inner_au: 1.0,
            removal_outer_au: 1e5,
            snapshot_interval_days: 1e9,
            hybrid_changeover_hill: 3.0,
            bs_epsilon: 1e-11,
        };
        let integrator = WhmIntegrator::new();

        // Rotating-frame, barycentric, nondimensional Jacobi constant.
        let jacobi = |particle: &StateVector, p9_state: &StateVector| -> f64 {
            let theta = p9_state.pos.y.atan2(p9_state.pos.x);
            let (s, c) = theta.sin_cos();
            let rot = |v: &Vector3<f64>| Vector3::new(c * v.x + s * v.y, -s * v.x + c * v.y, v.z);

            // Barycentric position/velocity (heliocentric frame, Sun at 0).
            let r_b = particle.pos - mu * p9_state.pos;
            let v_b = particle.vel - mu * p9_state.vel;

            // Nondimensionalize: lengths by a9, velocities by a9 * n9.
            let rho = rot(&r_b) / a9;
            let v_rot_frame = rot(&v_b) / (a9 * n9);
            // Rotating-frame velocity: u = v - ω ẑ × ρ with ω = 1.
            let u = v_rot_frame - Vector3::new(-rho.y, rho.x, 0.0);

            jacobi_constant(&StateVector::new(rho, u), mu, 1.0)
        };

        let c0 = jacobi(&particles[0], &bodies[0].state);
        let mut max_drift = 0.0_f64;
        // ~10 P9 orbits.
        let n_steps = (10.0 * TWO_PI / (n9 * core_config.dt)).ceil() as usize;
        for _ in 0..n_steps {
            integrator.step(
                &mut bodies,
                &mut particles,
                &mut active,
                core_config.dt,
                &core_config,
            );
            let c = jacobi(&particles[0], &bodies[0].state);
            max_drift = max_drift.max(((c - c0) / c0).abs());
        }
        assert!(
            max_drift < 1e-5,
            "Jacobi constant drift {max_drift:.2e} too large"
        );
    }

    /// Paper-scale configuration (6000 particles, 4 Gyr). Asserts the
    /// headline 38% high-inclination-excursion fraction within a loose
    /// band. Run only deliberately: this is hours of CPU.
    #[test]
    #[ignore]
    fn test_paper_scale_high_i_fraction() {
        let result = run_simulation(&SimConfig::default_paper(), 42, true);
        assert!(
            (0.1..0.7).contains(&result.high_i_fraction),
            "high-i fraction = {:.2} (paper: 0.38)",
            result.high_i_fraction
        );
    }
}
