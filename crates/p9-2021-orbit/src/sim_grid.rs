//! Model 1 of Brown & Batygin (2021), Section 2: the N-body simulation grid.
//!
//! The paper integrates scattered-disk test particles (a in 150-500 AU,
//! q in 30-50 AU, i in 0-25 deg) for 4 Gyr under the four giant planets plus
//! a Planet Nine drawn from a manually-constructed 121-point grid over
//! (m9, a9, e9, i9), with m9 in {3, 4, 5, 6, 7, 10, 12, 20} Earth masses.
//! Particles are removed at r < 4.5 AU or r > 10,000 AU; element snapshots
//! after the first Gyr are pooled as samples of the quasi-steady-state
//! distribution.
//!
//! Documented deviations from the paper (scale reductions, repo convention):
//!
//! - **Grid geometry.** The paper's exact 121-point manual grid is not
//!   published. We use a *regular* grid spanning the same extent
//!   (a9 in 300-800 AU, e9 in 0.1-0.5, i9 in 5-25 deg, the paper's 8 masses),
//!   sized by the [`GridScale`] knob: `Paper` builds a 121-point surrogate
//!   (8 masses x 15 (a9, e9, i9) combinations + 1 central point), `Reduced`
//!   and `Smoke` are coarse sub-grids that run in test time.
//! - **Perturber treatment.** The paper integrates all four giant planets
//!   directly with a 300 d timestep. We use the workspace-standard reduced
//!   configuration: Jupiter+Saturn+Uranus absorbed into the orbit-averaged
//!   J2 field with monopole boost (`ExtraForce::J2Jsu`), Neptune and Planet
//!   Nine direct, dt = 3000 d (P_Neptune/20, the WHM resolution criterion
//!   for the innermost direct body). This is valid for q >= 30 AU, which the
//!   initial conditions guarantee.
//! - **Inner removal boundary.** The paper removes particles at r < 4.5 AU
//!   after close encounters with the directly-integrated gas giants. With
//!   JSU absorbed into the averaged J2 field there are no close encounters
//!   to scatter particles that deep; the r < 4.5 AU boundary is kept (it is
//!   the configured `removal_inner_au`) but in practice particles leave
//!   through perihelion detachment or the r > 10,000 AU boundary. Particles
//!   that random-walk inside q ~ 30 AU are outside the J2 field's validity
//!   for a few orbits before removal — a documented cost of the averaged
//!   treatment.
//! - **Duration / population.** `Reduced` integrates 12 grid points x 200
//!   particles x 10 Myr (machinery validation; the secular clustering that
//!   the inference keys on is only partially developed). `Paper` is the
//!   121 x 16,800 x 4 Gyr configuration, run only behind `#[ignore]`
//!   (estimated ~3,000 CPU-days, see REPRODUCTION_NOTES.md).

use rand::{Rng, SeedableRng};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use p9_core::constants::{DEG2RAD, GM_SUN, GYR_DAYS, TWO_PI, YEAR_DAYS};
use p9_core::forces::ExtraForce;
use p9_core::initial_conditions::planets;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{
    cartesian_to_elements, OrbitalElements, P9Params, SimConfig as CoreSimConfig, StateVector,
};
use p9_core::units::{au, days, degrees, earth_masses, radians, Angle, Length, Mass, Time};

/// The paper's Planet Nine mass grid (Earth masses), Brown & Batygin (2021)
/// Section 2.
pub const P9_MASSES_PAPER: [f64; 8] = [3.0, 4.0, 5.0, 6.0, 7.0, 10.0, 12.0, 20.0];

/// Extent of the (m9, a9, e9, i9) grid. Also used as the uniform-prior box
/// for the MCMC stage (Model 4).
pub const M9_RANGE: (f64, f64) = (3.0, 20.0);
/// Semi-major axis extent of the grid (AU).
pub const A9_RANGE: (f64, f64) = (300.0, 800.0);
/// Eccentricity extent of the grid.
pub const E9_RANGE: (f64, f64) = (0.1, 0.5);
/// Inclination extent of the grid (degrees).
pub const I9_RANGE_DEG: (f64, f64) = (5.0, 25.0);

/// Grid size / fidelity knob.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum GridScale {
    /// 4 points x 30 particles x 0.1 Myr — runs in default (debug) test time.
    /// Exercises the full machinery; the dynamics has not had time to do
    /// anything, so the resulting likelihood surface is noise.
    Smoke,
    /// 12 points x 200 particles x 10 Myr — measured ~67 s wall / 12.9
    /// CPU-min (release, 64 cores), behind `#[ignore]`. Secular clustering
    /// is only partially developed.
    Reduced,
    /// 121-point surrogate of the paper grid x 16,800 particles x 4 Gyr.
    /// Estimated ~3,000 CPU-days at the measured reduced-scale throughput
    /// (~265 ns per particle-step); `#[ignore]`d entry point only.
    Paper,
}

/// One Planet Nine parameter point of the simulation grid.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct GridPoint {
    /// Planet Nine mass (Earth masses)
    pub m9: f64,
    /// Planet Nine semi-major axis (AU)
    pub a9: f64,
    /// Planet Nine eccentricity
    pub e9: f64,
    /// Planet Nine inclination (degrees)
    pub i9_deg: f64,
}

impl GridPoint {
    /// Normalized coordinates in the grid box (each dimension mapped to
    /// [0, 1] over the grid extent). This is the GP input space (Model 3).
    pub fn normalized(&self) -> [f64; 4] {
        [
            (self.m9 - M9_RANGE.0) / (M9_RANGE.1 - M9_RANGE.0),
            (self.a9 - A9_RANGE.0) / (A9_RANGE.1 - A9_RANGE.0),
            (self.e9 - E9_RANGE.0) / (E9_RANGE.1 - E9_RANGE.0),
            (self.i9_deg - I9_RANGE_DEG.0) / (I9_RANGE_DEG.1 - I9_RANGE_DEG.0),
        ]
    }

    /// Planet Nine parameters with the orientation angles fixed at zero.
    ///
    /// The recorded particle angles are *relative* to Planet Nine
    /// (delta_varpi, delta_Omega), and the initial disk is uniform in its
    /// own angles, so the absolute orientation of Planet Nine is a gauge
    /// choice; omega = Omega = M = 0 is used for every grid point.
    pub fn p9_params(&self) -> P9Params {
        P9Params {
            mass_earth: self.m9,
            a: self.a9,
            e: self.e9,
            i: self.i9_deg * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        }
    }

    // ---- typed boundary accessors (dimension-checked views of the f64 fields) ----

    /// Planet Nine mass as a [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.m9)
    }
    /// Planet Nine semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a9)
    }
    /// Planet Nine inclination as an [`Angle`].
    pub fn inclination(&self) -> Angle {
        degrees(self.i9_deg)
    }
}

/// Configuration for the simulation grid (Model 1).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimGridConfig {
    pub scale: GridScale,
    /// Base RNG seed; grid point k uses seed `seed + k`.
    pub seed: u64,
    /// Test particles per grid point.
    pub n_particles: usize,
    /// Integration time (Gyr). Paper: 4.0.
    pub t_gyr: f64,
    /// Timestep (days); 3000 d = P_Neptune/20 with Neptune innermost direct.
    pub dt_days: f64,
    /// Fraction of the run discarded before snapshots are pooled.
    /// Paper: 1 Gyr of 4 Gyr = 0.25.
    pub burn_fraction: f64,
    /// Interval between pooled element snapshots (days).
    pub snapshot_interval_days: f64,
    /// Initial disk: semi-major axis range (AU). Paper: 150-500.
    pub disk_a_range: (f64, f64),
    /// Initial disk: perihelion range (AU). Paper: 30-50.
    pub disk_q_range: (f64, f64),
    /// Initial disk: inclination range (degrees), uniform. Paper: 0-25.
    pub disk_i_max_deg: f64,
    /// Removal boundaries (AU). Paper: 4.5 / 10,000.
    pub r_remove_inner: f64,
    pub r_remove_outer: f64,
}

impl SimGridConfig {
    fn base(scale: GridScale, n_particles: usize, t_gyr: f64, snapshot_myr: f64) -> Self {
        Self {
            scale,
            seed: 2021,
            n_particles,
            t_gyr,
            dt_days: 3000.0,
            burn_fraction: 0.25,
            snapshot_interval_days: snapshot_myr * 1e6 * YEAR_DAYS,
            disk_a_range: (150.0, 500.0),
            disk_q_range: (30.0, 50.0),
            disk_i_max_deg: 25.0,
            r_remove_inner: 4.5,
            r_remove_outer: 10_000.0,
        }
    }

    /// Smoke scale: 4 points x 30 particles x 0.1 Myr. Default-test budget.
    pub fn smoke() -> Self {
        Self::base(GridScale::Smoke, 30, 1e-4, 0.02)
    }

    /// Reduced scale (default): 12 points x 200 particles x 10 Myr.
    pub fn reduced() -> Self {
        Self::base(GridScale::Reduced, 200, 0.01, 1.0)
    }

    /// Paper-scale surrogate: 121 points x 16,800 particles x 4 Gyr.
    /// `#[ignore]`d entry point only — estimated CPU-days.
    pub fn paper() -> Self {
        Self::base(GridScale::Paper, 16_800, 4.0, 10.0)
    }

    // ---- typed boundary accessors (dimension-checked views of the f64 fields) ----

    /// Integration time as a [`Time`].
    pub fn integration_time(&self) -> Time {
        days(self.t_gyr * GYR_DAYS)
    }
    /// Integration timestep as a [`Time`].
    pub fn timestep(&self) -> Time {
        days(self.dt_days)
    }
    /// Interval between pooled element snapshots as a [`Time`].
    pub fn snapshot_interval(&self) -> Time {
        days(self.snapshot_interval_days)
    }
    /// Maximum initial-disk inclination as an [`Angle`].
    pub fn max_disk_inclination(&self) -> Angle {
        degrees(self.disk_i_max_deg)
    }
    /// Inner removal boundary as a [`Length`].
    pub fn inner_removal_radius(&self) -> Length {
        au(self.r_remove_inner)
    }
    /// Outer removal boundary as a [`Length`].
    pub fn outer_removal_radius(&self) -> Length {
        au(self.r_remove_outer)
    }

    /// The grid points for this scale (deterministic, no RNG involved).
    ///
    /// - `Smoke`: m9 {5, 10} x a9 {400, 700} at (e9, i9) = (0.3, 15).
    /// - `Reduced`: m9 {5, 10} x a9 {350, 550, 750} x (e9, i9) in
    ///   {(0.2, 10), (0.4, 20)} — e9 and i9 co-vary so all four GP input
    ///   dimensions have spread at 12 points (documented confounding).
    /// - `Paper`: the 8 paper masses x 15 (a9, e9, i9) combinations
    ///   (a9 {300, 425, 550, 675, 800} x e9 {0.1, 0.3, 0.5}, i9 cycling
    ///   {5, 15, 25} over combinations) + 1 central point
    ///   (6 M_E, 550 AU, 0.3, 15 deg) = 121 points, matching the paper's
    ///   count over the same extent. The paper's actual manual grid is not
    ///   published.
    pub fn grid_points(&self) -> Vec<GridPoint> {
        match self.scale {
            GridScale::Smoke => {
                let mut pts = Vec::new();
                for &m9 in &[5.0, 10.0] {
                    for &a9 in &[400.0, 700.0] {
                        pts.push(GridPoint {
                            m9,
                            a9,
                            e9: 0.3,
                            i9_deg: 15.0,
                        });
                    }
                }
                pts
            }
            GridScale::Reduced => {
                let mut pts = Vec::new();
                for &m9 in &[5.0, 10.0] {
                    for &a9 in &[350.0, 550.0, 750.0] {
                        for &(e9, i9_deg) in &[(0.2, 10.0), (0.4, 20.0)] {
                            pts.push(GridPoint { m9, a9, e9, i9_deg });
                        }
                    }
                }
                pts
            }
            GridScale::Paper => {
                let a9_values = [300.0, 425.0, 550.0, 675.0, 800.0];
                let e9_values = [0.1, 0.3, 0.5];
                let i9_values = [5.0, 15.0, 25.0];
                let mut pts = Vec::new();
                for &m9 in P9_MASSES_PAPER.iter() {
                    for (j, &a9) in a9_values.iter().enumerate() {
                        for (k, &e9) in e9_values.iter().enumerate() {
                            // Latin-square assignment (j + k) mod 3: every
                            // (e9, i9) pair occurs (5x per mass), so i9 is NOT
                            // aliased with e9. The previous combo % 3 indexing
                            // locked i9 to e9 in all 121 points, leaving the
                            // GP/MCMC with zero data to separate the two
                            // dimensions.
                            pts.push(GridPoint {
                                m9,
                                a9,
                                e9,
                                i9_deg: i9_values[(j + k) % i9_values.len()],
                            });
                        }
                    }
                }
                // Central point to reach the paper's 121-point count.
                pts.push(GridPoint {
                    m9: 6.0,
                    a9: 550.0,
                    e9: 0.3,
                    i9_deg: 15.0,
                });
                pts
            }
        }
    }
}

/// One pooled element sample of a surviving test particle, relative to
/// Planet Nine. Angles in radians; `delta_varpi`/`delta_omega_big` wrapped
/// to [0, 2 pi).
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ElementSample {
    pub a: f64,
    pub e: f64,
    pub i: f64,
    /// varpi_particle - varpi_P9, wrapped to [0, 2 pi)
    pub delta_varpi: f64,
    /// Omega_particle - Omega_P9, wrapped to [0, 2 pi)
    pub delta_omega_big: f64,
}

impl ElementSample {
    /// Semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }
    /// Inclination as an [`Angle`].
    pub fn inclination(&self) -> Angle {
        radians(self.i)
    }
    /// Apsidal offset Δϖ relative to Planet Nine as an [`Angle`].
    pub fn delta_varpi_angle(&self) -> Angle {
        radians(self.delta_varpi)
    }
    /// Nodal offset ΔΩ relative to Planet Nine as an [`Angle`].
    pub fn delta_omega_big_angle(&self) -> Angle {
        radians(self.delta_omega_big)
    }
}

/// Result of integrating one grid point.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GridPointResult {
    pub point: GridPoint,
    /// RNG seed used for this point's initial conditions (recorded output).
    pub seed: u64,
    pub n_initial: usize,
    /// Particles still active at the end of the run.
    pub n_surviving: usize,
    /// Pooled post-burn snapshot samples of surviving particles.
    pub samples: Vec<ElementSample>,
}

/// Generate the paper's initial scattered disk: a ~ U(150, 500) AU,
/// q ~ U(30, 50) AU (resampled when q > a, which cannot occur for these
/// ranges), i ~ U(0, 25 deg), angles uniform.
pub fn generate_disk<R: Rng>(config: &SimGridConfig, rng: &mut R) -> Vec<StateVector> {
    let mut particles = Vec::with_capacity(config.n_particles);
    while particles.len() < config.n_particles {
        let a = rng.gen_range(config.disk_a_range.0..config.disk_a_range.1);
        let q = rng.gen_range(config.disk_q_range.0..config.disk_q_range.1);
        if q >= a {
            continue;
        }
        let elements = OrbitalElements {
            a,
            e: 1.0 - q / a,
            i: rng.gen_range(0.0..config.disk_i_max_deg) * DEG2RAD,
            omega: rng.gen_range(0.0..TWO_PI),
            omega_big: rng.gen_range(0.0..TWO_PI),
            mean_anomaly: rng.gen_range(0.0..TWO_PI),
        };
        particles.push(elements.to_state_vector(GM_SUN));
    }
    particles
}

fn wrap_two_pi(angle: f64) -> f64 {
    angle.rem_euclid(TWO_PI)
}

/// Integrate one grid point and pool post-burn element snapshots.
///
/// Physics: WHM (democratic heliocentric), Neptune + Planet Nine direct,
/// Jupiter+Saturn+Uranus as the orbit-averaged J2 field (`ExtraForce::J2Jsu`),
/// removal at the configured radii (see module docs for the deviations from
/// the paper's all-four-direct treatment).
pub fn run_grid_point(config: &SimGridConfig, point: &GridPoint, seed: u64) -> GridPointResult {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    let mut particles = generate_disk(config, &mut rng);
    let n = particles.len();
    let mut active = vec![true; n];

    let mut bodies = vec![planets::neptune_j2000(), point.p9_params().to_body()];

    let core_config = CoreSimConfig {
        dt: config.dt_days,
        t_start: 0.0,
        t_end: config.t_gyr * GYR_DAYS,
        removal_inner_au: config.r_remove_inner,
        removal_outer_au: config.r_remove_outer,
        snapshot_interval_days: config.snapshot_interval_days,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    // Sequential per-step kicks: the grid is already parallel over points
    // (`run_grid`), and at a few hundred particles per point the per-step
    // rayon fan-out costs more than the kick itself.
    let mut integrator = WhmIntegrator::with_extra_forces(vec![ExtraForce::J2Jsu]);
    integrator.parallel = false;

    let t_total = config.t_gyr * GYR_DAYS;
    let n_steps = (t_total / config.dt_days).ceil() as usize;
    let snap_every = (config.snapshot_interval_days / config.dt_days)
        .ceil()
        .max(1.0) as usize;
    let burn_steps = (config.burn_fraction * n_steps as f64).floor() as usize;

    let mut samples = Vec::new();

    for step in 1..=n_steps {
        integrator.step(
            &mut bodies,
            &mut particles,
            &mut active,
            config.dt_days,
            &core_config,
        );

        if step > burn_steps && (step % snap_every == 0 || step == n_steps) {
            // P9 elements at this snapshot (its varpi/Omega precess under
            // the J2 field and Neptune).
            let p9_elem = cartesian_to_elements(&bodies[1].state, GM_SUN);
            let p9_varpi = wrap_two_pi(p9_elem.omega + p9_elem.omega_big);
            for (idx, state) in particles.iter().enumerate() {
                if !active[idx] {
                    continue;
                }
                let elem = cartesian_to_elements(state, GM_SUN);
                if elem.a <= 0.0 || elem.e >= 1.0 {
                    continue; // unbound at this instant; not a disk sample
                }
                samples.push(ElementSample {
                    a: elem.a,
                    e: elem.e,
                    i: elem.i,
                    delta_varpi: wrap_two_pi(elem.omega + elem.omega_big - p9_varpi),
                    delta_omega_big: wrap_two_pi(elem.omega_big - p9_elem.omega_big),
                });
            }
        }
    }

    let n_surviving = active.iter().filter(|&&a| a).count();

    GridPointResult {
        point: *point,
        seed,
        n_initial: n,
        n_surviving,
        samples,
    }
}

/// Run the whole grid (parallel over grid points; point k seeded with
/// `config.seed + k`, so results are deterministic and independent of the
/// scheduling order).
pub fn run_grid(config: &SimGridConfig) -> Vec<GridPointResult> {
    let points = config.grid_points();
    points
        .par_iter()
        .enumerate()
        .map(|(k, point)| run_grid_point(config, point, config.seed + k as u64))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::{degree, radian};
    use uom::si::length::astronomical_unit;
    use uom::si::mass::kilogram;
    use uom::si::time::day;

    #[test]
    fn typed_grid_point_accessors_match_f64() {
        let gp = GridPoint {
            m9: 6.2,
            a9: 380.0,
            e9: 0.2,
            i9_deg: 16.0,
        };
        assert_relative_eq!(
            gp.mass().get::<kilogram>(),
            earth_masses(gp.m9).get::<kilogram>(),
            epsilon = 1.0
        );
        assert_relative_eq!(
            gp.semi_major_axis().get::<astronomical_unit>(),
            gp.a9,
            epsilon = 1e-9
        );
        assert_relative_eq!(gp.inclination().get::<degree>(), gp.i9_deg, epsilon = 1e-12);
    }

    #[test]
    fn typed_config_accessors_match_f64() {
        let c = SimGridConfig::reduced();
        assert_relative_eq!(
            c.integration_time().get::<day>(),
            c.t_gyr * GYR_DAYS,
            epsilon = 1e-6
        );
        assert_relative_eq!(c.timestep().get::<day>(), c.dt_days, epsilon = 1e-6);
        assert_relative_eq!(
            c.snapshot_interval().get::<day>(),
            c.snapshot_interval_days,
            epsilon = 1e-6
        );
        assert_relative_eq!(
            c.max_disk_inclination().get::<degree>(),
            c.disk_i_max_deg,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            c.inner_removal_radius().get::<astronomical_unit>(),
            c.r_remove_inner,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            c.outer_removal_radius().get::<astronomical_unit>(),
            c.r_remove_outer,
            epsilon = 1e-9
        );
    }

    #[test]
    fn typed_element_sample_accessors_match_f64() {
        let s = ElementSample {
            a: 250.0,
            e: 0.3,
            i: 0.4,
            delta_varpi: 1.2,
            delta_omega_big: 2.4,
        };
        assert_relative_eq!(
            s.semi_major_axis().get::<astronomical_unit>(),
            s.a,
            epsilon = 1e-9
        );
        assert_relative_eq!(s.inclination().get::<radian>(), s.i, epsilon = 1e-12);
        assert_relative_eq!(
            s.delta_varpi_angle().get::<radian>(),
            s.delta_varpi,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            s.delta_omega_big_angle().get::<radian>(),
            s.delta_omega_big,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_grid_point_counts() {
        assert_eq!(SimGridConfig::smoke().grid_points().len(), 4);
        assert_eq!(SimGridConfig::reduced().grid_points().len(), 12);
        // The paper-scale surrogate matches the paper's 121-point count.
        assert_eq!(SimGridConfig::paper().grid_points().len(), 121);
    }

    #[test]
    fn test_paper_grid_covers_extent_and_masses() {
        let pts = SimGridConfig::paper().grid_points();
        for &m in P9_MASSES_PAPER.iter() {
            assert!(pts.iter().any(|p| p.m9 == m), "mass {m} missing");
        }
        for p in &pts {
            assert!(p.a9 >= A9_RANGE.0 && p.a9 <= A9_RANGE.1);
            assert!(p.e9 >= E9_RANGE.0 && p.e9 <= E9_RANGE.1);
            assert!(p.i9_deg >= I9_RANGE_DEG.0 && p.i9_deg <= I9_RANGE_DEG.1);
            let n = p.normalized();
            assert!(n.iter().all(|&x| (0.0..=1.0).contains(&x)));
        }
    }

    #[test]
    fn test_disk_generation_ranges() {
        let config = SimGridConfig::smoke();
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let disk = generate_disk(&config, &mut rng);
        assert_eq!(disk.len(), config.n_particles);
        for state in &disk {
            let elem = cartesian_to_elements(state, GM_SUN);
            assert!(elem.a > 149.0 && elem.a < 501.0, "a = {}", elem.a);
            let q = elem.a * (1.0 - elem.e);
            assert!(q > 29.0 && q < 51.0, "q = {q}");
            assert!(elem.i >= 0.0 && elem.i <= 25.5 * DEG2RAD);
        }
    }

    #[test]
    fn test_grid_point_run_is_deterministic() {
        let mut config = SimGridConfig::smoke();
        config.n_particles = 10;
        let point = config.grid_points()[0];
        let r1 = run_grid_point(&config, &point, 99);
        let r2 = run_grid_point(&config, &point, 99);
        assert_eq!(r1.samples.len(), r2.samples.len());
        for (s1, s2) in r1.samples.iter().zip(&r2.samples) {
            assert_eq!(s1.a, s2.a);
            assert_eq!(s1.delta_varpi, s2.delta_varpi);
        }
        // Different seed -> different initial conditions -> different samples.
        let r3 = run_grid_point(&config, &point, 100);
        assert!(r1
            .samples
            .iter()
            .zip(&r3.samples)
            .any(|(s1, s3)| s1.a != s3.a));
    }

    #[test]
    fn test_smoke_grid_point_produces_samples() {
        let config = SimGridConfig::smoke();
        let point = config.grid_points()[0];
        let result = run_grid_point(&config, &point, config.seed);
        assert_eq!(result.n_initial, config.n_particles);
        // 0.1 Myr cannot remove the whole disk.
        assert!(result.n_surviving > 0);
        assert!(!result.samples.is_empty());
        for s in &result.samples {
            assert!(s.a > 0.0 && s.e < 1.0);
            assert!((0.0..TWO_PI).contains(&s.delta_varpi));
            assert!((0.0..TWO_PI).contains(&s.delta_omega_big));
        }
    }

    /// Reduced-scale grid: measured ~67 s wall / 12.9 CPU-min in release
    /// mode on 64 cores.
    /// Run with `cargo test --release -p p9-2021-orbit -- --ignored`.
    #[test]
    #[ignore]
    fn test_reduced_grid_runs() {
        let config = SimGridConfig::reduced();
        let results = run_grid(&config);
        assert_eq!(results.len(), 12);
        for r in &results {
            assert!(
                r.n_surviving > 0,
                "all particles lost at {:?} (seed {})",
                r.point,
                r.seed
            );
            assert!(!r.samples.is_empty());
        }
    }

    #[test]
    fn paper_grid_i9_not_aliased_with_e9() {
        use std::collections::HashSet;
        let pts = SimGridConfig::paper().grid_points();
        // The 121-point surrogate must exercise every (e9, i9) combination;
        // the old indexing produced only the 3 diagonal pairs.
        let pairs: HashSet<(u64, u64)> = pts
            .iter()
            .map(|p| ((p.e9 * 100.0) as u64, p.i9_deg as u64))
            .collect();
        assert!(
            pairs.len() >= 9,
            "only {} distinct (e9, i9) pairs in the Paper grid",
            pairs.len()
        );
        assert_eq!(pts.len(), 121);
    }
}
