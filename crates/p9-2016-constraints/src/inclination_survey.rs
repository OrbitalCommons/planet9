//! 3D inclination exploration from Brown & Batygin (2016) Section 3.
//!
//! Surveys the (i₉, ω₉) parameter space at fixed a₉=700 AU, e₉=0.6, m=10 M_Earth.
//! - i₉ ∈ {1, 10, 20, 30, 60, 90, 120, 150°}
//! - ω₉ ∈ (0-360°, 30° steps)
//!
//! For each combination, [`run_inclination_point`] runs the inclined
//! scattered-disk simulation (p9-core `scattered_disk_sim`, giants as the
//! J2-averaged quadrupole, P9 direct) and evaluates whether the surviving
//! population's orbital-pole statistics match the observed KBO sample.

use p9_core::constants::*;
use p9_core::initial_conditions::scattered_disk_sim::{run_scattered_disk, DiskSimConfig};
use p9_core::types::P9Params;
use p9_core::units::{radians, Angle};
use rayon::prelude::*;

use crate::parameter_grid::SurveyRunConfig;

/// Inclination survey grid point.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct InclinationGridPoint {
    /// Planet Nine inclination (radians)
    pub i_p9: f64,
    /// Planet Nine argument of perihelion (radians)
    pub omega_p9: f64,
}

impl InclinationGridPoint {
    /// Planet Nine inclination as a typed [`Angle`].
    pub fn inclination(&self) -> Angle {
        radians(self.i_p9)
    }

    /// Planet Nine argument of perihelion as a typed [`Angle`].
    pub fn argument_of_perihelion(&self) -> Angle {
        radians(self.omega_p9)
    }
}

/// Generate the inclination survey grid.
pub fn inclination_grid() -> Vec<InclinationGridPoint> {
    let inclinations_deg = [1.0, 10.0, 20.0, 30.0, 60.0, 90.0, 120.0, 150.0];
    let omega_steps = 12; // 0 to 360 in 30 degree steps

    let mut grid = Vec::new();

    for &i_deg in &inclinations_deg {
        for j in 0..omega_steps {
            let omega_deg = j as f64 * 30.0;
            grid.push(InclinationGridPoint {
                i_p9: i_deg * DEG2RAD,
                omega_p9: omega_deg * DEG2RAD,
            });
        }
    }

    grid
}

/// Convert an inclination grid point to P9 parameters.
pub fn grid_point_to_p9(point: &InclinationGridPoint) -> P9Params {
    P9Params {
        mass_earth: 10.0,
        a: 700.0,
        e: 0.6,
        i: point.i_p9,
        omega: point.omega_p9,
        omega_big: 100.0 * DEG2RAD,
        mean_anomaly: 0.0,
    }
}

/// Result of an inclination survey evaluation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct InclinationResult {
    pub point: InclinationGridPoint,
    /// Pole angle dispersion of the surviving particles (degrees)
    pub pole_angle_mean: f64,
    /// RMS spread of pole angles (degrees)
    pub pole_angle_rms: f64,
    /// Confinement probability
    pub confinement_prob: f64,
    /// Whether this point is accepted
    pub accepted: bool,
}

/// Observed orbital-pole statistics of the vetted ETNO sample
/// (`p9_core::data::etno`): returns (mean pole tilt from the ecliptic pole,
/// RMS angular scatter of the poles about their mean direction), both in
/// degrees. These are the comparison values for the inclination survey,
/// replacing previously invented thresholds (20°, 6.2°).
pub fn observed_pole_stats() -> (f64, f64) {
    let poles: Vec<(f64, f64, f64)> = p9_core::data::etno::BROWN_2017_SAMPLE
        .iter()
        .map(|k| pole_direction(k.i_deg * DEG2RAD, k.omega_big_deg * DEG2RAD))
        .collect();

    let mut mean = (0.0, 0.0, 0.0);
    for p in &poles {
        mean.0 += p.0;
        mean.1 += p.1;
        mean.2 += p.2;
    }
    let mag = (mean.0 * mean.0 + mean.1 * mean.1 + mean.2 * mean.2).sqrt();
    let mean = (mean.0 / mag, mean.1 / mag, mean.2 / mag);

    let tilt = pole_separation_deg(mean, (0.0, 0.0, 1.0));
    let rms = (poles
        .iter()
        .map(|&p| pole_separation_deg(mean, p).powi(2))
        .sum::<f64>()
        / poles.len() as f64)
        .sqrt();

    (tilt, rms)
}

/// Acceptance criteria for the inclination survey, recast as a comparison
/// against the observed ETNO sample (`observed_pole_stats`; tilt ≈ 13.8°,
/// RMS ≈ 15.8°):
/// 1. the simulated mean pole tilt is consistent with the observed mean
///    tilt to within the observed sample scatter;
/// 2. the simulated poles scatter no more than the observed RMS (the
///    population is at least as confined as the data);
/// 3. confinement probability > 0.5 — **assumption**: a majority floor not
///    published in the paper.
pub fn evaluate_inclination_acceptance(
    pole_angle_mean: f64,
    pole_angle_rms: f64,
    confinement_prob: f64,
) -> bool {
    let (obs_tilt, obs_rms) = observed_pole_stats();
    (pole_angle_mean - obs_tilt).abs() <= obs_rms
        && pole_angle_rms <= obs_rms
        && confinement_prob > 0.5
}

/// Run one (i₉, ω₉) survey point: inclined scattered-disk simulation with the
/// grid point's Planet Nine, then orbital-pole statistics of the surviving
/// distant population against the observed ETNO comparison values.
///
/// `confinement_prob` is the anti-aligned Δϖ fraction of the qualifying
/// survivors — the same confinement statistic the planar Section-2 survey
/// uses, evaluated here alongside the 3D pole criteria.
pub fn run_inclination_point(
    point: &InclinationGridPoint,
    run: &SurveyRunConfig,
    seed: u64,
) -> InclinationResult {
    let mut config = DiskSimConfig::inclined_nominal();
    config.p9 = grid_point_to_p9(point);
    config.n_particles = run.n_particles;
    config.t_total = run.t_total;
    config.dt = run.dt;
    config.snapshot_interval = run.t_total;

    let snapshots = run_scattered_disk(&config, seed);
    let last = snapshots.last().expect("at least the initial snapshot");

    // Distant survivors in the clustering band.
    let survivors: Vec<_> = last
        .elements
        .iter()
        .filter(|el| el.a >= run.a_min_cluster && el.a <= run.a_max_cluster)
        .collect();

    if survivors.is_empty() {
        return InclinationResult {
            point: *point,
            pole_angle_mean: f64::NAN,
            pole_angle_rms: f64::NAN,
            confinement_prob: 0.0,
            accepted: false,
        };
    }

    // Mean pole direction and scatter of the surviving population.
    let poles: Vec<(f64, f64, f64)> = survivors
        .iter()
        .map(|el| pole_direction(el.i, el.omega_big))
        .collect();
    let mut mean = (0.0, 0.0, 0.0);
    for p in &poles {
        mean.0 += p.0;
        mean.1 += p.1;
        mean.2 += p.2;
    }
    let mag = (mean.0 * mean.0 + mean.1 * mean.1 + mean.2 * mean.2)
        .sqrt()
        .max(1e-12);
    let mean = (mean.0 / mag, mean.1 / mag, mean.2 / mag);
    let pole_angle_mean = pole_separation_deg(mean, (0.0, 0.0, 1.0));
    let pole_angle_rms = (poles
        .iter()
        .map(|&p| pole_separation_deg(mean, p).powi(2))
        .sum::<f64>()
        / poles.len() as f64)
        .sqrt();

    // Anti-aligned Δϖ confinement fraction among the survivors.
    let anti = survivors
        .iter()
        .filter(|el| {
            let dv = (el.omega + el.omega_big - last.varpi_p9).rem_euclid(std::f64::consts::TAU);
            (dv - std::f64::consts::PI).abs() < std::f64::consts::FRAC_PI_2
        })
        .count();
    let confinement_prob = anti as f64 / survivors.len() as f64;

    InclinationResult {
        point: *point,
        pole_angle_mean,
        pole_angle_rms,
        confinement_prob,
        accepted: evaluate_inclination_acceptance(
            pole_angle_mean,
            pole_angle_rms,
            confinement_prob,
        ),
    }
}

/// Run the full (i₉, ω₉) survey in parallel with deterministic per-point seeds.
pub fn run_inclination_grid(
    grid: &[InclinationGridPoint],
    run: &SurveyRunConfig,
    base_seed: u64,
) -> Vec<InclinationResult> {
    grid.par_iter()
        .enumerate()
        .map(|(k, point)| run_inclination_point(point, run, base_seed.wrapping_add(k as u64)))
        .collect()
}

/// Compute the pole angle of an orbit relative to the ecliptic.
///
/// The orbital pole is at (sin(i)sin(Ω), -sin(i)cos(Ω), cos(i)) in ecliptic coords.
/// The pole angle is arccos(cos(i)) = i for inclination.
/// But the paper's "pole angle" refers to the angular distance between the
/// orbital pole and a reference direction.
pub fn pole_direction(i: f64, omega_big: f64) -> (f64, f64, f64) {
    let sin_i = i.sin();
    (sin_i * omega_big.sin(), -sin_i * omega_big.cos(), i.cos())
}

/// Compute angular separation between two pole directions (degrees).
pub fn pole_separation_deg(pole1: (f64, f64, f64), pole2: (f64, f64, f64)) -> f64 {
    let dot = pole1.0 * pole2.0 + pole1.1 * pole2.1 + pole1.2 * pole2.2;
    dot.clamp(-1.0, 1.0).acos() * RAD2DEG
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reduced_scale_runner_produces_populated_results() {
        // The Section-3 survey now has a real runner: a low-inclination P9
        // at the nominal orientation, run at reduced scale, yields finite
        // pole statistics from actual surviving particles (this was the
        // "no runner / InclinationResult never constructed" gap).
        use crate::parameter_grid::SurveyRunConfig;
        let point = InclinationGridPoint {
            i_p9: 20.0 * DEG2RAD,
            omega_p9: 150.0 * DEG2RAD,
        };
        let run = SurveyRunConfig::reduced();
        let res = run_inclination_point(&point, &run, 2016);
        assert!(
            res.pole_angle_mean.is_finite(),
            "no survivors at reduced scale"
        );
        assert!(res.pole_angle_rms >= 0.0);
        assert!((0.0..=1.0).contains(&res.confinement_prob));
        // Acceptance is a computed verdict either way, not a constant.
        let (_obs_tilt, obs_rms) = observed_pole_stats();
        assert!(obs_rms > 0.0);
    }

    #[test]
    fn test_inclination_grid_size() {
        let grid = inclination_grid();
        // 8 inclinations × 12 omega values = 96
        assert_eq!(grid.len(), 96);
    }

    #[test]
    fn test_grid_point_conversion() {
        let point = InclinationGridPoint {
            i_p9: 30.0 * DEG2RAD,
            omega_p9: 150.0 * DEG2RAD,
        };
        let p9 = grid_point_to_p9(&point);

        assert!((p9.i - 30.0 * DEG2RAD).abs() < 1e-10);
        assert!((p9.omega - 150.0 * DEG2RAD).abs() < 1e-10);
        assert!((p9.mass_earth - 10.0).abs() < 0.1);
    }

    #[test]
    fn test_grid_point_typed_angles_match_f64() {
        use approx::assert_relative_eq;
        let point = InclinationGridPoint {
            i_p9: 30.0 * DEG2RAD,
            omega_p9: 150.0 * DEG2RAD,
        };
        assert_relative_eq!(
            (point.inclination() / radians(1.0)).value,
            point.i_p9,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (point.argument_of_perihelion() / radians(1.0)).value,
            point.omega_p9,
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_pole_direction() {
        // Ecliptic orbit: pole at (0, 0, 1)
        let (x, y, z) = pole_direction(0.0, 0.0);
        assert!((z - 1.0).abs() < 1e-10);
        assert!(x.abs() < 1e-10);
        assert!(y.abs() < 1e-10);
    }

    #[test]
    fn test_pole_separation() {
        let pole1 = pole_direction(0.0, 0.0);
        let pole2 = pole_direction(30.0 * DEG2RAD, 0.0);
        let sep = pole_separation_deg(pole1, pole2);
        assert!(
            (sep - 30.0).abs() < 0.1,
            "Separation should be ~30°: {:.1}",
            sep
        );
    }

    #[test]
    fn test_observed_pole_stats() {
        // Mean pole tilt ≈ 13.8° (sample inclinations 12–30° with widely
        // spread nodes) and RMS scatter ≈ 15.8°.
        let (tilt, rms) = observed_pole_stats();
        assert!((tilt - 13.8).abs() < 0.5, "mean tilt = {tilt:.1}°");
        assert!((rms - 15.8).abs() < 0.5, "pole RMS = {rms:.1}°");
    }

    #[test]
    fn test_acceptance_criteria() {
        let (obs_tilt, obs_rms) = observed_pole_stats();
        // Matching the observed sample: accepted
        assert!(evaluate_inclination_acceptance(obs_tilt, obs_rms, 0.7));
        // Mean tilt inconsistent with the observed plane: rejected
        assert!(!evaluate_inclination_acceptance(
            obs_tilt + obs_rms + 1.0,
            obs_rms,
            0.7
        ));
        // Scatter larger than observed: rejected
        assert!(!evaluate_inclination_acceptance(
            obs_tilt,
            obs_rms + 1.0,
            0.7
        ));
        // Confinement below the (assumed) majority floor: rejected
        assert!(!evaluate_inclination_acceptance(obs_tilt, obs_rms, 0.3));
    }
}
