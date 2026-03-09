//! 3D inclination exploration from Brown & Batygin (2016) Section 3.
//!
//! Surveys the (i₉, ω₉) parameter space at fixed a₉=700 AU, e₉=0.6, m=10 M_Earth.
//! - i₉ ∈ {1, 10, 20, 30, 60, 90, 120, 150°}
//! - ω₉ ∈ (0-360°, 30° steps)
//!
//! For each combination, runs a scattered disk simulation and evaluates
//! whether the resulting population matches the observed KBO clustering.

use p9_core::constants::*;
use p9_core::types::P9Params;

/// Inclination survey grid point.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct InclinationGridPoint {
    /// Planet Nine inclination (radians)
    pub i_p9: f64,
    /// Planet Nine argument of perihelion (radians)
    pub omega_p9: f64,
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

/// Acceptance criteria for inclination survey:
/// 1. Mean pole angle > 20°
/// 2. RMS spread < 6.2° (matching observed 6 KBOs)
/// 3. Confinement probability > 0.5
pub fn evaluate_inclination_acceptance(
    pole_angle_mean: f64,
    pole_angle_rms: f64,
    confinement_prob: f64,
) -> bool {
    pole_angle_mean > 20.0 && pole_angle_rms < 6.2 && confinement_prob > 0.5
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
    fn test_acceptance_criteria() {
        assert!(evaluate_inclination_acceptance(25.0, 5.0, 0.7));
        assert!(!evaluate_inclination_acceptance(15.0, 5.0, 0.7)); // pole angle too low
        assert!(!evaluate_inclination_acceptance(25.0, 8.0, 0.7)); // spread too high
        assert!(!evaluate_inclination_acceptance(25.0, 5.0, 0.3)); // confinement too low
    }
}
