//! Orbital pole analysis for TNO clustering.
//!
//! Pole vector: (i*cos(Omega), i*sin(Omega)) on the polar plane.
//! Stable objects are misaligned 5-10 deg from Laplace plane;
//! unstable objects are 0-2 deg from Laplace plane.

use p9_core::constants::GM_SUN;
use serde::{Deserialize, Serialize};

/// Orbital pole in the (i*cos(Omega), i*sin(Omega)) representation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrbitalPole {
    pub name: &'static str,
    /// i * cos(Omega) component (rad)
    pub x: f64,
    /// i * sin(Omega) component (rad)
    pub y: f64,
}

impl OrbitalPole {
    pub fn from_elements(name: &'static str, i: f64, omega_big: f64) -> Self {
        Self {
            name,
            x: i * omega_big.cos(),
            y: i * omega_big.sin(),
        }
    }

    /// Magnitude (inclination in rad).
    pub fn inclination(&self) -> f64 {
        (self.x * self.x + self.y * self.y).sqrt()
    }

    /// Position angle (longitude of ascending node in rad).
    pub fn omega_big(&self) -> f64 {
        self.y.atan2(self.x)
    }
}

/// Misalignment of an orbital pole from a reference direction (rad).
pub fn misalignment(pole: &OrbitalPole, ref_x: f64, ref_y: f64) -> f64 {
    let ref_mag = (ref_x * ref_x + ref_y * ref_y).sqrt();
    let pole_mag = pole.inclination();
    if ref_mag < 1e-10 || pole_mag < 1e-10 {
        return 0.0;
    }

    let cos_angle = (pole.x * ref_x + pole.y * ref_y) / (pole_mag * ref_mag);
    cos_angle.clamp(-1.0, 1.0).acos()
}

/// Compute the mean pole direction for a set of poles.
pub fn mean_pole(poles: &[OrbitalPole]) -> (f64, f64) {
    if poles.is_empty() {
        return (0.0, 0.0);
    }
    let n = poles.len() as f64;
    let mean_x = poles.iter().map(|p| p.x).sum::<f64>() / n;
    let mean_y = poles.iter().map(|p| p.y).sum::<f64>() / n;
    (mean_x, mean_y)
}

/// Scatter (standard deviation) of pole positions.
pub fn pole_scatter(poles: &[OrbitalPole]) -> f64 {
    if poles.len() < 2 {
        return 0.0;
    }
    let (mx, my) = mean_pole(poles);
    let n = poles.len() as f64;
    let var = poles
        .iter()
        .map(|p| (p.x - mx).powi(2) + (p.y - my).powi(2))
        .sum::<f64>()
        / n;
    var.sqrt()
}

/// Nodal precession rate from giant planets (rad/yr).
///
/// Omega_dot = -(3/4) * sqrt(GM_sun/a^3) * cos(i) / (1-e^2)^2 * sum_j(m_j/M_sun)(a_j/a)^2
pub fn nodal_precession_rate(a: f64, e: f64, i: f64) -> f64 {
    let n = (GM_SUN / (a * a * a)).sqrt(); // rad/day
    let eta_sq = (1.0 - e * e) * (1.0 - e * e);

    let giants: [(f64, f64); 4] = [
        (9.548e-4, 5.203),
        (2.858e-4, 9.537),
        (4.366e-5, 19.189),
        (5.151e-5, 30.070),
    ];

    let sum: f64 = giants.iter().map(|(m, aj)| m * (aj / a).powi(2)).sum();
    let rate_day = -0.75 * n * i.cos() / eta_sq * sum;
    rate_day * 365.25 // Convert to rad/yr
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_orbital_pole_from_elements() {
        let pole = OrbitalPole::from_elements("test", 30.0 * DEG2RAD, 0.0);
        assert!((pole.x - 30.0 * DEG2RAD).abs() < 1e-10);
        assert!(pole.y.abs() < 1e-10);
    }

    #[test]
    fn test_pole_inclination() {
        let pole = OrbitalPole::from_elements("test", 30.0 * DEG2RAD, 45.0 * DEG2RAD);
        assert!((pole.inclination() - 30.0 * DEG2RAD).abs() < 1e-10);
    }

    #[test]
    fn test_misalignment_zero() {
        let pole = OrbitalPole {
            name: "test",
            x: 0.3,
            y: 0.0,
        };
        let mis = misalignment(&pole, 0.3, 0.0);
        assert!(mis < 1e-10);
    }

    #[test]
    fn test_misalignment_perpendicular() {
        let pole = OrbitalPole {
            name: "test",
            x: 0.3,
            y: 0.0,
        };
        let mis = misalignment(&pole, 0.0, 0.3);
        assert!(
            (mis - std::f64::consts::FRAC_PI_2).abs() < 0.01,
            "misalignment = {} deg",
            mis / DEG2RAD
        );
    }

    #[test]
    fn test_mean_pole() {
        let poles = vec![
            OrbitalPole {
                name: "a",
                x: 0.1,
                y: 0.0,
            },
            OrbitalPole {
                name: "b",
                x: 0.3,
                y: 0.0,
            },
        ];
        let (mx, my) = mean_pole(&poles);
        assert!((mx - 0.2).abs() < 1e-10);
        assert!(my.abs() < 1e-10);
    }

    #[test]
    fn test_nodal_precession_negative() {
        // Prograde orbit (i < 90) should precess retrograde (negative rate)
        let rate = nodal_precession_rate(300.0, 0.7, 20.0 * DEG2RAD);
        assert!(
            rate < 0.0,
            "nodal precession rate should be negative for prograde"
        );
    }
}
