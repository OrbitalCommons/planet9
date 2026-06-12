//! Orbital pole analysis for TNO clustering.
//!
//! A pole is represented by (i, Ω); the (i·cosΩ, i·sinΩ) plane is kept
//! for plotting, but misalignments are computed as true spherical angular
//! distances between pole unit vectors
//!
//!   n̂ = (sin i sin Ω, −sin i cos Ω, cos i),
//!
//! (via p9-core's `coords::sky::angular_distance`, i.e. starfield's
//! `Cartesian3::angular_distance`), not as 2D distances in the
//! (i cosΩ, i sinΩ) plane (which conflates a 90° node difference at i = 17°
//! with a 90° pole separation). The reference Laplace plane is *computed*
//! from the giant planets' total orbital angular momentum (J2000 states from
//! p9-core) — at ETNO distances the Laplace plane asymptotes to the
//! invariable plane.

use nalgebra::Vector3;
use p9_core::constants::GM_SUN;
use p9_core::coords::sky::angular_distance;
use p9_core::initial_conditions::planets::{giant_planets_at, giant_planets_j2000, Time};
use p9_core::types::{cartesian_to_elements, MassiveBody};
use serde::{Deserialize, Serialize};

/// Orbital pole, stored as inclination + node with the planar projection
/// kept for plotting.
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

    /// Unit vector of the orbit normal in ecliptic coordinates:
    /// n̂ = (sin i sin Ω, −sin i cos Ω, cos i).
    pub fn unit_vector(&self) -> Vector3<f64> {
        let i = self.inclination();
        let node = self.omega_big();
        Vector3::new(i.sin() * node.sin(), -i.sin() * node.cos(), i.cos())
    }
}

/// True spherical angular distance between two orbital poles (rad):
///
///   cos d = cos i₁ cos i₂ + sin i₁ sin i₂ cos(Ω₁ − Ω₂)
pub fn misalignment(pole: &OrbitalPole, reference: &OrbitalPole) -> f64 {
    angular_distance(&pole.unit_vector(), &reference.unit_vector())
}

/// Pole of the plane normal to the total orbital angular momentum
/// L = Σ mᵢ rᵢ×vᵢ of the given bodies (heliocentric states).
fn angular_momentum_pole(bodies: &[MassiveBody]) -> OrbitalPole {
    let mut l = Vector3::zeros();
    for body in bodies {
        l += body.mass * body.state.pos.cross(&body.state.vel);
    }
    let n = l.normalize();
    let i = n.z.clamp(-1.0, 1.0).acos();
    // n̂ = (sin i sin Ω, −sin i cos Ω, cos i) → Ω = atan2(n_x, −n_y)
    let omega_big = n.x.atan2(-n.y);
    OrbitalPole::from_elements("Laplace plane", i, omega_big)
}

/// Reference pole of the giant-planet (Laplace/invariable) plane,
/// computed from the total orbital angular momentum L = Σ mᵢ rᵢ×vᵢ of
/// the four giants at J2000 (p9-core DE440 states). Inclination to the
/// ecliptic ≈ 1.58°.
pub fn laplace_pole() -> OrbitalPole {
    angular_momentum_pole(&giant_planets_j2000())
}

/// Ephemeris-backed variant of [`laplace_pole`]: the same angular-momentum
/// pole, but from DE-kernel giant-planet states at an arbitrary epoch.
/// Because the giants carry ~98% of the planetary angular momentum, the
/// pole is epoch-stable to a small fraction of a degree.
///
/// Requires a cached DE kernel (never downloads); `Err` on a hermetic
/// machine — fall back to [`laplace_pole`].
pub fn laplace_pole_at(t: &Time) -> Result<OrbitalPole, String> {
    Ok(angular_momentum_pole(&giant_planets_at(t)?))
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

/// Nodal precession rate from the giant planets (rad/yr):
///
/// Omega_dot = -(3/4) √(GM_sun/a³) cos(i)/(1−e²)² Σ_j (m_j/M_sun)(a_j/a)²
///
/// Giant-planet masses and semi-major axes come from the p9-core J2000
/// states (no local table).
pub fn nodal_precession_rate(a: f64, e: f64, i: f64) -> f64 {
    let n = (GM_SUN / (a * a * a)).sqrt(); // rad/day
    let eta_sq = (1.0 - e * e) * (1.0 - e * e);

    let sum: f64 = giant_planets_j2000()
        .iter()
        .map(|body| {
            let a_j = cartesian_to_elements(&body.state, GM_SUN).a;
            body.mass * (a_j / a).powi(2)
        })
        .sum();

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
    fn test_unit_vector_zero_inclination_is_z() {
        let pole = OrbitalPole::from_elements("flat", 0.0, 1.234);
        let n = pole.unit_vector();
        assert!(n[0].abs() < 1e-12 && n[1].abs() < 1e-12);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_misalignment_zero() {
        let p1 = OrbitalPole::from_elements("a", 0.3, 0.7);
        let p2 = OrbitalPole::from_elements("b", 0.3, 0.7);
        assert!(misalignment(&p1, &p2) < 1e-10);
    }

    #[test]
    fn test_misalignment_spherical_not_planar() {
        // i = 0.3 rad, nodes 90° apart: the planar (i cosΩ, i sinΩ)
        // metric wrongly gives 90°; the true spherical distance is
        // acos(cos²i) ≈ 0.4215 rad ≈ 24.2°.
        let p1 = OrbitalPole::from_elements("a", 0.3, 0.0);
        let p2 = OrbitalPole::from_elements("b", 0.3, std::f64::consts::FRAC_PI_2);
        let mis = misalignment(&p1, &p2);
        let expected = (0.3_f64.cos().powi(2)).acos();
        assert!(
            (mis - expected).abs() < 1e-9,
            "misalignment = {:.4} rad, expected {:.4}",
            mis,
            expected
        );
        assert!(mis < std::f64::consts::FRAC_PI_2 * 0.5);
    }

    #[test]
    fn test_misalignment_antipodal_inclinations() {
        // Same node, inclinations 10° and 30° → distance exactly 20°.
        let p1 = OrbitalPole::from_elements("a", 10.0 * DEG2RAD, 1.0);
        let p2 = OrbitalPole::from_elements("b", 30.0 * DEG2RAD, 1.0);
        let mis = misalignment(&p1, &p2);
        assert!((mis - 20.0 * DEG2RAD).abs() < 1e-9);
    }

    #[test]
    fn test_laplace_pole_near_invariable_plane() {
        // The giant-planet angular-momentum plane is inclined ≈ 1.58° to
        // the ecliptic.
        let lap = laplace_pole();
        let i_deg = lap.inclination() / DEG2RAD;
        assert!(
            i_deg > 1.0 && i_deg < 2.2,
            "Laplace-plane inclination = {i_deg:.2}°, expected ≈ 1.58°"
        );
    }

    #[test]
    fn test_laplace_pole_matches_ephemeris_variant() {
        use p9_core::initial_conditions::planets::Timescale;

        let ts = Timescale::default();
        // J2000: identical inputs up to kernel choice (de421 fallback
        // differs from the pinned DE440 states by ≲ 3.5e-6 AU).
        let at_j2000 = match laplace_pole_at(&ts.tdb_jd(p9_core::constants::J2000)) {
            Ok(p) => p,
            Err(_) => {
                eprintln!("skipping: no DE kernel in starfield cache");
                return;
            }
        };
        let pinned = laplace_pole();
        let sep = misalignment(&pinned, &at_j2000) / DEG2RAD;
        assert!(sep < 0.01, "J2000 pole separation = {sep:.4}°");

        // Other epochs: the invariable-plane direction is epoch-stable to
        // a small fraction of a degree (osculating exchanges only).
        let at_2017 = laplace_pole_at(&ts.utc((2017, 1, 1))).unwrap();
        let sep = misalignment(&pinned, &at_2017) / DEG2RAD;
        assert!(sep < 0.2, "2017 pole separation = {sep:.4}°");
        let i_deg = at_2017.inclination() / DEG2RAD;
        assert!(i_deg > 1.0 && i_deg < 2.2, "inclination = {i_deg:.2}°");
    }

    #[test]
    fn test_misalignment_from_laplace_plane() {
        // A typical ETNO pole (i ≈ 20°) sits ~18–22° from the computed
        // Laplace plane, depending on the node.
        let lap = laplace_pole();
        let etno = OrbitalPole::from_elements("etno", 20.0 * DEG2RAD, 2.0);
        let mis = misalignment(&etno, &lap) / DEG2RAD;
        assert!(mis > 17.0 && mis < 23.0, "misalignment = {mis:.1}°");
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
