//! Planet Nine heliocentric geometry along the Brown & Batygin orbit.
//!
//! The Cassini constraint is a function of WHERE P9 is on its orbit, i.e. of
//! its true anomaly `ν`. This module turns a `ν` into the full heliocentric
//! position vector (and distance) using `p9_core`'s element→Cartesian
//! conversion, reusing the workspace `P9Params` and the shared Kepler solver.

use nalgebra::Vector3;

use p9_core::constants::{GM_SUN, TWO_PI};
use p9_core::types::{elements_to_cartesian, OrbitalElements, P9Params, StateVector};

/// Heliocentric geometry of Planet Nine at one true anomaly.
#[derive(Debug, Clone, Copy)]
pub struct P9Geometry {
    /// True anomaly along the orbit (radians).
    pub true_anomaly: f64,
    /// Heliocentric position (AU, ecliptic J2000).
    pub position: Vector3<f64>,
    /// Heliocentric distance `r(ν) = a(1-e²)/(1 + e cos ν)` (AU).
    pub distance: f64,
}

/// Convert a true anomaly to the mean anomaly for an elliptic orbit
/// (`0 ≤ e < 1`). Inverse of the Kepler chain `M → E → ν`.
pub fn true_to_mean_anomaly(e: f64, nu: f64) -> f64 {
    // Eccentric anomaly E from true anomaly ν.
    let ea = ((1.0 - e) / (1.0 + e)).sqrt() * (nu / 2.0).tan();
    let ea = 2.0 * ea.atan();
    // Mean anomaly M = E - e sin E, wrapped to [0, 2π).
    let m = ea - e * ea.sin();
    m.rem_euclid(TWO_PI)
}

/// Heliocentric position of P9 at true anomaly `nu` (radians) on the orbit
/// described by `params`. The `mean_anomaly` field of `params` is overridden:
/// the orbit's orientation (a, e, i, ω, Ω) is what matters, and we sweep ν.
pub fn p9_position_at_true_anomaly(params: &P9Params, nu: f64) -> P9Geometry {
    let elements = OrbitalElements {
        a: params.a,
        e: params.e,
        i: params.i,
        omega_big: params.omega_big,
        omega: params.omega,
        mean_anomaly: true_to_mean_anomaly(params.e, nu),
    };
    let state: StateVector = elements_to_cartesian(&elements, GM_SUN);
    let position = state.pos;
    P9Geometry {
        true_anomaly: nu,
        position,
        distance: position.norm(),
    }
}

/// Analytic heliocentric distance `r(ν)` for the conic (AU), independent of the
/// Cartesian path; used as a cross-check.
pub fn helio_distance_at_true_anomaly(params: &P9Params, nu: f64) -> f64 {
    let p = params.a * (1.0 - params.e * params.e);
    p / (1.0 + params.e * nu.cos())
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn bb_orbit() -> P9Params {
        crate::brown_batygin_orbit()
    }

    #[test]
    fn true_mean_roundtrip_through_core_kepler() {
        let e = 0.6;
        for deg in [10.0, 73.0, 108.0, 200.0, 333.0] {
            let nu = deg * PI / 180.0;
            let m = true_to_mean_anomaly(e, nu);
            // Re-solve M → E → ν via the shared core solver and compare.
            let ea = p9_core::types::solve_kepler(e, m);
            let sin_nu = (1.0 - e * e).sqrt() * ea.sin() / (1.0 - e * ea.cos());
            let cos_nu = (ea.cos() - e) / (1.0 - e * ea.cos());
            let nu_back = sin_nu.atan2(cos_nu).rem_euclid(TWO_PI);
            assert_relative_eq!(nu_back, nu.rem_euclid(TWO_PI), epsilon = 1e-9);
        }
    }

    #[test]
    fn perihelion_and_aphelion_distances() {
        let p = bb_orbit();
        let peri = p9_position_at_true_anomaly(&p, 0.0);
        let apo = p9_position_at_true_anomaly(&p, PI);
        assert_relative_eq!(peri.distance, p.a * (1.0 - p.e), epsilon = 1e-9);
        assert_relative_eq!(apo.distance, p.a * (1.0 + p.e), epsilon = 1e-9);
    }

    #[test]
    fn cartesian_distance_matches_conic_formula() {
        let p = bb_orbit();
        for deg in [0.0, 45.0, 108.0, 180.0, 270.0] {
            let nu = deg * PI / 180.0;
            let g = p9_position_at_true_anomaly(&p, nu);
            let r_conic = helio_distance_at_true_anomaly(&p, nu);
            assert_relative_eq!(g.distance, r_conic, epsilon = 1e-7);
        }
    }
}
