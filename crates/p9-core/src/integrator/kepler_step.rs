//! Kepler drift step: advance each body along its two-body Keplerian orbit.
//!
//! Uses the universal variable formulation (Battin 1999) for robust propagation
//! of elliptic, parabolic, and hyperbolic orbits.

use crate::types::StateVector;

/// Advance a state vector along its Keplerian orbit for time dt.
/// Uses Stumpff function formulation with Newton-Raphson iteration on the
/// universal variable s.
///
/// `gm` is the central body GM in AU^3/day^2.
pub fn kepler_drift(state: &StateVector, dt: f64, gm: f64) -> StateVector {
    let r0 = state.pos;
    let v0 = state.vel;
    let r0_mag = r0.norm();
    let v0_mag2 = v0.norm_squared();

    if r0_mag < 1e-30 {
        return *state;
    }

    // Specific orbital energy -> semi-major axis
    let energy = 0.5 * v0_mag2 - gm / r0_mag;
    let alpha = -2.0 * energy / gm; // = 1/a (positive for elliptic)

    let r0_dot_v0 = r0.dot(&v0);

    // Initial guess for universal variable s
    let mut s = if alpha > 0.0 {
        // Elliptic: s ~ n*dt where n = sqrt(gm * alpha^3)
        let n = (gm * alpha.powi(3)).sqrt();
        dt * n / (1.0 + 0.85 * (r0_dot_v0 / r0_mag).abs() * dt * n)
    } else if alpha < -1e-15 {
        // Hyperbolic
        let a = 1.0 / alpha;
        let sign_dt = dt.signum();
        sign_dt
            * (-a).sqrt()
            * ((-2.0 * gm * alpha * dt)
                / (r0_dot_v0 + sign_dt * (-gm / alpha).sqrt() * (1.0 - r0_mag * alpha)))
                .abs()
                .ln()
    } else {
        // Near-parabolic
        dt / r0_mag
    };

    // Newton-Raphson iteration on the universal Kepler equation
    // f(s) = r0*s + r0_dot_v0/sqrt(gm) * s^2 * c2(alpha*s^2)
    //        + (1 - r0*alpha) * s^3 * c3(alpha*s^2) - sqrt(gm)*dt = 0
    let sqrt_gm = gm.sqrt();

    for _ in 0..50 {
        let psi = alpha * s * s;
        let (c2, c3) = stumpff_c2c3(psi);

        let s2 = s * s;
        let s3 = s2 * s;

        let f_val = r0_mag * s + r0_dot_v0 / sqrt_gm * s2 * c2 + (1.0 - r0_mag * alpha) * s3 * c3
            - sqrt_gm * dt;

        // F'(χ) = r(χ) — the distance at the current iterate
        let df =
            r0_mag + r0_dot_v0 / sqrt_gm * s * (1.0 - psi * c3) + (1.0 - r0_mag * alpha) * s2 * c2;

        let ds = -f_val / df;
        s += ds;

        if ds.abs() < 1e-15 * s.abs().max(1.0) {
            break;
        }
    }

    // Compute f, g, fdot, gdot Lagrange coefficients
    let psi = alpha * s * s;
    let (c2, c3) = stumpff_c2c3(psi);

    let s2 = s * s;
    let s3 = s2 * s;

    let r_mag =
        r0_mag + r0_dot_v0 / sqrt_gm * s * (1.0 - psi * c3) + (1.0 - r0_mag * alpha) * s2 * c2;

    let f = 1.0 - s2 * c2 / r0_mag;
    let g = dt - s3 * c3 / sqrt_gm;

    let pos = f * r0 + g * v0;

    let f_dot = sqrt_gm / (r_mag * r0_mag) * s * (alpha * s2 * c3 - 1.0);
    let g_dot = 1.0 - s2 * c2 / r_mag;

    let vel = f_dot * r0 + g_dot * v0;

    StateVector::new(pos, vel)
}

/// Stumpff functions c2(psi) and c3(psi).
/// c2(psi) = (1 - cos(sqrt(psi))) / psi     for psi > 0
/// c3(psi) = (sqrt(psi) - sin(sqrt(psi))) / psi^(3/2)  for psi > 0
fn stumpff_c2c3(psi: f64) -> (f64, f64) {
    if psi > 1e-6 {
        let sqrt_psi = psi.sqrt();
        let c2 = (1.0 - sqrt_psi.cos()) / psi;
        let c3 = (sqrt_psi - sqrt_psi.sin()) / (psi * sqrt_psi);
        (c2, c3)
    } else if psi < -1e-6 {
        let sqrt_neg_psi = (-psi).sqrt();
        let c2 = (1.0 - sqrt_neg_psi.cosh()) / psi;
        let c3 = (sqrt_neg_psi.sinh() - sqrt_neg_psi) / (-psi * sqrt_neg_psi);
        (c2, c3)
    } else {
        // Taylor series for small psi
        (1.0 / 2.0 - psi / 24.0, 1.0 / 6.0 - psi / 120.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GM_SUN;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_circular_orbit() {
        // Circular orbit at 1 AU: v = sqrt(GM/r)
        let v_circ = (GM_SUN / 1.0_f64).sqrt();
        let state = StateVector::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, v_circ, 0.0));

        // One full period
        let period = 2.0 * std::f64::consts::PI / (GM_SUN).sqrt();
        let final_state = kepler_drift(&state, period, GM_SUN);

        assert_relative_eq!(final_state.pos.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(final_state.pos.y, 0.0, epsilon = 1e-10);
        assert_relative_eq!(final_state.pos.z, 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_energy_conservation() {
        // Eccentric orbit
        let state = StateVector::new(
            Vector3::new(1.0, 0.0, 0.0),
            Vector3::new(0.0, 0.02, 0.002), // AU/day
        );

        let energy =
            |s: &StateVector| -> f64 { 0.5 * s.vel.norm_squared() - GM_SUN / s.pos.norm() };

        let e0 = energy(&state);
        let mut s = state;
        for _ in 0..10_000 {
            s = kepler_drift(&s, 10.0, GM_SUN);
        }
        let e1 = energy(&s);

        assert_relative_eq!(e0, e1, epsilon = 1e-12);
    }

    #[test]
    fn test_half_orbit() {
        // Start at perihelion of e=0.5 orbit at 1 AU perihelion
        // a = q/(1-e) = 2 AU
        let a: f64 = 2.0;
        let e = 0.5;
        let v_peri = (GM_SUN * (1.0 + e) / (a * (1.0 - e))).sqrt();
        let state = StateVector::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, v_peri, 0.0));

        let period = 2.0 * std::f64::consts::PI * (a.powi(3) / GM_SUN).sqrt();
        let half = kepler_drift(&state, period / 2.0, GM_SUN);

        // At aphelion: x should be -Q = -a(1+e) = -3
        assert_relative_eq!(half.pos.x, -3.0, epsilon = 1e-8);
        assert_relative_eq!(half.pos.y, 0.0, epsilon = 1e-8);
    }
}
