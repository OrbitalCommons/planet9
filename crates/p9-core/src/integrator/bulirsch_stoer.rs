//! Bulirsch-Stoer integrator for high-accuracy integration during close encounters.
//!
//! Uses the modified midpoint method with Richardson extrapolation to achieve
//! arbitrary accuracy. This is the fallback integrator when a particle enters
//! a planet's Hill sphere and the symplectic integrator would fail.
//!
//! Reference: Press et al. Numerical Recipes, Stoer & Bulirsch (1980)

use nalgebra::Vector3;

use crate::constants::GM_SUN;
use crate::types::{MassiveBody, StateVector};

/// Maximum subdivision sequence for Bulirsch-Stoer extrapolation.
/// Bulirsch-Stoer sequence: 2, 4, 6, 8, 12, 16, 24, 32, ...
const BS_SEQUENCE: [usize; 12] = [2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128];

/// Maximum number of columns in the extrapolation tableau.
const MAX_COLUMNS: usize = 12;

/// Compute gravitational acceleration on a test particle from the Sun + massive bodies.
fn acceleration(pos: &Vector3<f64>, bodies: &[MassiveBody]) -> Vector3<f64> {
    let r = pos.norm();
    let mut accel = if r > 1e-30 {
        -GM_SUN * pos / (r * r * r)
    } else {
        Vector3::zeros()
    };

    for body in bodies {
        let dr = pos - body.state.pos;
        let dr_mag = dr.norm();
        if dr_mag > 1e-30 {
            accel -= body.gm * dr / (dr_mag * dr_mag * dr_mag);
        }
    }

    accel
}

/// Modified midpoint method: integrate from t to t+H using n substeps.
/// Returns the final state vector.
fn modified_midpoint(
    state: &StateVector,
    bodies: &[MassiveBody],
    h_total: f64,
    n_steps: usize,
) -> StateVector {
    let h = h_total / n_steps as f64;

    // First step: Euler
    let a0 = acceleration(&state.pos, bodies);
    let mut z_prev = state.pos;
    let mut z_curr = state.pos + h * state.vel;
    let mut v_prev = state.vel;
    let mut v_curr = state.vel + h * a0;

    // General steps: leapfrog
    for _ in 1..n_steps {
        let a = acceleration(&z_curr, bodies);
        let z_next = z_prev + 2.0 * h * v_curr;
        let v_next = v_prev + 2.0 * h * a;

        z_prev = z_curr;
        z_curr = z_next;
        v_prev = v_curr;
        v_curr = v_next;
    }

    // Final smoothing step
    let a_final = acceleration(&z_curr, bodies);
    let pos = 0.5 * (z_prev + z_curr + h * v_curr);
    let vel = 0.5 * (v_prev + v_curr + h * a_final);

    StateVector::new(pos, vel)
}

/// Perform one Bulirsch-Stoer step on a single test particle.
///
/// Adaptively subdivides the step until the extrapolation converges to within
/// the specified tolerance `epsilon`.
///
/// Returns (new_state, actual_dt_used). If the step cannot converge,
/// it returns the best estimate with the full dt.
pub fn bs_step(
    state: &StateVector,
    bodies: &[MassiveBody],
    dt: f64,
    epsilon: f64,
) -> (StateVector, f64) {
    // Extrapolation tableau (rational or polynomial)
    // We use polynomial extrapolation for simplicity.
    let mut pos_tableau = [[Vector3::<f64>::zeros(); MAX_COLUMNS]; MAX_COLUMNS];
    let mut vel_tableau = [[Vector3::<f64>::zeros(); MAX_COLUMNS]; MAX_COLUMNS];

    for k in 0..MAX_COLUMNS {
        let n = BS_SEQUENCE[k];
        let result = modified_midpoint(state, bodies, dt, n);

        pos_tableau[k][0] = result.pos;
        vel_tableau[k][0] = result.vel;

        if k > 0 {
            // Richardson extrapolation (must go forward: each j depends on j-1)
            for j in 1..=k {
                let ratio = (BS_SEQUENCE[k] as f64 / BS_SEQUENCE[k - j] as f64).powi(2);
                let factor = 1.0 / (ratio - 1.0);

                pos_tableau[k][j] = pos_tableau[k][j - 1]
                    + factor * (pos_tableau[k][j - 1] - pos_tableau[k - 1][j - 1]);
                vel_tableau[k][j] = vel_tableau[k][j - 1]
                    + factor * (vel_tableau[k][j - 1] - vel_tableau[k - 1][j - 1]);
            }

            // Check convergence
            let pos_err = (pos_tableau[k][k] - pos_tableau[k][k - 1]).norm();
            let vel_err = (vel_tableau[k][k] - vel_tableau[k][k - 1]).norm();

            let pos_scale = pos_tableau[k][k].norm().max(1.0);
            let vel_scale = vel_tableau[k][k].norm().max(1e-6);

            if pos_err / pos_scale < epsilon && vel_err / vel_scale < epsilon {
                return (StateVector::new(pos_tableau[k][k], vel_tableau[k][k]), dt);
            }
        }
    }

    // Did not converge — return best estimate
    let k = MAX_COLUMNS - 1;
    (StateVector::new(pos_tableau[k][k], vel_tableau[k][k]), dt)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_bs_circular_orbit() {
        let a = 1.0;
        let v = (GM_SUN / a).sqrt();
        let state = StateVector::new(Vector3::new(a, 0.0, 0.0), Vector3::new(0.0, v, 0.0));

        let period = 2.0 * std::f64::consts::PI * (a * a * a / GM_SUN).sqrt();
        let n_steps = 100;
        let dt = period / n_steps as f64;

        let mut s = state;
        for _ in 0..n_steps {
            let (new_s, _) = bs_step(&s, &[], dt, 1e-12);
            s = new_s;
        }

        // Should return to start within high accuracy
        assert_relative_eq!(s.pos.x, a, epsilon = 1e-8);
        assert_relative_eq!(s.pos.y, 0.0, epsilon = 1e-8);
    }

    #[test]
    fn test_bs_energy_conservation() {
        let state = StateVector::new(Vector3::new(1.0, 0.0, 0.0), Vector3::new(0.0, 0.02, 0.002));

        let energy = |s: &StateVector| 0.5 * s.vel.norm_squared() - GM_SUN / s.pos.norm();

        let e0 = energy(&state);
        let mut s = state;
        for _ in 0..100 {
            let (new_s, _) = bs_step(&s, &[], 50.0, 1e-12);
            s = new_s;
        }
        let e1 = energy(&s);

        let de = ((e1 - e0) / e0).abs();
        assert!(de < 1e-10, "BS energy error: {:.2e}", de);
    }
}
