//! Interaction Hamiltonian kick: apply gravitational perturbation velocity changes.
//!
//! In the Wisdom-Holman splitting, the kick step computes accelerations from the
//! interaction Hamiltonian (everything not captured by the Kepler step) and
//! applies velocity impulses for dt/2 (leapfrog).
//!
//! The interaction Hamiltonian includes:
//! - Planet-planet gravitational interactions
//! - Planet-particle gravitational interactions
//! - Indirect term (barycentric correction)
//! - Optionally: J2/J4 oblateness

use nalgebra::Vector3;

use crate::types::{MassiveBody, StateVector};

/// Compute the gravitational acceleration on a body at position `r` due to
/// a perturbing body with GM `gm_pert` at position `r_pert`.
///
/// This is the direct + indirect term:
///   a = -gm_pert * [ (r - r_pert)/|r - r_pert|^3 + r_pert/|r_pert|^3 ]
///
/// The second term is the indirect (Sun-centered) acceleration.
#[inline]
pub fn perturbation_acceleration(
    r: &Vector3<f64>,
    r_pert: &Vector3<f64>,
    gm_pert: f64,
) -> Vector3<f64> {
    let dr = r - r_pert;
    let dr_mag = dr.norm();
    let r_pert_mag = r_pert.norm();

    if dr_mag < 1e-30 || r_pert_mag < 1e-30 {
        return Vector3::zeros();
    }

    let dr3 = dr_mag.powi(3);
    let rp3 = r_pert_mag.powi(3);

    // Direct term + indirect (Sun recoil) term
    -gm_pert * (dr / dr3 + r_pert / rp3)
}

/// Apply a kick (velocity impulse) to all massive bodies due to mutual interactions.
/// Each body receives acceleration from all other bodies.
pub fn kick_bodies(bodies: &mut [MassiveBody], dt: f64) {
    let n = bodies.len();
    if n < 2 {
        return;
    }

    // Compute accelerations first, then apply (to avoid order dependence)
    let mut accels = vec![Vector3::<f64>::zeros(); n];

    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            accels[i] +=
                perturbation_acceleration(&bodies[i].state.pos, &bodies[j].state.pos, bodies[j].gm);
        }
    }

    for (body, accel) in bodies.iter_mut().zip(accels.iter()) {
        body.state.vel += dt * accel;
    }
}

/// Apply a kick to all test particles due to massive bodies.
/// Test particles do not affect each other or the massive bodies.
pub fn kick_particles(
    particles: &mut [StateVector],
    active: &[bool],
    bodies: &[MassiveBody],
    dt: f64,
) {
    for (i, particle) in particles.iter_mut().enumerate() {
        if !active[i] {
            continue;
        }

        let mut accel = Vector3::zeros();
        for body in bodies {
            accel += perturbation_acceleration(&particle.pos, &body.state.pos, body.gm);
        }

        particle.vel += dt * accel;
    }
}

/// Apply a kick to all test particles from massive bodies, using rayon for parallelism.
pub fn kick_particles_parallel(
    particles: &mut [StateVector],
    active: &[bool],
    bodies: &[MassiveBody],
    dt: f64,
) {
    use rayon::prelude::*;

    particles
        .par_iter_mut()
        .zip(active.par_iter())
        .for_each(|(particle, &is_active)| {
            if !is_active {
                return;
            }
            let mut accel = Vector3::zeros();
            for body in bodies {
                accel += perturbation_acceleration(&particle.pos, &body.state.pos, body.gm);
            }
            particle.vel += dt * accel;
        });
}

/// J2 oblateness acceleration.
/// Adds the acceleration due to the J2 (and optionally J4) oblateness of a body
/// centered at the origin with equatorial radius R.
///
/// The J2 potential is:
///   Φ_J2 = -GM/r * (R/r)^2 * J2 * P2(sin(φ))
/// where φ is the latitude.
pub fn j2_acceleration(
    r: &Vector3<f64>,
    gm: f64,
    radius: f64,
    j2: f64,
    j4: Option<f64>,
) -> Vector3<f64> {
    let r_mag = r.norm();
    if r_mag < 1e-30 {
        return Vector3::zeros();
    }

    let r2 = r_mag * r_mag;
    let r5 = r2 * r2 * r_mag;
    let r7 = r5 * r2;
    let rr = radius * radius;

    let z2 = r.z * r.z;
    let z2_r2 = z2 / r2;

    // J2 acceleration components
    let factor_j2 = 1.5 * gm * j2 * rr;
    let ax_j2 = -factor_j2 * r.x / r5 * (1.0 - 5.0 * z2_r2);
    let ay_j2 = -factor_j2 * r.y / r5 * (1.0 - 5.0 * z2_r2);
    let az_j2 = -factor_j2 * r.z / r5 * (3.0 - 5.0 * z2_r2);

    let mut accel = Vector3::new(ax_j2, ay_j2, az_j2);

    // J4 acceleration if provided
    if let Some(j4_val) = j4 {
        let r4 = rr * rr;
        let r9 = r7 * r2;
        let z4_r4 = z2_r2 * z2_r2;

        let factor_j4 = 0.625 * gm * j4_val * r4;
        let bracket = 3.0 - 42.0 * z2_r2 + 63.0 * z4_r4;
        let bracket_z = 15.0 - 70.0 * z2_r2 + 63.0 * z4_r4;

        accel.x += factor_j4 * r.x / r9 * bracket;
        accel.y += factor_j4 * r.y / r9 * bracket;
        accel.z += factor_j4 * r.z / r9 * bracket_z;
    }

    accel
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::GM_SUN;
    use approx::assert_relative_eq;

    #[test]
    fn test_perturbation_self_consistency() {
        // Two bodies on opposite sides of the Sun
        let r = Vector3::new(5.0, 0.0, 0.0);
        let r_pert = Vector3::new(-5.0, 0.0, 0.0);
        let gm = 1e-4;

        let accel = perturbation_acceleration(&r, &r_pert, gm);

        // Direct term: points from r toward r_pert (negative x)
        // Indirect term: points from origin toward -r_pert (positive x)
        // Direct dominates since |r-r_pert| = 10 < |r_pert| = 5... wait, opposite:
        // direct: -gm * (r-r_pert)/|r-r_pert|^3 = -gm * (10,0,0)/1000
        // indirect: -gm * r_pert/|r_pert|^3 = -gm * (-5,0,0)/125
        // So: ax = -gm*(10/1000 + (-5)/125) = -gm*(0.01 - 0.04) = -gm*(-0.03) = 0.03*gm
        assert!(accel.x > 0.0); // Net acceleration is toward the perturber... no, let's check sign
                                // Actually: 0.01 - 0.04 = -0.03, then -gm * (-0.03) = +0.03*gm
                                // Positive x means away from perturber... that's the indirect term dominating.
                                // This is correct — the indirect term accounts for the Sun being accelerated toward the perturber.
        assert_relative_eq!(accel.x, gm * 0.03, epsilon = 1e-10);
    }

    #[test]
    fn test_j2_equatorial_vs_polar() {
        // J2 acceleration should be different at pole vs equator
        let gm = GM_SUN;
        let r_eq = Vector3::new(1.0, 0.0, 0.0);
        let r_pole = Vector3::new(0.0, 0.0, 1.0);

        let a_eq = j2_acceleration(&r_eq, gm, 0.00465, 0.01, None);
        let a_pole = j2_acceleration(&r_pole, gm, 0.00465, 0.01, None);

        // At equator (z=0), J2 makes gravity stronger (pulls inward)
        // At pole (z=r), J2 makes gravity weaker (pushes outward)
        assert!(a_eq.x < 0.0); // Pulled inward at equator
        assert!(a_pole.z > 0.0); // Pushed outward at pole
    }
}
