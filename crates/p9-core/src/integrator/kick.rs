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

/// Direct gravitational acceleration only (no indirect term).
///
/// In the democratic heliocentric splitting (Duncan, Levison & Lee 1998) the
/// Sun-recoil physics is carried by the separate solar-drift substep, so the
/// kick must apply only the direct term.
#[inline]
pub fn direct_acceleration(r: &Vector3<f64>, r_pert: &Vector3<f64>, gm_pert: f64) -> Vector3<f64> {
    let dr = r - r_pert;
    let dr_mag = dr.norm();
    if dr_mag < 1e-30 {
        return Vector3::zeros();
    }
    -gm_pert * dr / dr_mag.powi(3)
}

/// Direct-only mutual kick for massive bodies (democratic heliocentric scheme).
pub fn kick_bodies_direct(bodies: &mut [MassiveBody], dt: f64) {
    let n = bodies.len();
    if n < 2 {
        return;
    }

    let mut accels = vec![Vector3::<f64>::zeros(); n];
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            accels[i] +=
                direct_acceleration(&bodies[i].state.pos, &bodies[j].state.pos, bodies[j].gm);
        }
    }

    for (body, accel) in bodies.iter_mut().zip(accels.iter()) {
        body.state.vel += dt * accel;
    }
}

/// Direct-only kick on test particles from massive bodies (democratic heliocentric).
pub fn kick_particles_direct(
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
            accel += direct_acceleration(&particle.pos, &body.state.pos, body.gm);
        }
        particle.vel += dt * accel;
    }
}

/// Direct-only particle kick using rayon for parallelism.
pub fn kick_particles_direct_parallel(
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
                accel += direct_acceleration(&particle.pos, &body.state.pos, body.gm);
            }
            particle.vel += dt * accel;
        });
}

/// Apply a kick (velocity impulse) to all massive bodies due to mutual interactions.
/// Each body receives acceleration from all other bodies.
///
/// Direct + indirect form, for heliocentric-velocity schemes. The
/// democratic-heliocentric integrator uses `kick_bodies_direct` instead.
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
        let z4_r4 = z2_r2 * z2_r2;

        let factor_j4 = 0.625 * gm * j4_val * r4;
        let bracket = 3.0 - 42.0 * z2_r2 + 63.0 * z4_r4;
        let bracket_z = 15.0 - 70.0 * z2_r2 + 63.0 * z4_r4;

        accel.x += factor_j4 * r.x / r7 * bracket;
        accel.y += factor_j4 * r.y / r7 * bracket;
        accel.z += factor_j4 * r.z / r7 * bracket_z;
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

    /// Perturbing potential U_n = (gm/r) J_n (R/r)^n P_n(sin φ) whose negative
    /// gradient `j2_acceleration` implements (verified for J2 against the
    /// closed form used in the function body).
    fn oblateness_potential(r: &Vector3<f64>, gm: f64, radius: f64, j2: f64, j4: f64) -> f64 {
        let r_mag = r.norm();
        let u = (r.z * r.z) / (r_mag * r_mag);
        let rr = radius * radius;
        let u2_term = gm * j2 * rr / (r_mag * r_mag * r_mag) * (3.0 * u - 1.0) / 2.0;
        let u4_term = gm * j4 * rr * rr / r_mag.powi(5) * (35.0 * u * u - 30.0 * u + 3.0) / 8.0;
        u2_term + u4_term
    }

    /// The acceleration must equal -∇U for the J2+J4 perturbing potential;
    /// checked by central differences over a grid of radii and latitudes.
    /// A wrong radial power in either term fails this immediately.
    #[test]
    fn test_j2_j4_acceleration_matches_potential_gradient() {
        let gm = GM_SUN;
        let radius = 0.0342; // effective JSU ring radius scale, AU
        let j2 = 1.0e-3;
        let j4 = -4.0e-5;
        let h = 1e-6;

        for &r_mag in &[2.0, 11.0, 36.0, 120.0] {
            for &lat_deg in &[0.0, 25.0, 60.0, 88.0_f64] {
                let lat = lat_deg.to_radians();
                let pos = Vector3::new(r_mag * lat.cos(), 0.0, r_mag * lat.sin());

                let accel = j2_acceleration(&pos, gm, radius, j2, Some(j4));

                for axis in 0..3 {
                    let mut hi = pos;
                    let mut lo = pos;
                    hi[axis] += h;
                    lo[axis] -= h;
                    let grad = (oblateness_potential(&hi, gm, radius, j2, j4)
                        - oblateness_potential(&lo, gm, radius, j2, j4))
                        / (2.0 * h);
                    assert_relative_eq!(accel[axis], -grad, max_relative = 1e-6, epsilon = 1e-24);
                }
            }
        }
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
