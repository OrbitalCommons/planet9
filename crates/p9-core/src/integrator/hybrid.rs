//! Mercury6-style hybrid integrator.
//!
//! Uses the Wisdom-Holman symplectic integrator for the majority of particles,
//! switching to Bulirsch-Stoer for any particle that enters a planet's Hill sphere.
//!
//! The switching criterion is:
//!   |r_particle - r_planet| < changeover * r_Hill
//!
//! where r_Hill = a * (m_planet / 3*M_sun)^(1/3) and changeover is typically 3.0.
//!
//! Reference: Chambers (1999), "A hybrid symplectic integrator..."

use crate::constants::GM_SUN;
use crate::integrator::bulirsch_stoer;
use crate::integrator::kepler_step::kepler_drift;
use crate::integrator::kick;
use crate::types::{MassiveBody, SimConfig, StateVector};

/// Determine which particles need Bulirsch-Stoer integration.
/// Returns a boolean mask (true = needs BS).
pub fn find_close_encounters(
    particles: &[StateVector],
    active: &[bool],
    bodies: &[MassiveBody],
    changeover_hill: f64,
) -> Vec<bool> {
    let mut needs_bs = vec![false; particles.len()];

    for (i, particle) in particles.iter().enumerate() {
        if !active[i] {
            continue;
        }
        for body in bodies {
            let dr = (particle.pos - body.state.pos).norm();
            let r_planet = body.state.pos.norm();
            let r_hill = r_planet * (body.mass / 3.0).cbrt();
            if dr < changeover_hill * r_hill {
                needs_bs[i] = true;
                break;
            }
        }
    }

    needs_bs
}

/// Perform one hybrid step: WHM for most particles, BS for close encounters.
pub fn hybrid_step(
    bodies: &mut [MassiveBody],
    particles: &mut [StateVector],
    active: &mut [bool],
    dt: f64,
    config: &SimConfig,
) {
    let half_dt = 0.5 * dt;

    // Check which particles have close encounters BEFORE the step
    let needs_bs = find_close_encounters(particles, active, bodies, config.hybrid_changeover_hill);

    // === KICK dt/2 for ALL particles ===
    kick::kick_bodies(bodies, half_dt);
    kick::kick_particles_parallel(particles, active, bodies, half_dt);

    // === DRIFT dt ===
    // Bodies always use Kepler drift
    for body in bodies.iter_mut() {
        let gm_total = GM_SUN + body.gm;
        body.state = kepler_drift(&body.state, dt, gm_total);
    }

    // Particles: Kepler for most, BS for close encounters
    for (i, particle) in particles.iter_mut().enumerate() {
        if !active[i] {
            continue;
        }

        if needs_bs[i] {
            // Bulirsch-Stoer for this particle
            let (new_state, _) = bulirsch_stoer::bs_step(particle, bodies, dt, config.bs_epsilon);
            *particle = new_state;
        } else {
            // Standard Kepler drift
            *particle = kepler_drift(particle, dt, GM_SUN);
        }
    }

    // === KICK dt/2 ===
    kick::kick_bodies(bodies, half_dt);
    kick::kick_particles_parallel(particles, active, bodies, half_dt);

    // === BOUNDARY CHECK ===
    for (i, particle) in particles.iter().enumerate() {
        if !active[i] {
            continue;
        }
        let r = particle.pos.norm();
        if r < config.removal_inner_au || r > config.removal_outer_au {
            active[i] = false;
        }
    }
}
