//! Wisdom-Holman Mapping (WHM) symplectic integrator.
//!
//! Implements the second-order leapfrog scheme:
//!   1. Kick for dt/2 (interaction Hamiltonian)
//!   2. Drift for dt (Kepler step for each body)
//!   3. Kick for dt/2 (interaction Hamiltonian)
//!
//! This preserves a shadow Hamiltonian close to the true Hamiltonian,
//! giving bounded energy error over arbitrarily long integrations.
//!
//! Reference: Wisdom & Holman (1991), Chambers (1999)

use crate::constants::GM_SUN;
use crate::integrator::kepler_step::kepler_drift;
use crate::integrator::kick;
use crate::types::{MassiveBody, SimConfig, StateVector};

/// The WHM integrator state. Holds references to the system being integrated.
pub struct WhmIntegrator {
    /// Use parallel particle kicks
    pub parallel: bool,
}

impl WhmIntegrator {
    pub fn new() -> Self {
        Self { parallel: true }
    }

    /// Perform a single WHM step.
    ///
    /// `bodies`: massive bodies (planets, Planet Nine) in heliocentric coordinates
    /// `particles`: test particle state vectors
    /// `active`: which particles are still active
    /// `dt`: timestep in days
    /// `j2_bodies`: optional indices of bodies that contribute J2/J4 as a central force
    ///              (used when some planets are orbit-averaged rather than directly integrated)
    pub fn step(
        &self,
        bodies: &mut [MassiveBody],
        particles: &mut [StateVector],
        active: &mut [bool],
        dt: f64,
        config: &SimConfig,
    ) {
        let half_dt = 0.5 * dt;

        // === KICK dt/2 ===
        // Body-body interactions
        kick::kick_bodies(bodies, half_dt);

        // Body-particle interactions
        if self.parallel && particles.len() > 100 {
            kick::kick_particles_parallel(particles, active, bodies, half_dt);
        } else {
            kick::kick_particles(particles, active, bodies, half_dt);
        }

        // === DRIFT dt (Kepler step) ===
        // Each body drifts in its Kepler orbit around the Sun
        for body in bodies.iter_mut() {
            let gm_total = GM_SUN + body.gm;
            body.state = kepler_drift(&body.state, dt, gm_total);
        }

        // Each particle drifts in its Kepler orbit around the Sun
        for (i, particle) in particles.iter_mut().enumerate() {
            if !active[i] {
                continue;
            }
            *particle = kepler_drift(particle, dt, GM_SUN);
        }

        // === KICK dt/2 ===
        kick::kick_bodies(bodies, half_dt);

        if self.parallel && particles.len() > 100 {
            kick::kick_particles_parallel(particles, active, bodies, half_dt);
        } else {
            kick::kick_particles(particles, active, bodies, half_dt);
        }

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
}

/// Compute total energy of the system (for diagnostics).
/// Returns (kinetic, potential, total).
pub fn system_energy(
    bodies: &[MassiveBody],
    particles: &[StateVector],
    active: &[bool],
) -> (f64, f64, f64) {
    let mut kinetic = 0.0;
    let mut potential = 0.0;

    // Body kinetic energy (relative to Sun)
    for body in bodies {
        kinetic += 0.5 * body.mass * body.state.vel.norm_squared();
    }

    // Body-Sun potential
    for body in bodies {
        let r = body.state.pos.norm();
        if r > 0.0 {
            potential -= GM_SUN * body.mass / r;
        }
    }

    // Body-body potential
    for i in 0..bodies.len() {
        for j in (i + 1)..bodies.len() {
            let dr = (bodies[i].state.pos - bodies[j].state.pos).norm();
            if dr > 0.0 {
                potential -= bodies[i].gm * bodies[j].mass / dr;
            }
        }
    }

    // Particle kinetic energy (massless, but track for diagnostics using unit mass)
    for (i, p) in particles.iter().enumerate() {
        if !active[i] {
            continue;
        }
        kinetic += 0.5 * p.vel.norm_squared(); // unit mass
        let r = p.pos.norm();
        if r > 0.0 {
            potential -= GM_SUN / r;
            for body in bodies {
                let dr = (p.pos - body.state.pos).norm();
                if dr > 0.0 {
                    potential -= body.gm / dr;
                }
            }
        }
    }

    (kinetic, potential, kinetic + potential)
}

/// Compute the Jacobi constant for the circular restricted three-body problem.
/// Useful for testing. `mu` is the mass ratio m2/(m1+m2).
pub fn jacobi_constant(state: &StateVector, mu: f64, omega: f64) -> f64 {
    let x = state.pos.x;
    let y = state.pos.y;

    let r1 = ((x + mu).powi(2) + y.powi(2) + state.pos.z.powi(2)).sqrt();
    let r2 = ((x - 1.0 + mu).powi(2) + y.powi(2) + state.pos.z.powi(2)).sqrt();

    let v2 = state.vel.norm_squared();

    // Jacobi constant: C_J = -2E - 2 * omega * L_z in rotating frame
    // Or equivalently: C_J = (x^2 + y^2)*omega^2 + 2*(1-mu)/r1 + 2*mu/r2 - v^2
    omega * omega * (x * x + y * y) + 2.0 * (1.0 - mu) / r1 + 2.0 * mu / r2 - v2
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    #[test]
    fn test_whm_two_body_energy_conservation() {
        // Test planet on circular orbit at 5 AU (Jupiter-like)
        let a = 5.2;
        let gm_planet = 1e-3 * GM_SUN; // ~Jupiter mass
        let v = (GM_SUN / a).sqrt();

        let mut bodies = vec![MassiveBody {
            name: "Jupiter".to_string(),
            gm: gm_planet,
            mass: gm_planet / GM_SUN,
            state: StateVector::new(Vector3::new(a, 0.0, 0.0), Vector3::new(0.0, v, 0.0)),
            radius_au: 4.78e-4,
            j2: None,
            j4: None,
        }];

        let mut particles = vec![];
        let mut active = vec![];

        let config = SimConfig {
            dt: 100.0,
            t_start: 0.0,
            t_end: 1e6,
            removal_inner_au: 0.1,
            removal_outer_au: 1e6,
            snapshot_interval_days: 1e6,
            hybrid_changeover_hill: 3.0,
            bs_epsilon: 1e-12,
        };

        let integrator = WhmIntegrator::new();
        let (_, _, e0) = system_energy(&bodies, &particles, &active);

        // Integrate for 1000 steps
        for _ in 0..1000 {
            integrator.step(&mut bodies, &mut particles, &mut active, config.dt, &config);
        }

        let (_, _, e1) = system_energy(&bodies, &particles, &active);

        // Symplectic integrator should have bounded energy error
        let de = ((e1 - e0) / e0).abs();
        assert!(de < 1e-6, "Energy drift too large: {:.2e}", de);
    }

    #[test]
    fn test_whm_particle_on_circular_orbit() {
        // No planets, just a particle on a circular orbit
        let a = 30.0; // Neptune distance
        let v = (GM_SUN / a).sqrt();

        let mut bodies = vec![];
        let mut particles = vec![StateVector::new(
            Vector3::new(a, 0.0, 0.0),
            Vector3::new(0.0, v, 0.0),
        )];
        let mut active = vec![true];

        let config = SimConfig {
            dt: 300.0,
            t_start: 0.0,
            t_end: 1e8,
            removal_inner_au: 0.1,
            removal_outer_au: 1e6,
            snapshot_interval_days: 1e8,
            hybrid_changeover_hill: 3.0,
            bs_epsilon: 1e-12,
        };

        let integrator = WhmIntegrator::new();

        // One orbital period of Neptune
        let period = 2.0 * std::f64::consts::PI * (a * a * a / GM_SUN).sqrt();
        let n_steps = (period / config.dt).round() as usize;

        for _ in 0..n_steps {
            integrator.step(&mut bodies, &mut particles, &mut active, config.dt, &config);
        }

        // Should return close to starting position
        assert_relative_eq!(particles[0].pos.x, a, epsilon = 0.1);
        assert_relative_eq!(particles[0].pos.y, 0.0, epsilon = 0.1);
    }
}
