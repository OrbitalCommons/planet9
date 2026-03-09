//! Democratic heliocentric coordinate system for the Wisdom-Holman integrator.
//!
//! In this splitting, the Hamiltonian is decomposed as:
//!   H = H_Kepler + H_interaction
//!
//! H_Kepler: each body (planet or test particle) orbits the Sun in a Keplerian ellipse
//! H_interaction: planet-planet and planet-particle gravitational perturbations
//!
//! Coordinates are heliocentric positions with barycentric momenta (velocities).
//! This is the coordinate system used by mercury6 (Chambers 1999).

use nalgebra::Vector3;

use crate::types::{MassiveBody, StateVector};

/// Convert from heliocentric positions + velocities to democratic heliocentric.
/// In DH coords:
///   - Positions remain heliocentric
///   - Momenta (velocities) become barycentric:
///     v_DH_i = v_helio_i - v_bary
/// where v_bary = sum(m_i * v_i) / M_total (including the Sun at rest in heliocentric frame).
///
/// Returns the barycentric velocity of the system.
pub fn helio_to_democratic(
    bodies: &mut [MassiveBody],
    particles: &mut [StateVector],
    m_sun: f64,
) -> Vector3<f64> {
    // Compute total momentum in heliocentric frame
    let mut total_momentum = Vector3::zeros();
    let mut total_mass = m_sun; // Sun is at rest in heliocentric frame

    for b in bodies.iter() {
        total_momentum += b.mass * b.state.vel;
        total_mass += b.mass;
    }
    // Test particles are massless, don't contribute to momentum

    let v_bary = total_momentum / total_mass;

    // Transform velocities
    for b in bodies.iter_mut() {
        b.state.vel -= v_bary;
    }
    for p in particles.iter_mut() {
        p.vel -= v_bary;
    }

    v_bary
}

/// Convert from democratic heliocentric back to heliocentric.
pub fn democratic_to_helio(
    bodies: &mut [MassiveBody],
    particles: &mut [StateVector],
    v_bary: &Vector3<f64>,
) {
    for b in bodies.iter_mut() {
        b.state.vel += v_bary;
    }
    for p in particles.iter_mut() {
        p.vel += v_bary;
    }
}

/// Compute the barycentric correction velocity from current body states.
pub fn compute_v_bary(bodies: &[MassiveBody], m_sun: f64) -> Vector3<f64> {
    let mut total_momentum = Vector3::zeros();
    let mut total_mass = m_sun;

    for b in bodies {
        total_momentum += b.mass * b.state.vel;
        total_mass += b.mass;
    }

    total_momentum / total_mass
}
