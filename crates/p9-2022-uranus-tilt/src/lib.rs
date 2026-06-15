//! Reproduction of Lu & Laughlin (2022), arXiv:2207.11823
//! "The Primordial Spin-Orbit Resonance Origin of Uranus' Obliquity"
//!
//! The paper proposes that Uranus' present-day ~98° obliquity was excited not
//! by a giant impact but by capture into a *secular spin–orbit resonance*
//! driven by the outward migration of a distant Planet Nine. As Planet Nine
//! migrates outward, the forced nodal-precession frequency of Uranus' orbit
//! sweeps slowly through the planet's spin-axis precession frequency `α`. When
//! the commensurability `α cos θ ≈ |g|` is crossed adiabatically the spin axis
//! is captured into a Cassini state whose obliquity grows as the resonance
//! migrates, reaching ~90–98°.
//!
//! Three pieces are reproduced from real quantities (no hard-coded answers):
//!
//!   * [`spin_axis`] — the spin-axis precession constant
//!     α = (3 n² / 2 ω)·(J₂ + q)/(λ + l) for Uranus, including the satellite
//!     contributions q, l. Pinned against the paper's α ≈ 0.045 arcsec/yr.
//!   * [`cassini`] — the Colombo-top Cassini-state geometry: the four Cassini
//!     states, the critical ratio (α/g)_crit where states 1 and 4 merge, and
//!     the resonant obliquity at given α/g.
//!   * [`resonance_capture`] — adiabatic sweep of the orbital precession
//!     frequency through resonance, driving the obliquity to ~98°, with the
//!     monotone dependence on the controlling parameter α.

pub mod cassini;
pub mod resonance_capture;
pub mod spin_axis;

use p9_core::constants::DEG2RAD;

/// Uranus' observed obliquity (deg). IAU value; the paper's target is 98°.
pub const URANUS_OBLIQUITY_DEG: f64 = 97.77;

/// Uranus' observed obliquity (rad).
pub fn uranus_obliquity_rad() -> f64 {
    URANUS_OBLIQUITY_DEG * DEG2RAD
}

/// Lu & Laughlin (2022) headline spin-axis precession constant for Uranus
/// today (arcsec/yr). Reference label only; [`spin_axis::alpha_arcsec_per_yr`]
/// computes it.
pub const ALPHA_PRESENT_ARCSEC_PER_YR: f64 = 0.045;
