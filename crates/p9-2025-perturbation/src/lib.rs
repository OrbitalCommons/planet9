//! Reproduction of Belyakov & Batygin (2025)
//! "The Perturbation Theory Approach to Stability in the Scattered Disk"
//!
//! Extends the perturbative scattered disk model to octupole order and beyond.
//! The 1:j and 3:j resonances from the octupole expansion do not individually
//! set new stability limits, but increasingly higher-index resonance chains
//! (2:j, 3:j, 4:j) explain chaotic behavior for Neptune-proximate orbits.
//!
//! Three stability boundary regimes:
//! - Fully chaotic layer: q = 13.37 + 0.55 * a^0.64
//! - 1:j resonance diffusion barrier: q = -63.1 + 19.4 * ln(a)
//! - Fine comb boundary (weak, *bounded* chaos): q = 10.28 + 6.01 * a^0.34
//!
//! Resonance widths use the single sourced p9-core Hansen/Chirikov
//! machinery; the 2:j chain is exactly consistent with
//! `p9_core::analysis::resonance::chirikov_overlap_parameter`.

pub mod chirikov;
pub mod hansen;
pub mod plots;
pub mod resonance;
pub mod spherical_harmonics;
pub mod stability_boundary;
