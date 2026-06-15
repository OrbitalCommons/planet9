//! Reproduction of Becker, Adams et al. (2017),
//! "Evaluating the dynamical stability of outer solar system objects in the
//! presence of Planet Nine" (arXiv:1706.06609).
//!
//! HEADLINE: distant TNOs are *not* locked in a single Planet Nine mean-motion
//! resonance for their lifetimes. They HOP between adjacent resonances
//! (metastable resonance sticking) while keeping approximate apsidal
//! anti-alignment, because neighbouring P9 resonances OVERLAP. Remaining in one
//! resonance is not a requirement for stability.
//!
//! This crate computes — with no hard-coded answers — the chain of exterior P9
//! n:1 (and n:2) mean-motion resonance semimajor axes near a representative
//! distant TNO, the libration width of each, and the Chirikov resonance-overlap
//! parameter between neighbours. It shows that beyond a computed semimajor axis
//! the resonances overlap (K ≳ 1) → a hopping band where an object cannot stay
//! in one resonance, while at smaller a they are isolated (K < 1).
//!
//! Reuses `p9_core::analysis::resonance::resonance_semi_major_axis` for the
//! resonance locations (the single workspace Kepler relation) and cross-checks
//! the overlap scaling against the validated
//! `p9_core::analysis::resonance::chirikov_overlap_parameter`. The distant-TNO
//! sample comes from `p9_core::data::etno::BROWN_2017_SAMPLE`.

pub mod chain;
pub mod tno;

pub use chain::{
    first_order_chain, hopping_threshold_au, mass_ratio, n_over_1_chain, n_over_2_chain,
    neighbour_overlap, overlap_profile, published, OverlapLink, P9Resonance,
};
pub use tno::{
    classify, classify_sample, hopping_threshold_au_nominal, nearest_n1_resonance, HoppingState,
};
