//! Reproduction of Khain, Becker, Adams et al. (2020),
//! "Dynamical classification of distant TNOs in the presence of Planet Nine"
//! (resonance hopping II; arXiv:2010.02234).
//!
//! Headline: distant trans-Neptunian objects sort into dynamical CLASSES —
//! resonant (stuck in a single Planet Nine mean-motion resonance),
//! metastable / hopping (transiently sticking to a sequence of P9 resonances
//! and abruptly transitioning between commensurabilities), and
//! detached / scattering (not bound to any resonance). Resonance *sticking*
//! is most common near low-order n:1 resonances, where the resonances are
//! well separated; resonance *hopping* dominates where neighbouring P9
//! resonances overlap (Chirikov K ≳ 1), which happens at larger semimajor
//! axis because the n:1 resonance spacing collapses there.
//!
//! What this crate computes (no hard-coded answers): a seeded synthetic
//! distant-TNO population in semimajor axis, each object classified by
//!   (a) proximity to the nearest P9 n:1 / n:2 resonance, and
//!   (b) whether the neighbouring resonances OVERLAP (resonance-overlap K),
//! reusing `p9_core::analysis::resonance` and the resonance catalog /
//! `identify_resonance` machinery from `p9-2018-resonance`.
//!
//! It then reports the fraction of the population in each class and shows
//! that the hopping fraction RISES with semimajor axis (more overlap), while
//! clean resonant sticking dominates at smaller a.

pub mod classification;
pub mod population;
pub mod published;

pub use classification::{
    classify, resonance_overlap_k, Class, Classification, ResonanceLandscape,
};
pub use population::{class_fractions, synthetic_population, ClassFractions, Tno};
