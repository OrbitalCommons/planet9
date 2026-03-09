//! Reproduction of Bailey, Brown & Batygin (2018)
//! "Feasibility of a Resonance-Based Planet Nine Search"
//!
//! Characterizes mean-motion resonance (MMR) distributions between Planet Nine
//! and scattered Kuiper Belt objects. Shows that high-order resonances dominate,
//! and the resonance-based constraint on a₉ dissolves into a broad plateau
//! when the full resonance spectrum is considered.
//!
//! Key finding: P < 5% that all 6 observed KBOs reside in N/1 or N/2 resonances.
//!
//! TODO: Full 4 Gyr planar integration with J2 approximation.
//! Currently uses analytical resonance identification and probability analysis.

pub mod plots;
pub mod probability_analysis;
pub mod resonance_catalog;
pub mod simulation;
