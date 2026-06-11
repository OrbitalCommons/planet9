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
//! The planar integration applies the giant-planet J2 quadrupole
//! (`p9_core::forces::ExtraForce::J2Jsu`) in the WHM kick and classifies
//! resonance membership by libration of recorded resonant-angle time series
//! (`p9_core::analysis::resonance`); the 4 Gyr paper-scale sweep sits
//! behind an `#[ignore]`d test.

pub mod plots;
pub mod probability_analysis;
pub mod resonance_catalog;
pub mod simulation;
