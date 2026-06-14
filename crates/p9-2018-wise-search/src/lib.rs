//! Reproduction of Meisner, Bromley, Kenyon & Anderson (2018)
//! "Searching for Planet Nine with Coadded WISE and NEOWISE-Reactivation
//! Images" (arXiv:1712.04950).
//!
//! A wide-area WISE W1/W2 thermal-infrared shift-and-stack search for a slow
//! moving point source detects no Planet Nine. The W1 (3.4 µm) co-add depth
//! (W1 ≈ 16.5 Vega, `p9_core::analysis::surveys`) over a near-all-sky
//! footprint excludes a region of the Planet Nine mass–distance parameter
//! space — a thermal-IR analog of the optical ZTF/DES/PS1 exclusions.
//!
//! What this crate computes (no hard-coded answers):
//! - [`thermal_model`]: the W1 apparent magnitude of a P9 of given mass and
//!   heliocentric distance, from a documented thermal-emission + reflected
//!   sunlight model converted to the WISE W1 Vega magnitude system.
//! - [`survey_model`]: the W1 depth (from the shared survey table) and the
//!   near-all-sky footprint with a masked galactic-plane confusion band.
//! - [`detectability`]: the maximum heliocentric distance at which a
//!   given-mass P9 reaches the W1 depth (a finite, computed AU value that
//!   decreases with decreasing mass).
//! - [`exclusion`]: the fraction of a seeded reference P9 population WISE
//!   would have detected (in (0,1), increasing with assumed mass).
//!
//! Published reference numbers are kept as labelled constants
//! ([`survey_model::W1_DEPTH_REFERENCE`] = 16.5; the excluded distance range
//! of ~hundreds of AU is documented in the detectability tests).

pub mod detectability;
pub mod exclusion;
pub mod sky;
pub mod survey_model;
pub mod thermal_model;
