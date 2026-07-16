//! Reproduction of Hadden, Li, Payne & Holman (2018)
//! "Chaotic dynamics of trans-Neptunian objects perturbed by Planet Nine"
//! (arXiv:1712.06547).
//!
//! HEADLINE: a distant Planet Nine renders large-`a` ETNO orbits *chaotic*
//! through the *overlap* of neighbouring mean-motion / secular resonances. The
//! chaos sets in above a critical eccentricity (equivalently, above a critical
//! semi-major axis for fixed perihelion); regular, non-chaotic orbits survive
//! at low eccentricity. The web of overlapped j:1 resonances widens with
//! Planet Nine's mass.
//!
//! This crate computes that picture from first principles using the shared
//! `p9-core` machinery — never hard-coding the answers:
//!
//! - [`chaos`] builds a resonance-overlap (Chirikov) chaos indicator over the
//!   (a, e) plane. The overlap parameter K = (sum of neighbouring resonance
//!   half-widths) / (spacing) crosses unity along a computed boundary in e at
//!   each a. Below it (low e) the resonances are isolated and the motion is
//!   regular; above it they overlap and the motion is chaotic.
//! - [`width`] computes the secular/resonant libration width and shows it
//!   grows monotonically with Planet Nine's mass.
//! - The resonance locations are cross-checked against the analytic Kepler
//!   relation a_j = a9 j^{-2/3} to 1e-9.
//!
//! REUSE: the resonance locations come from
//! `p9_core::analysis::resonance::resonance_semi_major_axis`, the overlap
//! zone from `p9_core::analysis::resonance::overlap_zone_width_fraction`,
//! and the width machinery from `p9_core::analysis::hansen`. Nothing here
//! re-derives those.
//!
//! Published claims are kept as labelled reference constants in [`reference`];
//! residuals against them are documented in the tests.

pub mod chaos;
pub mod reference;
pub mod width;
