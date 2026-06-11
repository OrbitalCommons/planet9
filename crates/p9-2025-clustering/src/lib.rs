//! Reproduction of Pichierri & Batygin (2025)
//! "Measuring the Degree of Clustering and Diffusion of Trans-Neptunian Objects"
//!
//! End-to-end pipeline (`pipeline::run_pipeline`): TNO clones sampled
//! from the vetted `p9_core::data::etno` table with correlated (a, e)
//! draws → p9-core WHM integration (Neptune direct, Jupiter/Saturn/Uranus
//! via the orbit-averaged J2 field, optional Planet Nine) → semi-major-
//! axis diffusion from an MSD(τ) fit → three-class stability
//! classification → clustering statistics (small-n-corrected Rayleigh
//! test plus a seeded bias-resampling Monte Carlo against a documented
//! longitude-dependent selection weighting) → EM fit of a two-component
//! von Mises mixture. All results are computed; the paper's reported
//! numbers are kept only as labelled reference values.
//!
//! Paper-reported values for comparison:
//! - Stable/metastable objects cluster at varpi* ~ 50 deg with sigma ~ 0.5 rad
//! - Unstable objects show bimodal distribution at ~25 deg and ~205 deg
//! - D_a,crit = 10^{-3} AU^2/yr for stability classification
//! - 25 clones per TNO, 4 Gyr integrations (reduced-scale smoke runs by
//!   default; the paper scale runs behind `#[ignore]`)

pub mod clone_generation;
pub mod clustering;
pub mod diffusion;
pub mod orbital_poles;
pub mod pipeline;
pub mod plots;
pub mod stability;
