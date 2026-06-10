//! Reproduction of Pichierri & Batygin (2025)
//! "Measuring the Degree of Clustering and Diffusion of Trans-Neptunian Objects"
//!
//! Numerically integrates TNO clones and measures orbital diffusion to classify
//! objects as stable, metastable, or unstable. Measures clustering of longitude
//! of perihelion and orbital poles as functions of distance and stability class.
//!
//! Key results:
//! - Stable/metastable objects cluster at varpi* ~ 50 deg with sigma ~ 0.5 rad
//! - Unstable objects show bimodal distribution at ~25 deg and ~205 deg
//! - D_a,crit = 10^{-3} AU^2/yr for stability classification
//! - 25 clones per TNO, 4 Gyr integrations

pub mod clone_generation;
pub mod clustering;
pub mod diffusion;
pub mod orbital_poles;
pub mod plots;
pub mod stability;
