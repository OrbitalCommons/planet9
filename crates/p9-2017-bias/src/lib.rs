//! Reproduction of Brown (2017)
//! "Observational Bias and the Clustering of Distant Eccentric Kuiper Belt Objects"
//!
//! Implements bias-corrected statistical tests for orbital clustering of distant KBOs.
//! The core result: combined ϖ + ω clustering has only 0.025% probability of
//! arising by chance, even accounting for observational bias.
//!
//! The bias model is a documented simplification of the paper's MPC
//! pointing/depth catalog: detectability from H magnitudes and the
//! reflected-light r⁻⁴ law, longitudinal galactic-plane suppression at the
//! ecliptic crossings (λ ≈ 95°/275°), and an ecliptic-latitude coverage
//! penalty. See `bias_function::BiasParams` for the assumptions.

pub mod bias_function;
pub mod clustering_test;
pub mod kbo_sample;
pub mod plots;
