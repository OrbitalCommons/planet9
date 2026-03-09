//! Reproduction of Brown (2017)
//! "Observational Bias and the Clustering of Distant Eccentric Kuiper Belt Objects"
//!
//! Implements bias-corrected statistical tests for orbital clustering of distant KBOs.
//! The core result: combined ϖ + ω clustering has only 0.025% probability of
//! arising by chance, even accounting for observational bias.
//!
//! TODO: Full bias function computation requires MPC discovery catalog
//! (survey pointings, depths, detection efficiencies). Currently uses
//! simplified perihelion-distance bias model from the paper.

pub mod bias_function;
pub mod clustering_test;
pub mod kbo_sample;
pub mod plots;
