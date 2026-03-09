//! Reproduction of Belyakov, Bernardinelli & Brown (2022)
//! "Limits on the Detection of Planet Nine in the Dark Energy Survey"
//!
//! Evaluates the Dark Energy Survey's sensitivity to Planet Nine by injecting
//! synthetic orbits into the DES cadence and measuring recovery efficiency.
//! Of 11709 synthetic orbits that cross the DES footprint, 10187 (87%) are
//! successfully recovered. When combined with ZTF, DES contributes an
//! additional 5% unique exclusion for a combined 61.2% exclusion of the
//! prior parameter space.

pub mod color_models;
pub mod plots;
pub mod recovery_analysis;
pub mod survey_model;
