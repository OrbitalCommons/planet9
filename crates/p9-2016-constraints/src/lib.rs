//! Reproduction of Brown & Batygin (2016)
//! "Observational Constraints on the Orbit and Location of Planet Nine"
//!
//! This crate performs the systematic parameter survey to determine which
//! Planet Nine orbital parameters are consistent with the observed KBO clustering:
//! 1. Planar parameter grid: (a₉, e₉, m₉) survey (Section 2)
//! 2. 3D inclination exploration: (i₉, ω₉) survey (Section 3)
//! 3. Detection limit assessment: brightness and survey coverage (Section 4)

pub mod clustering_metric;
pub mod detection_limits;
pub mod inclination_survey;
pub mod parameter_grid;
pub mod plots;
