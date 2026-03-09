//! Reproduction of Batygin, Adams, Brown & Becker (2019)
//! "The Planet Nine Hypothesis"
//!
//! Comprehensive review consolidating evidence from all prior papers.
//! Key new contribution: revised parameter estimates from a 1,794-simulation
//! ensemble (1,134 semi-averaged + 660 fully resolved).
//!
//! Revised parameters: m₉ ~ 5-10 M_Earth, a₉ ~ 400-800 AU,
//! e₉ ~ 0.2-0.5, i₉ ~ 15-25°
//!
//! TODO: Full semi-averaged integrator (J2+Neptune+P9) and N-body parameter grid.
//! Currently implements the parameter survey framework and detection prospects.

pub mod detection_prospects;
pub mod parameter_survey;
pub mod plots;
pub mod revised_parameters;
