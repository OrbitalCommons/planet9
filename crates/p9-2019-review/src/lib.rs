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
//! The critical semi-major axis a_c is computed with real secular machinery
//! (p9-core numerical Gauss-ring averaging + the giants' J2 field); the
//! paper's N-body-ensemble statistics (f_ϖ, η, μ, σ) are not emulated —
//! they require the full simulation suite and are omitted rather than
//! replaced by scaling guesses.

pub mod detection_prospects;
pub mod parameter_survey;
pub mod plots;
pub mod revised_parameters;
