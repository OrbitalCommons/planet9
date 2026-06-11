//! Reproduction of Brown & Batygin (2019)
//! "Orbital Clustering in the Distant Solar System"
//!
//! Extends the 2017 bias analysis to 14 KBOs with a > 230 AU, analyzing
//! clustering in both longitude of perihelion AND orbital pole position
//! simultaneously using Poincaré action-angle variables.
//!
//! Key result: combined clustering probability = 0.2% (99.8% confidence).
//! OSSOS survey is insensitive to the clustering signal.
//!
//! The Monte Carlo null is survey-bias weighted (galactic-plane suppression
//! of perihelion longitudes + declination coverage), a simplified
//! single-factor stand-in for the paper's full MPC-metadata bias maps; the
//! uniform null is kept as a labeled baseline.

pub mod clustering_analysis;
pub mod kbo_sample;
pub mod ossos_comparison;
pub mod plots;
pub mod poincare_variables;
