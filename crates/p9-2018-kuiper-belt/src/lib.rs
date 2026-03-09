//! Reproduction of Khain, Batygin & Brown (2018)
//! "The Generation of the Distant Kuiper Belt by Planet Nine
//! from an Initially Broad Perihelion Distribution"
//!
//! Compares narrow vs broad initial perihelion distributions to show that
//! only a broad distribution (q up to ~300 AU) produces both aligned and
//! anti-aligned populations — the bimodal structure observed in the real
//! distant Kuiper Belt.
//!
//! Planet Nine parameters: a₉ = 700 AU, e₉ = 0.6, i₉ = 0°, m₉ = 10 M_Earth.
//!
//! TODO: Full 4 Gyr integration with Mercury6-style hybrid integrator.
//! Currently uses WHM with the same integration framework as other subcrates.

pub mod plots;
pub mod population_analysis;
pub mod simulation;
