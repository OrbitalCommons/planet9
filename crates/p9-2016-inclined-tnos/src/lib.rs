//! Reproduction of Batygin & Brown (2016)
//! "Generation of Highly Inclined Trans-Neptunian Objects by Planet Nine"
//!
//! Simulates the two-step process that generates highly inclined TNOs:
//! 1. Planet Nine pumps inclination via Kozai-Lidov oscillations
//! 2. Neptune scattering reduces semi-major axis, decoupling from P9
//!
//! Uses full N-body integration with all four giant planets.

pub mod known_objects;
pub mod kozai_lidov;
pub mod plots;
pub mod simulation;
