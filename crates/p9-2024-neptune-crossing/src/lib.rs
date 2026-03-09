//! Reproduction of Batygin, Morbidelli, Brown & Nesvorny (2024)
//! "Generation of Low-Inclination, Neptune-Crossing Trans-Neptunian Objects by Planet Nine"
//!
//! Examines 17 well-characterized multi-opposition TNOs with a > 100 AU, i < 40 deg,
//! and q < 30 AU (Neptune-crossing). The P9-inclusive model yields zeta = -7.9 (p = 0.41),
//! consistent with observations, while the P9-free null model gives zeta = -16.5
//! (p = 0.0034), rejected at approximately 5 sigma significance.

pub mod hypothesis_test;
pub mod observed_tnos;
pub mod plots;
pub mod simulation;
