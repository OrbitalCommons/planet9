//! Reproduction of Batygin, Morbidelli, Brown & Nesvorny (2024)
//! "Generation of Low-Inclination, Neptune-Crossing Trans-Neptunian Objects by Planet Nine"
//!
//! Examines 17 well-characterized multi-opposition TNOs with a > 100 AU,
//! i < 40 deg, and q < 30 AU (Neptune-crossing). The P9-inclusive model
//! yields zeta = -7.9 (KS p = 0.41), consistent with observations, while the
//! P9-free null model gives zeta = -16.5 (KS p = 0.0034, i.e. ~2.7 sigma
//! one-sided). The paper's "~5 sigma" statement refers to zeta = -16.5
//! lying ~5 null standard deviations (sigma = 1.8) below the Monte-Carlo
//! null mean <zeta> = -7.2 — a different statistic from the KS p-value.

pub mod hypothesis_test;
pub mod observed_tnos;
pub mod plots;
pub mod simulation;
