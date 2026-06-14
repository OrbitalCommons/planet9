//! Two 2025 distant-TNO discoveries that STRESS the Planet Nine clustering
//! hypothesis.
//!
//! Both papers report extreme trans-Neptunian objects (ETNOs) whose
//! longitude of perihelion ϖ = ω + Ω lies *outside* the clustered group that
//! motivates Planet Nine. Each is an object Planet Nine should have herded
//! into apsidal alignment, yet it is not aligned — so adding it to the
//! clustered sample *weakens* the clustering significance.
//!
//! 1. Cheng, Yang, et al. 2025, "2017 OF201" (arXiv:2505.15806): a very
//!    distant (a ≈ 830 AU, q ≈ 45 AU, e ≈ 0.95) dwarf-planet-sized ETNO whose
//!    ϖ is far from the clustered ETNO mean.
//! 2. Chen, Kavelaars, et al. 2025, "Ammonite / 2023 KQ14"
//!    (arXiv:2508.02162): a fourth "sednoid" (q ≈ 66 AU, a ≈ 252 AU, i ≈ 11°)
//!    whose ϖ is misaligned with the other three sednoids (Sedna, 2012 VP113,
//!    Leleakuhonua), reducing the apparent sednoid apsidal alignment.
//!
//! Everything here is computed from real orbital elements reusing
//! `p9_core::analysis::circular` and the `p9_core::data::etno` machinery; the
//! published headline numbers (830 AU; "fourth misaligned sednoid") are kept
//! only as labelled reference constants. See [`discoveries`] for element
//! provenance and [`stress`] for the clustering-degradation analysis.

pub mod discoveries;
pub mod stress;

pub use discoveries::{ammonite, of201, NewDiscovery};
pub use stress::{clustering_stat, ClusteringStat, StressResult};
