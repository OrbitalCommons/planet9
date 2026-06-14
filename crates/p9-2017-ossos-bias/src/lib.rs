//! Reproduction of Shankman et al. (2017)
//! "OSSOS VI. Striking Biases in the Detection of Large Semimajor Axis
//! Trans-Neptunian Objects" (arXiv:1706.05348).
//!
//! Headline result: wide-field surveys impose striking selection biases (in
//! discovery longitude, perihelion direction, ecliptic latitude, and
//! magnitude/H) such that a detected large-a TNO sample drawn from an
//! **intrinsically isotropic** parent population (uniform longitude of
//! perihelion ϖ) shows **apparent clustering** in ϖ. OSSOS's own large-a
//! detections are consistent with no intrinsic clustering once these biases
//! are modelled; the observed apsidal clustering can be explained by
//! detection bias.
//!
//! This crate reproduces the mechanism directly, with no hard-coded answer:
//!
//! 1. [`population`] generates a large seeded (`StdRng`) intrinsic population
//!    with ϖ drawn uniformly — verified to have intrinsic R̄ ≈ 0 and a
//!    non-significant Rayleigh test.
//! 2. [`detection`] applies a wide-field survey model: bright enough near
//!    perihelion (reflected-light r⁻⁴ law and a survey depth, both from
//!    `p9_core::analysis`), within a declination/ecliptic-latitude
//!    acceptance, and inside covered discovery-longitude windows — with the
//!    sky geometry from `p9_core::coords::sky`.
//! 3. [`bias`] runs the experiment and reports the mean resultant length R̄
//!    and Rayleigh p-value (`p9_core::analysis::circular`) of ϖ for the parent
//!    and the detected subsample.
//!
//! The reproduced quantity is the **bias-induced R̄** of the detected sample:
//! starting from an isotropic input, the survey model alone lifts the
//! detected ϖ mean resultant length from ~0 to ~0.3–0.6 with a highly
//! significant Rayleigh p — apparent clustering manufactured by bias.

pub mod bias;
pub mod detection;
pub mod population;

pub use bias::{run_default, run_experiment, BiasExperiment, ClusteringStats};
pub use detection::{detected_subsample, DetectionOutcome, SurveyModel};
pub use population::{generate_population, PopulationParams, SyntheticEtno};
