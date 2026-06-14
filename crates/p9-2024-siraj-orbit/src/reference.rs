//! Published reference values from Siraj, Chyba & Tremaine (2024),
//! "Orbit of a Possible Planet X" (arXiv:2410.18170).
//!
//! These are *labelled reference constants only* — they are not used to
//! short-circuit any computation. The inferred MAP is recovered numerically by
//! [`crate::posterior::ForcingPosterior::map_estimate`]; these values are the
//! published numbers the regression tests pin the computation against, and the
//! single scale-setting calibration in [`crate::posterior::calibrated_posterior`].
//!
//! The headline of the paper is that this independent inference yields a
//! *markedly different* perturber than Brown & Batygin (2021): a lower mass, a
//! smaller semi-major axis, and a lower inclination. The two solutions are
//! incompatible (see [`crate::tension`]).

/// Inferred perturber mass (Earth masses).
pub const SIRAJ_2024_MASS_EARTH: f64 = 4.4;

/// Inferred perturber semi-major axis (AU).
pub const SIRAJ_2024_A_AU: f64 = 290.0;

/// Inferred perturber inclination (degrees).
pub const SIRAJ_2024_I_DEG: f64 = 6.8;

/// Approximate 1σ uncertainty on the inferred mass (Earth masses). Used as the
/// Siraj-side scatter when computing the tension against Brown & Batygin 2021.
pub const SIRAJ_2024_MASS_SIGMA_EARTH: f64 = 1.0;

/// Approximate 1σ uncertainty on the inferred semi-major axis (AU).
pub const SIRAJ_2024_A_SIGMA_AU: f64 = 30.0;

/// Approximate 1σ uncertainty on the inferred inclination (degrees).
pub const SIRAJ_2024_I_SIGMA_DEG: f64 = 1.7;
