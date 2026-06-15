//! Reproduction of Siraj, Chyba & Tremaine (2026), "Measuring Apsidal
//! Clustering" (arXiv:2604.25990).
//!
//! The paper introduces a conditional-likelihood estimator for apsidal
//! (longitude-of-perihelion ϖ) clustering of the extreme trans-Neptunian
//! objects (ETNOs). The estimator is a likelihood-ratio test of an apsidally
//! clustered model against an isotropic null:
//!
//! - **H1** — ϖ is drawn from a von Mises distribution
//!     f(ϖ | µ, κ) = exp(κ·cos(ϖ − µ)) / (2π·I₀(κ))
//!   with (µ, κ) fitted by maximum likelihood from the sample
//!   (µ = circular mean, κ = A⁻¹(R̄)).
//! - **H0** — ϖ is uniform on the circle, density 1/2π (κ = 0).
//!
//! The log-likelihood ratio Λ = ln L(H1) − ln L(H0) is converted to a
//! Gaussian-equivalent significance via Wilks' theorem (see [`sigma`]).
//!
//! Headline result reproduced: the stable-ETNO sample of n = 21 objects gives
//! an apsidal-clustering significance of ≈2.7σ. Expanding the sample to n = 25
//! by adding 4 less-aligned stable objects DROPS the significance to ≈1.9σ —
//! the SAME estimator on a larger, diluted sample yields lower significance.
//!
//! Circular summary statistics (circular mean, mean resultant length) and the
//! von Mises helpers (`bessel_i0`, `kappa_from_r_bar`) come from
//! `p9_core::analysis::circular`.

pub mod estimator;
pub mod samples;
pub mod significance;

pub use estimator::{log_likelihood_ratio, von_mises_fit, ApsidalFit};
pub use samples::{stable_21, stable_25};
pub use significance::{p_value_to_sigma, sigma};

/// Published apsidal-clustering significance for the stable n = 21 sample
/// (Siraj, Chyba & Tremaine 2026). Reference constant only.
pub const PUBLISHED_SIGMA_21: f64 = 2.7;

/// Published apsidal-clustering significance for the diluted n = 25 sample
/// (Siraj, Chyba & Tremaine 2026). Reference constant only.
pub const PUBLISHED_SIGMA_25: f64 = 1.9;
