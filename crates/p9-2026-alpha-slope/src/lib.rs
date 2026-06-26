//! Reproduction of Eldadi & Loeb (2026), "A Framework for Applying the
//! Loeb–Turner α-Slope Test to Archival Photometry of Trans-Neptunian Objects"
//! (arXiv:2605.17098).
//!
//! # The d⁻⁴ vs d⁻² physics
//!
//! Reflected sunlight reaches Earth after two inverse-square legs: the incident
//! solar flux at the body falls as d⁻² (Sun→body), and the reflected light then
//! falls as Δ⁻² on the way to the observer (body→Earth). Near opposition Δ ≈ d,
//! so the observed reflected flux scales as
//!
//!   F_refl ∝ d⁻⁴.
//!
//! Self-luminous emission — thermal radiation or any internal luminosity —
//! travels only the body→Earth leg, so for a fixed intrinsic brightness it
//! scales as
//!
//!   F_self ∝ d⁻².
//!
//! Fitting the power-law index α in F ∝ dᵅ across epochs (the slope of log₁₀F
//! against log₁₀d) is therefore a nature/technosignature diagnostic
//! (Loeb & Turner 2012): **α ≈ −4 ⇒ reflected sunlight (natural)**,
//! **α ≈ −2 ⇒ self-luminous (anomalous)**.
//!
//! # The Q1–Q6 funnel
//!
//! The paper builds a candidate bin for every (KBO × observatory × band)
//! combination in the MPC archive and passes each through a six-criterion
//! eligibility pipeline, then classifies the survivors by their fitted α:
//!
//! ```text
//!   8557  candidate bins        →  1089 pass Q1–Q3  →  186 pass Q4–Q6
//!                                                        ├── 53  reflected (α≈−4)
//!                                                        ├── 24  self-luminous (α≈−2)
//!                                                        └── 109 anomalous / other
//! ```
//!
//! Two diagnostics: Pluto (the brightest TNO) has 22 (observatory × band) bins
//! and recovers the reflected α ≈ −4 in *none* of them under a single
//! instrument+band — the archive can't cleanly run the test even on the easiest
//! target — and all 24 self-luminous detections come from Pan-STARRS PS1/PS2
//! alone, a calibration systematic rather than real self-luminosity.
//!
//! # What is reproduced
//!
//! - [`slope`]: the α-slope physics — a least-squares fit of α in F ∝ dᵅ
//!   ([`fit_alpha`]) and the Loeb–Turner classifier ([`classify_alpha`]).
//! - [`funnel`]: the published bin counts as `pub const` reference constants and
//!   a [`Funnel::is_consistent`] internal-consistency check, plus the Pluto and
//!   PS1/PS2 diagnostics.
//!
//! # p9-core reuse
//!
//! The synthetic photometry that exercises the regression is generated from
//! `p9_core::analysis::thermal::{reflected_flux_jy, thermal_flux_jy}` — the
//! canonical reflected (d⁻⁴) and thermal (d⁻²) flux laws — so the recovered
//! α ≈ −4 / α ≈ −2 are a verification of the shared core physics, not of a
//! local re-derivation. See `tests/headline.rs`.

pub mod funnel;
pub mod slope;

pub use funnel::{
    Funnel, ANOMALOUS, CANDIDATE_BINS, PASS_Q1_Q3, PASS_Q4_Q6, PLUTO_BINS,
    PLUTO_REFLECTED_RECOVERIES, PUBLISHED_FUNNEL, REFLECTED, SELF_LUMINOUS,
    SELF_LUMINOUS_ALL_PANSTARRS,
};
pub use slope::{
    classify_alpha, fit_alpha, magnitude_to_relative_flux, PhotometryPoint, SlopeClass, SlopeFit,
    ALPHA_REFLECTED, ALPHA_SELF_LUMINOUS, ALPHA_TOLERANCE,
};
