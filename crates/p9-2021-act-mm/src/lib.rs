//! Reproduction of Naess et al. (2021), "The Atacama Cosmology Telescope: A
//! search for Planet 9" (arXiv:2104.10264).
//!
//! A cold (~40 K) Planet Nine is hopeless in the near-IR (the WISE W1 thermal
//! flux sits far down the Wien tail; see `p9-2018-wise-search`) but is
//! comparatively bright at millimeter wavelengths, where ACT's 98/150/229 GHz
//! maps are sensitive. At 150 GHz a 40 K body is deep on the **Rayleigh–Jeans**
//! tail (hν/kT ≈ 0.18), so its flux scales as ν²·T rather than collapsing
//! exponentially. A blind ACT search over ~18,000 deg² found candidate
//! sources but no confirmed Planet Nine, setting a 95%-confidence point-source
//! limit of 4–12 mJy at 150 GHz.
//!
//! What this crate computes (no hard-coded answers):
//! - [`thermal_model`]: the mm-band blackbody flux density of a P9 of given
//!   effective temperature, radius (Neptune mass–radius) and distance, both
//!   exact-Planck and Rayleigh–Jeans, reusing the WISE crate's [`P9Thermal`].
//! - [`survey_model`]: the ACT search frequencies, the 4–12 mJy point-source
//!   limit, the ~18,000 deg² footprint and the 10-candidate count, as labelled
//!   reference constants.
//! - [`detectability`]: the maximum heliocentric distance at which a P9 reaches
//!   the ACT flux limit (closed-form 1/d² inversion, cross-checked by
//!   bisection), the mm-vs-WISE-W1 reach comparison, and an estimate of the
//!   expected number of spurious candidates over the footprint.
//!
//! [`P9Thermal`]: p9_2018_wise_search::thermal_model::P9Thermal

pub mod detectability;
pub mod survey_model;
pub mod thermal_model;
