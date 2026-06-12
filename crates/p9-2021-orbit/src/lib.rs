//! Reproduction of Brown & Batygin (2021)
//! "The Orbit of Planet Nine"
//!
//! The paper derives the Planet Nine posterior by MCMC over a likelihood
//! built from N-body simulations: key results m = 6.2 (+2.2/-1.3) Earth
//! masses, a = 380 (+140/-80) AU, i = 16 +/- 5 degrees, q = 300 (+85/-60)
//! AU, clustering significant at 99.6%.
//!
//! Two paths through this crate:
//!
//! 1. **From-scratch inference pipeline** (the paper's models 1-4):
//!    [`sim_grid`] integrates scattered-disk particles over a grid of
//!    (m9, a9, e9, i9) points, [`kde`] turns the surviving-particle
//!    distributions into conditional likelihoods of the observed ETNOs,
//!    [`gp`] emulates the likelihood surface with a Matern 3/2 Gaussian
//!    process, and [`mcmc`] (affine-invariant ensemble sampler) draws the
//!    posterior; [`pipeline`] wires them end to end. Default tests run the
//!    machinery at reduced scale; the published posterior itself is only
//!    reachable at the `#[ignore]`d paper scale (121 points x 16,800
//!    particles x 4 Gyr — estimated CPU-days, see REPRODUCTION_NOTES.md).
//!
//! 2. **Published-posterior summary** ([`posterior`]): split-normal
//!    resampling of the published medians/intervals (with a documented a-q
//!    correlation assumption). This is what the downstream survey crates
//!    (p9-2021-ztf, p9-2022-des, p9-2024-panstarrs) consume via
//!    [`reference_population`]; it represents the *published* posterior,
//!    not the reduced-scale pipeline output.
//!
//! [`statistical_measures`] recomputes the clustering significance from the
//! vetted ETNO table (0.989 computed vs 0.996 published, documented).

pub mod gp;
pub mod kde;
pub mod mcmc;
pub mod pipeline;
pub mod plots;
pub mod posterior;
pub mod reference_population;
pub mod sim_grid;
pub mod statistical_measures;
