//! Reproduction of Meisner, Bromley, Kenyon & Anderson (2016),
//! "Searching for Planet Nine with coadded WISE and NEOWISE data"
//! (arXiv:1611.00015).
//!
//! Meisner et al. build *inertial coadds* of WISE/NEOWISE W1 exposures binned
//! into ∼1-day intervals and shift-and-stack them along trial Planet Nine
//! motions. Coadding deepens the limiting magnitude by the √N stacking law,
//! extending the search from a single-exposure reach of ~430 AU to ~800 AU over
//! ∼2000 deg² (W1 < 16.66 at 90% completeness). The non-detection of a moving
//! source excludes a region of P9 (mass, distance) space. For a *cold* body the
//! paper notes W2 (4.6 µm) is more favorable than W1 (3.4 µm).
//!
//! This crate reuses the sibling `p9-2018-wise-search` aggressively — its
//! `P9Thermal` geometry/physics, its `WiseSurvey` footprint and logistic
//! efficiency, and its orbit→galactic-latitude sky mapping — and adds only the
//! genuinely new 2016 pieces:
//!
//! - [`coadd`]: the √N (1.25·log10 N) coadd depth gain over the single-frame
//!   WISE depth (the shared survey table's W1 = 16.5), with the published
//!   coadd W1 = 16.66 kept as a labelled reference.
//! - [`band`]: the same thermal+reflected physics evaluated in **both** W1 and
//!   W2, so the cold-body W2-over-W1 advantage is computed, not assumed.
//! - [`detectability`]: the maximum heliocentric distance a given-mass P9 is
//!   detectable in, against a coadd depth — increasing with mass, deeper for
//!   the coadd than single frames, and farther in W2 than W1 for a cold body.
//! - [`exclusion`]: the excluded fraction of a seeded reference P9 population —
//!   in (0,1), increasing with mass, larger for the coadd than single frames.
//!
//! **Honest residual (same finding as the sibling).** At W1 (3.4 µm) a cold
//! (~40 K) P9 lies so far down the Wien tail that its thermal flux is utterly
//! negligible and the W1 detectability is set by *reflected* sunlight. W2 is
//! less hostile (the Wien-tail factor improves by ~exp[(hc/kT)(1/λ₁−1/λ₂)]),
//! which is exactly the 2016 paper's point — but at 40 K even W2 thermal
//! emission remains weak, so both bands' cold-P9 reach is reflected-dominated
//! at these distances. The W2-over-W1 advantage this crate pins is real and in
//! the right direction; its absolute magnitude depends on the (uncertain) cold
//! effective temperature, kept as a documented model input.

pub mod band;
pub mod coadd;
pub mod detectability;
pub mod exclusion;
