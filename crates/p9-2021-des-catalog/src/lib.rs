//! Reproduction of Bernardinelli et al. (2021),
//! "A search for trans-Neptunian objects with the six-year Dark Energy Survey"
//! (arXiv:2109.03758).
//!
//! Headline: the full six-year DES yields a catalog of **815 TNOs** (461 newly
//! reported, plus a Centaur and an Oort-cloud comet) over the ~5000 deg²
//! southern footprint, complete to **r ≈ 23.8** mag. The catalog contains
//! **16 "extreme" TNOs** (a > 150 AU, q > 30 AU). There is **no Planet Nine**
//! detection; the extreme-TNO subset is consistent with *azimuthal isotropy*
//! once the DES selection function is accounted for, and the well-characterized
//! detection efficiency over the footprint feeds the clustering / exclusion
//! analyses.
//!
//! This crate computes — from real, seeded calculations, never hard-coded
//! answers:
//!
//!   * the orbital-angle **clustering statistics** of the DES extreme-TNO
//!     subset (mean resultant length R̄ and Rayleigh p, via
//!     `p9_core::analysis::circular`) for ϖ, ω and Ω, both against a flat
//!     isotropic null and against the **DES selection function** — showing the
//!     DES sample alone is consistent with its selection (`clustering`);
//!
//!   * the catalog's **distant-object detection completeness** over the
//!     footprint as a function of magnitude / heliocentric distance, anchored
//!     on the shared `p9_core::analysis::surveys` DES depth (r = 23.8) and the
//!     declared 5000 deg² footprint with δ and galactic cuts
//!     (`completeness`); and a Planet-Nine **footprint-coverage fraction** in
//!     (0, 1).
//!
//! The DES extreme-TNO table (`catalog`) is transcribed from the survey's
//! published elements (the same DES eTNOs reproduced in `p9-2020-des-isotropy`,
//! reused here as the a > 150 / q > 30 subset). We reuse the shared circular
//! statistics and survey-depth tables from `p9-core` rather than reimplementing
//! them.
//!
//! Reference constants (catalog size 815, 16 extreme TNOs, footprint, depth)
//! are carried as labelled `reference` items and asserted as identities;
//! everything else is computed and pinned with documented residuals.

pub mod catalog;
pub mod clustering;
pub mod completeness;

/// Published reference numbers from Bernardinelli et al. (2021), carried as
/// labelled constants (not used to back-fill any computed quantity).
pub mod reference {
    /// Total TNOs in the six-year DES catalog (Table / abstract).
    pub const N_TNOS: usize = 815;
    /// Newly reported TNOs.
    pub const N_NEW_TNOS: usize = 461;
    /// Extreme TNOs (a > 150 AU, q > 30 AU) in the catalog.
    pub const N_EXTREME_TNOS: usize = 16;
    /// Effective wide-survey footprint area (deg²).
    pub const FOOTPRINT_DEG2: f64 = 5000.0;
    /// Catalog completeness depth in r (mag).
    pub const DEPTH_R: f64 = 23.8;
    /// One-line conclusion.
    pub const PAPER_CONCLUSION: &str =
        "The 815-object six-year DES TNO catalog (complete to r ~ 23.8 over \
         ~5000 deg^2) contains 16 extreme TNOs whose orbital angles are \
         consistent with azimuthal isotropy under the DES selection function; \
         no Planet Nine is detected.";
}
