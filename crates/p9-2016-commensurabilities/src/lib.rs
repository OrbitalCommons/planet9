//! Reproduction of de la Fuente Marcos, de la Fuente Marcos & Aarseth (2016),
//! "Dynamical grouping / commensurabilities of extreme trans-Neptunian objects
//! and a trans-Plutonian planet" (arXiv:1604.06241).
//!
//! HEADLINE. The extreme trans-Neptunian objects (ETNOs) are not a random
//! population: their orbital *angles* — the longitude of the ascending node Ω
//! and the longitude of perihelion ϖ = ω + Ω — are grouped on the circle, and
//! their orbital *periods* sit close to small-integer mean-motion
//! commensurabilities with one another and with a putative distant planet.
//! Both features are hard to reconcile with a population whose angles and
//! period ratios are drawn uniformly, and point to one (or more) distant
//! perturbing planet shepherding the ETNOs into near-resonant, apsidally and
//! nodally aligned configurations.
//!
//! WHAT THIS CRATE COMPUTES (real quantities, no hard-coded answers):
//!
//!   * `clustering` — the grouping significance of Ω and ϖ for the vetted
//!     `p9_core::data::etno` sample, via `p9_core::analysis::circular` (mean
//!     resultant length R̄ and the small-n-corrected Rayleigh p-value), shown
//!     against a seeded uniform-angle control.
//!
//!   * `commensurability` — a commensurability statistic measuring how close
//!     the ETNO semi-major axes are to small-integer period ratios, both among
//!     themselves (pairwise) and against a scanned candidate distant planet
//!     (reusing `p9_core::analysis::resonance::resonance_semi_major_axis`),
//!     shown to exceed a seeded uniform-control distribution via Monte Carlo.
//!
//! HONEST FINDING (the headline residual of this crate). The two parts of the
//! paper's argument behave very differently against a seeded random control:
//!
//!   * The **angle grouping is real and strong.** Both ϖ and Ω of the observed
//!     sample are far more grouped than uniform random angles (ϖ Rayleigh
//!     p ≈ 0.02, R̄ ≈ 0.45; the MC exceedance over a uniform control is a few
//!     percent). This is the robust, distribution-free signal.
//!
//!   * The **commensurability statistic does NOT survive a random control.**
//!     The observed ETNO period ratios sit close to small-integer ratios in
//!     absolute terms (mean pairwise distance ≈ 0.04 at max integer 9), but a
//!     uniform random population drawn over the same semi-major-axis range is
//!     *at least as commensurate* ~97% of the time, and a planet scan finds a
//!     comparably commensurate a₉ for random sets too (exceedance ≈ 0.53 at
//!     max integer 9, worse at low order). This is exactly the degeneracy the
//!     sibling crate `p9-2018-resonance` (Bailey, Brown & Batygin 2018)
//!     quantifies: the small-integer ratio grid is dense enough that "near a
//!     commensurability" is not rare. We report the exceedance honestly rather
//!     than manufacture significance — the commensurability argument is a
//!     suggestive arrangement, not a statistically significant excess. The
//!     residual is documented in REPRODUCTION_NOTES.md.
//!
//! Published reference values are kept as labelled constants in `published`.
//!
//! This crate reuses, and never duplicates, `p9_core`: the ETNO table
//! (`data::etno`), circular statistics (`analysis::circular`) and the
//! resonance-location arithmetic (`analysis::resonance`).

pub mod clustering;
pub mod commensurability;

/// Published reference values from de la Fuente Marcos, de la Fuente Marcos &
/// Aarseth (2016) and the closely related de la Fuente Marcos & de la Fuente
/// Marcos (2014, 2016) grouping papers. Kept as labelled constants for
/// cross-checks; the crate computes its own statistics from the data and pins
/// the residuals honestly.
pub mod published {
    /// The core dynamically grouped ETNO set the 2016 paper integrates
    /// (Sedna, 2012 VP113, 2004 VN112, 2007 TG422, 2013 RF98 + companions).
    /// The number of objects in the grouped sample analysed.
    pub const N_GROUPED_ETNOS: usize = 6;

    /// The paper frames the grouping as significant at well above the 2σ
    /// level; the related de la Fuente Marcos & de la Fuente Marcos (2014)
    /// node/argument-of-perihelion grouping is quoted at the ~0.0001 level.
    /// Used only as a labelled scale for "the observed angles cluster far more
    /// than a uniform population", not asserted as an exact reproduction.
    pub const NODE_GROUPING_P_REFERENCE: f64 = 1.0e-4;

    /// The trans-Plutonian planet the commensurabilities point to is placed
    /// well beyond the ETNOs, in the several-hundred-AU range; the companion
    /// resonance-prediction papers localize it near ~660 AU. Labelled scale
    /// only — the commensurability statistic here scans a range rather than
    /// asserting a single value.
    pub const CANDIDATE_PLANET_A_AU: f64 = 700.0;
}
