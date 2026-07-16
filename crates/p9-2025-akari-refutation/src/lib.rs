//! Reproduction of Chen, Goto, Yamamura, Phan et al. (2025),
//! "A Far-Infrared Search for Planet Nine Using the AKARI All-Sky Survey"
//! (arXiv:2506.12854).
//!
//! This crate REFUTES the `p9-2025-iras-akari` reproduction of the
//! Phan et al. (2025) IRAS(1983)–AKARI(2006) candidate *pair*: the angular
//! motion implied by linking the two epochs is inconsistent with a single
//! bound distant-Solar-system object.
//!
//! ## The physical argument
//!
//! Chen et al. (§3.2) decompose the sky motion of a bound Planet Nine into
//! **parallax** (the geocentric reflex of Earth's orbit) and **proper
//! motion** (the object's own heliocentric orbital drift). Reproducing
//! their Cowan, Holder & Kaib (2016) numbers:
//!
//! - over **six months**, the parallactic swing is **10–25′** across the
//!   ~300–800 AU search range, while the proper motion is only **0.6–3.2′** —
//!   so on the half-year timescale *parallax dominates* (this is the basis of
//!   the Chen et al. AKARI-only search strategy).
//! - the parallax amplitude is `(1 AU / d)` rad; at **d = 600 AU** that is
//!   `1/600 rad ≈ 5.73′` (half-amplitude), `≈ 11.46′` peak-to-peak.
//!
//! The Phan et al. scheme instead links IRAS (mid-epoch 1983-06-24) to AKARI
//! (mid-epoch 2006-12-31), a **~23.5 yr baseline**. Crucially that baseline
//! is close to a *half-integer* number of years, so Earth returns to nearly
//! the *opposite* side of its orbit at the two mid-epochs: the parallax does
//! **not** cancel — but neither is it freely orientable. With Earth's
//! geometry fixed at the two real epochs, a bound object on a 500–700 AU
//! orbit produces a *specific, bounded* two-epoch separation. The Phan
//! candidate's 47.46′ separation requires an orbit and viewing geometry that
//! the bound-object model can only reach at one extreme of its envelope —
//! i.e. the "pair" is far more naturally explained as **two unrelated
//! infrared sources**. (A previous version also claimed the candidate's sky
//! DIRECTION mismatches the bound-object prediction; no position-angle
//! computation exists in this crate, so that claim was documentation only
//! and has been removed — computing the two-epoch parallax position angle
//! against the candidate's motion PA would be a genuine strengthening.)
//!
//! ## What is computed (no hard-coded answers)
//!
//! - [`expected_motion`]: the bound-object half-year parallax and proper
//!   motion as functions of distance, pinned to the Chen et al. 10–25′ /
//!   0.6–3.2′ ranges and the 5.73′ parallax half-amplitude at 600 AU.
//! - [`consistency`]: the bound-object two-epoch separation envelope at the
//!   *real* IRAS/AKARI epochs (parallax orientation fixed by Earth, only the
//!   object's orbit scanned), the Phan candidate's implied motion (reusing
//!   `p9-2025-iras-akari`'s encoded positions), and a quantitative
//!   inconsistency metric.
//!
//! Geometry and orbital helpers are reused from `p9-core` and
//! `p9-2025-iras-akari` rather than duplicated.

pub mod consistency;
pub mod expected_motion;

/// Published reference numbers from Chen et al. (2025), kept as labelled
/// constants so the reproductions below can be pinned against them.
pub mod published {
    /// Chen et al. §3.2: lower end of the bound-object six-month parallax
    /// range (arcmin), reached near the far end (~700–800 AU) of the search.
    pub const PARALLAX_6MO_MIN_ARCMIN: f64 = 10.0;
    /// Chen et al. §3.2: upper end of the six-month parallax range (arcmin),
    /// reached near the near end (~300 AU) of the search.
    pub const PARALLAX_6MO_MAX_ARCMIN: f64 = 25.0;
    /// Chen et al. §3.2: lower end of the bound-object six-month proper
    /// motion range (arcmin).
    pub const PROPER_MOTION_6MO_MIN_ARCMIN: f64 = 0.6;
    /// Chen et al. §3.2: upper end of the six-month proper motion range
    /// (arcmin).
    pub const PROPER_MOTION_6MO_MAX_ARCMIN: f64 = 3.2;
    /// Chen et al. distance range corresponding to the Phan scheme this work
    /// re-examines (AU).
    pub const DISTANCE_MIN_AU: f64 = 500.0;
    /// Upper bound of the Phan-scheme distance range (AU).
    pub const DISTANCE_MAX_AU: f64 = 700.0;
    /// Nominal IRAS–AKARI baseline quoted by Chen et al. ("separated by 23
    /// years"); the real mid-to-mid baseline is ~23.5 yr.
    pub const NOMINAL_BASELINE_YEARS: f64 = 23.0;
}
