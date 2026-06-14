//! Reproduction of Holman et al. (2025)
//! "A search for Planet Nine and distant solar system objects in
//! Pan-STARRS1, Part 1" (arXiv:2506.02144).
//!
//! An INDEPENDENT shift-and-stack / moving-object search of the
//! Pan-STARRS1 (PS1) 3-pi survey. By co-adding the multi-epoch PS1
//! exposures along trial Planet Nine motion vectors ("shift-and-stack"),
//! the search reaches fainter than any single PS1 exposure (r ~ 21.5
//! single-epoch, gaining ~1 mag in the stack), over the 3-pi footprint
//! (declination > -30 deg). No Planet Nine is found, so the search sets
//! exclusion limits across the P9 mass-distance parameter space.
//!
//! ## Relationship to `p9-2024-panstarrs`
//!
//! This crate is INDEPENDENT OF AND COMPLEMENTARY TO the existing
//! `p9-2024-panstarrs` crate (Brown, Holman & Batygin 2024), which models
//! the catalog-level moving-object linking search of the same PS1 3-pi
//! data. Holman et al. (2025) instead perform a pixel-level
//! shift-and-stack search that recovers fainter (hence more distant /
//! lower-mass) objects below the catalog detection threshold. Both share
//! the same survey depth anchor (PS1 3pi r = 21.5) and the same
//! declination > -30 deg footprint from `p9-core`; this crate reuses that
//! shared survey/photometry logic rather than duplicating it.

pub mod exclusion;
pub mod photometry;
pub mod survey_model;

/// Published reference figures from Holman et al. (2025), kept as labelled
/// constants (no computed quantity is hard-coded to these — they are pinned
/// against in tests as sanity anchors).
pub mod published {
    /// PS1 3-pi single-epoch r-band 50%-completeness depth (Chambers et al.
    /// 2016; shared `p9-core` survey table entry "PS1 3pi").
    pub const PS1_R_DEPTH: f64 = 21.5;

    /// Effective depth gain from co-adding the multi-epoch PS1 stack along a
    /// trial motion vector. Holman et al. report reaching ~1 mag deeper than
    /// a single exposure in the shift-and-stack search.
    pub const SHIFT_STACK_DEPTH_GAIN_MAG: f64 = 1.0;

    /// Southern declination limit of the PS1 3-pi footprint (degrees).
    pub const PS1_DEC_LIMIT_DEG: f64 = -30.0;

    /// Order-of-magnitude scale of the parameter-space fraction this single
    /// independent PS1 search excludes for a Brown & Batygin reference P9
    /// population: a few tens of percent (the search is sensitive but does
    /// not cover the whole prior). Used only to bound test assertions.
    pub const EXCLUDED_FRACTION_ORDER: f64 = 0.5;
}
