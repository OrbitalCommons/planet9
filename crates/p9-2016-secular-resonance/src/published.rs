//! Labelled reference constants from Beust (2016), arXiv:1605.02473.
//!
//! These are the paper's stated values, kept as named constants so the
//! computed quantities can be compared against them in tests without any of
//! them being used as an input or hard-coded answer.

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;

/// Nominal Planet Nine used in Beust's secular study, following Batygin &
/// Brown (2016): m₉ = 10 M⊕, a₉ = 700 AU, e₉ = 0.6. (Beust scans a range of
/// configurations; this is the fiducial.)
pub fn nominal_p9() -> P9Params {
    P9Params::nominal_2016()
}

/// The libration islands the paper reports are **apsidally anti-aligned**:
/// the center of the libration of Δϖ = ϖ − ϖ₉ sits at π (180°).
pub const ANTI_ALIGNED_CENTER_RAD: f64 = std::f64::consts::PI;
pub const ANTI_ALIGNED_CENTER_DEG: f64 = 180.0;

/// The islands occur at **high eccentricity** (the abstract's phrase). Beust's
/// portraits place the anti-aligned libration zone above e ≈ 0.6 for the
/// distant (a ≳ 200 AU) test particles relevant to the ETNO sample.
pub const HIGH_ECCENTRICITY_THRESHOLD: f64 = 0.6;

/// Representative distant-ETNO semi-major axis the paper's portraits target
/// (the clustered objects have a from ~250 to ~500 AU). Used only as a
/// portrait location, not as an answer.
pub const REPRESENTATIVE_ETNO_A_AU: f64 = 250.0;

/// Convenience: anti-aligned center in radians from a degree literal, so the
/// 180° figure and the π constant are pinned to be consistent.
pub fn anti_aligned_center_rad_from_deg() -> f64 {
    ANTI_ALIGNED_CENTER_DEG * DEG2RAD
}
