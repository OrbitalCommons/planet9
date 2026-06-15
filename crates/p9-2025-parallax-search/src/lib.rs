//! Reproduction of Socas-Navarro & Trujillo (2025),
//! *A targeted parallax search for Planet Nine* (arXiv:2504.05473).
//!
//! # The method
//!
//! Planet Nine, if it lies at a heliocentric distance `d` of order 500–700 AU,
//! has a tiny *orbital* (proper) motion but a **huge reflex parallax**. Earth's
//! 1 AU heliocentric displacement subtends, seen from `d`, an angle
//!
//! ```text
//!   p_half = (1 AU / d) * 206265"   (the half-angle, arcsec)
//! ```
//!
//! At 460 AU this is ≈ 7.5 *arcminutes* — four orders of magnitude larger than
//! any background star's parallax (a field star at even 10 pc has a parallax of
//! 0.1″; typical survey stars are at hundreds of pc, ≲ a few milliarcsec). So
//! observing the same field at two epochs separated by a chunk of Earth's orbit
//! (the paper uses consecutive nights up to a ~6-month / quadrature baseline)
//! makes Planet Nine *leap* across the field while the stellar background barely
//! stirs. Combined with P9's slow orbital proper motion (it does not run off in a
//! consistent direction the way a nearby fast-moving asteroid would), the large
//! reflex parallax cleanly separates a bound distant planet from the field.
//!
//! # What this crate computes (all real, none hard-coded)
//!
//! - [`parallax`] — the parallactic displacement of a body at distance `d` over a
//!   chosen epoch separation, in arcmin, reusing `p9_core::coords::parallax`. The
//!   amplitude falls as `1/d`; the full annual peak-to-peak swing and the
//!   half-angle ("parallax factor") are both exposed.
//! - [`discrimination`] — the contrast ratio against a background star of given
//!   parallax (≳ 10³–10⁴×), and the per-epoch astrometric precision required to
//!   detect the shift at a chosen SNR, which *grows linearly with `d`*.
//!
//! Published reference numbers from the paper are kept as labelled constants in
//! [`published`].

pub mod discrimination;
pub mod parallax;

/// Labelled reference numbers from Socas-Navarro & Trujillo (2025) and the
/// physical setup, kept separate from anything this crate computes.
pub mod published {
    /// Lower edge of the P9 heliocentric-distance prior the search targets (AU).
    pub const P9_DISTANCE_MIN_AU: f64 = 500.0;
    /// Upper edge of the P9 heliocentric-distance prior the search targets (AU).
    pub const P9_DISTANCE_MAX_AU: f64 = 700.0;
    /// A representative distance used in the worked examples (AU). The reflex
    /// half-parallax at this distance is ≈ 7.5 arcmin (see tests).
    pub const REFERENCE_DISTANCE_AU: f64 = 460.0;
    /// Reflex parallax half-angle at [`REFERENCE_DISTANCE_AU`], arcmin. This is
    /// the order-of-arcminutes headline number; the crate recomputes it.
    pub const HALF_PARALLAX_AT_REF_ARCMIN: f64 = 7.47;
    /// The targeted survey area (square degrees) imaged in the paper.
    pub const SURVEY_AREA_DEG2: f64 = 98.0;
    /// Mean r-band limiting magnitude of the search.
    pub const LIMITING_MAG_R: f64 = 21.3;
    /// A representative *background field star* parallax (arcsec). Stars used as
    /// astrometric reference in these fields sit at hundreds of pc, giving
    /// parallaxes of a few milliarcsec; 1 mas corresponds to 1 kpc.
    pub const TYPICAL_FIELD_STAR_PARALLAX_ARCSEC: f64 = 1.0e-3;
}
