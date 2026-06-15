//! Reproduction of Holman & Payne (2016), "Observational constraints on Planet
//! Nine: Cassini range observations of Saturn" (arXiv:1603.09008 part I and
//! arXiv:1604.03180 part II).
//!
//! # Headline result
//!
//! A Planet Nine on the Brown & Batygin (2016) orbit tugs gravitationally on
//! Saturn relative to the Sun. The Cassini spacecraft measured the Earth–Saturn
//! range to a precision of order tens of metres over the mission; a P9 close
//! enough and massive enough would imprint an anomalous range signal larger than
//! that precision and is therefore EXCLUDED, while a distant/light P9 produces a
//! sub-threshold signal and is ALLOWED. Holman & Payne map this into the
//! `(mass, heliocentric-distance, true-anomaly)` space: large regions are
//! excluded, and a localized zone — on the near side of the orbit, toward the
//! sky position `RA ≈ 40°, Dec ≈ −15°` (±~20°) — is favored, complementary to
//! and consistent with Fienga et al. (2016).
//!
//! # What this crate computes (real, no hard-coded answers)
//!
//! This is the SAME differential-acceleration / range-signal physics as the
//! sibling [`p9_2016_cassini_ranging`] crate, which we REUSE directly rather
//! than reimplement. On top of it this crate adds the Holman & Payne framing:
//!
//! 1. [`signal`]: reduce a P9 at `(mass, distance via a, true anomaly ν)` to a
//!    single Earth–Saturn range-residual AMPLITUDE (m), via the sibling's
//!    Cassini-arc range signature. Demonstrates the published SCALINGS — linear
//!    in P9 mass, ~`1/d³` in heliocentric distance.
//! 2. [`exclusion`]: compare that amplitude to a documented Cassini range
//!    precision ([`published::CASSINI_RANGE_PRECISION_M`]) to decide
//!    EXCLUDED / ALLOWED, and sweep the `(mass, distance)` plane into an
//!    exclusion map (the paper's mass–distance constraint).
//! 3. [`sky`]: the favored P9 sky direction on the near (post-perihelion) side
//!    of the B&B orbit, reported as ecliptic-derived RA/Dec, compared to the
//!    paper's `RA ≈ 40°, Dec ≈ −15°` preferred region.
//!
//! # Approximations (analytic proxy, not the full Cassini-data fit)
//!
//! Holman & Payne fit the actual Cassini range residual stream while floating
//! Saturn's orbit; we use the sibling's quasi-static analytic proxy (double
//! time-integration of the differential acceleration over Saturn's real
//! Cassini-epoch arc, projected onto the line of sight, with the low-order
//! polynomial the orbit fit absorbs removed). This reproduces the SHAPE of the
//! constraint (a near-side favored zone, perihelion-facing exclusion) and the
//! SCALINGS the tests pin; it does not reproduce the paper's absolute residual
//! amplitudes. See [`published`] for the labelled reference constants and the
//! test modules / `REPRODUCTION_NOTES.md` for the honest residual discussion.

pub mod exclusion;
pub mod signal;
pub mod sky;

use p9_core::types::P9Params;

// Reuse the sibling's orbit and geometry verbatim — the two papers analyze the
// same Brown & Batygin (2016) reference orbit and the same Sun–Saturn geometry.
pub use p9_2016_cassini_ranging::brown_batygin_orbit;
pub use p9_2016_cassini_ranging::geometry::{p9_position_at_true_anomaly, P9Geometry};

/// Planet Nine on the Brown & Batygin (2016) orbit used throughout this paper.
/// Alias of the sibling crate's [`brown_batygin_orbit`]; provided so callers can
/// build the orbit without naming the sibling crate.
pub fn holman_payne_orbit() -> P9Params {
    brown_batygin_orbit()
}

/// Published reference values from Holman & Payne (2016), kept as labelled
/// constants (NOT used to derive any computed quantity).
pub mod published {
    /// Planet Nine mass assumed by the Brown & Batygin orbit the paper analyzes
    /// (Earth masses).
    pub const P9_MASS_EARTH: f64 = 10.0;

    /// Representative Cassini Earth–Saturn range precision (metres). Cassini
    /// range residuals are at the tens-of-metres level over the mission; Holman
    /// & Payne (2016) treat an induced range signal above this floor as a
    /// detectable (hence excludable) perturbation. Used only as the
    /// detection threshold in [`crate::exclusion`], never to derive a signal.
    pub const CASSINI_RANGE_PRECISION_M: f64 = 30.0;

    /// Preferred current sky location of P9 from the Cassini analysis
    /// (arXiv:1604.03180): `RA ≈ 40°`, extending ~20° in all directions.
    pub const PREFERRED_RA_DEG: f64 = 40.0;

    /// Preferred current sky location of P9: `Dec ≈ −15°`, ±~20°.
    pub const PREFERRED_DEC_DEG: f64 = -15.0;

    /// Half-extent of the preferred sky region (degrees): "~20° in all
    /// directions" (arXiv:1604.03180).
    pub const PREFERRED_SKY_HALF_EXTENT_DEG: f64 = 20.0;

    /// Most-probable true anomaly of P9 on the near (post-perihelion) side of
    /// the B&B orbit that the range data favor (degrees). Consistent with the
    /// Fienga et al. (2016) `ν ≈ 117.8°` zone the companion crate pins.
    pub const PREFERRED_TRUE_ANOMALY_DEG: f64 = 117.8;

    /// Tolerance (deg) within which the analytic-proxy favored anomaly should
    /// land of the published preferred value.
    pub const TRUE_ANOMALY_TOLERANCE_DEG: f64 = 30.0;
}

pub use exclusion::{excluded, max_allowed_mass_earth, ExclusionMap, ExclusionVerdict};
pub use signal::{favored_true_anomaly, range_residual_m};
pub use sky::{favored_sky_position, SkyPosition};
