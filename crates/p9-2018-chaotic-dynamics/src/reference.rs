//! Published reference constants from Hadden, Li, Payne & Holman (2018),
//! "Chaotic dynamics of trans-Neptunian objects perturbed by Planet Nine"
//! (arXiv:1712.06547).
//!
//! These are *labelled reference values* for cross-checking the from-scratch
//! computation, never used as computed outputs. Residuals against them are
//! documented in the test modules.

/// Fiducial Planet Nine parameters used throughout the paper's resonance maps:
/// a₉ ≈ 500 AU, m₉ ≈ 10 M⊕, moderate eccentricity. The paper explores the
/// resonance web for "appropriate choices of planet mass and orbit"; the
/// 500 AU / 10 M⊕ point is the canonical illustration.
pub const FIDUCIAL_A9_AU: f64 = 500.0;
/// Fiducial Planet Nine mass (Earth masses).
pub const FIDUCIAL_M9_EARTH: f64 = 10.0;
/// Fiducial Planet Nine eccentricity.
pub const FIDUCIAL_E9: f64 = 0.2;

/// The paper's qualitative headline, stated as a checkable proposition rather
/// than a number: chaotic motion driven by mean-motion-resonance overlap
/// dominates the large-semi-major-axis, high-eccentricity ETNO region, while
/// low-eccentricity orbits remain regular. Encoded as a representative
/// (a, e) point the paper places firmly in the chaotic web: a large-a orbit
/// (close to Planet Nine's overlap zone) at high eccentricity.
pub const CHAOTIC_EXAMPLE_A_AU: f64 = 480.0;
pub const CHAOTIC_EXAMPLE_E: f64 = 0.85;

/// A representative low-eccentricity point the paper places in the regular
/// (non-overlapped) regime at the same semi-major axis.
pub const REGULAR_EXAMPLE_A_AU: f64 = 430.0;
pub const REGULAR_EXAMPLE_E: f64 = 0.05;
