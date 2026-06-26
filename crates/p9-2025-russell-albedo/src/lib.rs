//! Reproduction of Russell & White (2025), "The Radius, Composition, Albedo,
//! and Absolute Magnitude of Planet Nine Based on Exoplanets with Teq < 600 K
//! and the Planet Nine Reference Population 3.0" (arXiv:2507.22297).
//!
//! # What the paper does
//!
//! Starting from the Brown et al. (2024) mass prior for Planet Nine —
//! 6.6 (+2.6/−1.7) Earth masses — the authors use the population of cold
//! (Teq < 600 K) exoplanets to argue Planet Nine is most likely a **mini-
//! Neptune**: a rock/ice core under a thin H-He envelope (0.6%–3.5% by mass),
//! giving a radius of **2.0–2.6 R⊕**. They extrapolate the geometric V-band
//! albedo from the Solar System giant planets (Planet Nine being colder and
//! more reflective at the larger-radius end), obtaining **p = 0.33 → 0.47**.
//! Folding radius and albedo through the reflected-light absolute-magnitude
//! law gives **H(V) = −6.1 → −5.2**, and placing the planet at the most-likely
//! orbit of Reference Population 3.0 gives an apparent magnitude of
//! **+21.9 → +22.7**.
//!
//! # Physics reproduced here
//!
//! Pure reflected-light photometry, built directly on
//! `p9_core::analysis::photometry` (no reimplementation, no wrappers):
//!
//! - [`absolute_magnitude`](p9_core::analysis::photometry::absolute_magnitude)
//!   — H = 5 log10(1329 km / (√p · 2R_km)).
//! - [`apparent_magnitude`](p9_core::analysis::photometry::apparent_magnitude)
//!   — m = H + 5 log10(r·Δ), the 1/(r²Δ²) reflected-sunlight law.
//! - [`opposition_delta`](p9_core::analysis::photometry::opposition_delta)
//!   — Δ = r − 1 AU.
//! - [`EARTH_RADIUS_KM`](p9_core::constants::EARTH_RADIUS_KM) for the
//!   R⊕ → km conversion.
//!
//! # Composition endpoints and the albedo ordering
//!
//! [`composition`] pins the two mini-Neptune endpoints. The mapping that
//! reproduces the published H range is **larger radius ↔ higher albedo**:
//!
//! | endpoint              | R (R⊕) | albedo p | H-He frac | H(V)   |
//! |-----------------------|--------|----------|-----------|--------|
//! | [`SMALL_WARM`]        | 2.0    | 0.33     | 0.6%      | ≈ −5.2 |
//! | [`LARGE_COLD`]        | 2.6    | 0.47     | 3.5%      | ≈ −6.1 |
//!
//! so the **bigger, colder, brighter-albedo** endpoint is the **brighter**
//! (most negative H) one, and the smaller/warmer/darker endpoint is the faint
//! end of the range.
//!
//! # Most-likely distance
//!
//! The apparent magnitude depends on the current heliocentric distance. The
//! paper quotes m = +21.9..+22.7 at the most-likely orbit of Reference
//! Population 3.0. Solving the reflected-light law for the single distance that
//! places the bright endpoint near m = 21.9 and the faint endpoint near
//! m = 22.7 gives [`photometry::MOST_LIKELY_DISTANCE_AU`] ≈ 625 AU (a few-
//! hundred-AU heliocentric distance, consistent with the most-likely orbit's
//! outer reach).

pub mod composition;
pub mod photometry;

pub use composition::{Composition, LARGE_COLD, SMALL_WARM};
pub use photometry::{absolute_h, apparent_v, apparent_v_most_likely, MOST_LIKELY_DISTANCE_AU};

/// Published Planet Nine mass prior (Brown et al. 2024): central value in Earth
/// masses. Reference constant only.
pub const PUBLISHED_MASS_EARTH: f64 = 6.6;
/// Published mass prior upper 1σ offset (+2.6 M⊕). Reference constant only.
pub const PUBLISHED_MASS_PLUS: f64 = 2.6;
/// Published mass prior lower 1σ offset (−1.7 M⊕). Reference constant only.
pub const PUBLISHED_MASS_MINUS: f64 = 1.7;

/// Published mini-Neptune radius range, lower bound (R⊕). Reference only.
pub const PUBLISHED_RADIUS_MIN_EARTH: f64 = 2.0;
/// Published mini-Neptune radius range, upper bound (R⊕). Reference only.
pub const PUBLISHED_RADIUS_MAX_EARTH: f64 = 2.6;

/// Published H-He envelope mass fraction, lower bound (0.6%). Reference only.
pub const PUBLISHED_ENVELOPE_FRAC_MIN: f64 = 0.006;
/// Published H-He envelope mass fraction, upper bound (3.5%). Reference only.
pub const PUBLISHED_ENVELOPE_FRAC_MAX: f64 = 0.035;

/// Published geometric V-band albedo, lower bound. Reference only.
pub const PUBLISHED_ALBEDO_MIN: f64 = 0.33;
/// Published geometric V-band albedo, upper bound. Reference only.
pub const PUBLISHED_ALBEDO_MAX: f64 = 0.47;

/// Published absolute magnitude H(V), bright (most negative) bound.
/// Reference only.
pub const PUBLISHED_H_BRIGHT: f64 = -6.1;
/// Published absolute magnitude H(V), faint (least negative) bound.
/// Reference only.
pub const PUBLISHED_H_FAINT: f64 = -5.2;

/// Published apparent V magnitude, bright (brightest) bound. Reference only.
pub const PUBLISHED_APPARENT_BRIGHT: f64 = 21.9;
/// Published apparent V magnitude, faint (faintest) bound. Reference only.
pub const PUBLISHED_APPARENT_FAINT: f64 = 22.7;
