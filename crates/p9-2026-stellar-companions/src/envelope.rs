//! Tidal mass–distance envelope for a hidden distant companion.
//!
//! A distant perturber of mass `M` at heliocentric distance `d` raises a
//! tidal/quadrupole acceleration across the inner Solar System that scales as
//!
//! ```text
//!     a_tidal  ∝  G·M / d³
//! ```
//!
//! (the differential of the `G·M/d²` direct term over the size of the inner
//! system). Benakli (2026) calibrates the maximum allowed companion mass by
//! holding this tidal amplitude `G·M/d³` fixed at the level still permitted by
//! the Planet-Nine-like ephemeris constraint. Holding `M/d³` constant means the
//! envelope grows as the cube of distance:
//!
//! ```text
//!     M_max(d) = M_anchor · (d / d_anchor)³
//! ```
//!
//! The anchor `(M_anchor, d_anchor)` is the Planet-Nine dynamical point: this
//! crate uses `5 M⊕ at 500 AU`, the Brown & Batygin (2019) revised best fit
//! carried in `p9_core` as [`P9Params::revised_2019`](p9_core::types::P9Params).

use crate::{ANCHOR_DISTANCE_AU, ANCHOR_MASS_EARTH};

/// Maximum companion mass (Earth masses) permitted by the tidal envelope at a
/// heliocentric distance `d_au` (AU).
///
/// `M_max(d) = M_anchor · (d / d_anchor)³`, with the Planet-Nine-like anchor of
/// [`ANCHOR_MASS_EARTH`] at [`ANCHOR_DISTANCE_AU`]. The cubic growth is the
/// direct consequence of holding the tidal amplitude `G·M/d³` fixed.
pub fn mass_envelope_earth(d_au: f64) -> f64 {
    let ratio = d_au / ANCHOR_DISTANCE_AU;
    ANCHOR_MASS_EARTH * ratio * ratio * ratio
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn passes_through_the_anchor() {
        assert_relative_eq!(
            mass_envelope_earth(ANCHOR_DISTANCE_AU),
            ANCHOR_MASS_EARTH,
            epsilon = 1e-12
        );
    }

    #[test]
    fn grows_as_distance_cubed() {
        // Doubling distance must multiply the allowed mass by 2³ = 8, at any
        // base distance — the defining property of the envelope.
        for &d in &[300.0_f64, 550.0, 1234.0] {
            let r = mass_envelope_earth(2.0 * d) / mass_envelope_earth(d);
            assert_relative_eq!(r, 8.0, epsilon = 1e-9);
        }
    }

    #[test]
    fn monotonic_in_distance() {
        let mut last = 0.0;
        for d in (200..=2000).step_by(100) {
            let m = mass_envelope_earth(d as f64);
            assert!(m > last, "envelope must increase with distance");
            last = m;
        }
    }
}
