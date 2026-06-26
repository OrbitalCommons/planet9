//! Mini-Neptune composition endpoints for Planet Nine (Russell & White 2025).
//!
//! The paper's most-likely composition is a mini-Neptune: a rock/ice core under
//! a thin H-He envelope (0.6%–3.5% by mass), giving a radius of 2.0–2.6 Earth
//! radii for the Brown et al. (2024) mass prior of 6.6 (+2.6/−1.7) M⊕.
//!
//! The two endpoints are anchored to a geometric albedo extrapolated from the
//! Solar System giant planets: the smaller, warmer 2.0 R⊕ end is the *darker*
//! p = 0.33 case, and the larger, colder 2.6 R⊕ end is the *brighter* p = 0.47
//! case (more reflective clouds at lower equilibrium temperature). This
//! ordering — larger radius ↔ higher albedo — is what reproduces the published
//! absolute-magnitude range (verified in [`crate::photometry`]).

/// A single mini-Neptune composition endpoint: radius, geometric V-band albedo,
/// and the H-He envelope mass fraction that produces that radius.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Composition {
    /// Radius in Earth radii.
    pub radius_earth: f64,
    /// Geometric albedo in the V band (extrapolated from the giant planets).
    pub albedo: f64,
    /// H-He envelope fraction by mass (dimensionless, e.g. 0.006 = 0.6%).
    pub envelope_fraction: f64,
}

/// Smaller / warmer mini-Neptune endpoint: 2.0 R⊕, the *darker* (p = 0.33)
/// case with the thinner 0.6%-by-mass H-He envelope. This is the *faint*
/// (least negative H) end of the published absolute-magnitude range.
pub const SMALL_WARM: Composition = Composition {
    radius_earth: 2.0,
    albedo: 0.33,
    envelope_fraction: 0.006,
};

/// Larger / colder mini-Neptune endpoint: 2.6 R⊕, the *brighter* (p = 0.47)
/// case with the thicker 3.5%-by-mass H-He envelope. This is the *bright*
/// (most negative H) end of the published absolute-magnitude range.
pub const LARGE_COLD: Composition = Composition {
    radius_earth: 2.6,
    albedo: 0.47,
    envelope_fraction: 0.035,
};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn endpoints_bracket_the_published_radius_range() {
        assert_eq!(SMALL_WARM.radius_earth, 2.0);
        assert_eq!(LARGE_COLD.radius_earth, 2.6);
        assert!(SMALL_WARM.radius_earth < LARGE_COLD.radius_earth);
    }

    #[test]
    fn larger_radius_pairs_with_higher_albedo() {
        // The mapping that reproduces the H range: bigger, colder → brighter.
        assert!(LARGE_COLD.albedo > SMALL_WARM.albedo);
        assert_eq!(SMALL_WARM.albedo, 0.33);
        assert_eq!(LARGE_COLD.albedo, 0.47);
    }

    #[test]
    fn envelope_fractions_match_published_bounds() {
        // 0.6%–3.5% H-He by mass; thicker envelope inflates the radius.
        assert!((SMALL_WARM.envelope_fraction - 0.006).abs() < 1e-12);
        assert!((LARGE_COLD.envelope_fraction - 0.035).abs() < 1e-12);
        assert!(LARGE_COLD.envelope_fraction > SMALL_WARM.envelope_fraction);
    }
}
