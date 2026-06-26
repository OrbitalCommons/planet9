//! Solid-angle fraction of encounter orientations that satisfy the
//! inclination / apsidal-alignment constraints.
//!
//! The paper's N-body experiments show that a flyby only reproduces the
//! observed detached-ETNO population — **low inclinations (i < 30°) AND
//! primordial apsidal alignment** — for a narrow set of encounter geometries:
//!
//!   1. a near-**coplanar** encounter, inclination of the stellar orbit
//!      `i_* ≈ 0°` relative to the planetary plane, OR
//!   2. an **ecliptic-symmetric** encounter, `i_* ≈ 90°` with argument of
//!      perihelion `ω_* ≈ 0°`, so the perturbation is symmetric about the
//!      ecliptic and does not pump inclinations.
//!
//! Random encounters in a cluster arrive isotropically: the stellar orbit
//! orientation (i_*, ω_*) is uniform over the sphere. We compute the fraction
//! of that 4π solid angle that falls within a tolerance band of either special
//! geometry. This fraction is the geometric "filter" the encounter must pass.
//!
//! Orientation measure. For an isotropic orientation, `cos i_*` is uniform on
//! [−1, 1] and `ω_*` is uniform on [0, 2π). We integrate the indicator of the
//! two allowed bands against that measure.

use std::f64::consts::PI;

/// Half-width tolerance (radians) on the coplanar inclination band, `i_* ≈ 0`.
/// A few-degree window: outside it the flyby over-excites inclinations past the
/// observed i < 30° envelope. We take 10° as a generous acceptance band.
pub const COPLANAR_TOL_RAD: f64 = 10.0 * PI / 180.0;

/// Half-width tolerance (radians) on the ecliptic-symmetric band: `i_* ≈ 90°`
/// AND `ω_* ≈ 0°` (mod π). 10° in each of the two angles.
pub const SYMMETRIC_TOL_RAD: f64 = 10.0 * PI / 180.0;

/// Fraction of isotropic orientations within `tol` of the coplanar pole
/// (`i_* ≈ 0`, either prograde i_*≈0 or retrograde i_*≈180°).
///
/// For uniform `cos i`, the fraction with `i < tol` (the polar cap) is
/// `(1 − cos tol)`; counting both poles gives `2(1 − cos tol)` over the full
/// range of `cos i ∈ [−1, 1]` (normalised by 2).
pub fn coplanar_fraction(tol: f64) -> f64 {
    // Cap near cos i = +1: measure (1 − cos tol). Cap near cos i = −1: same.
    // Total measure of cos i is 2, so fraction = 2(1 − cos tol)/2 = (1 − cos tol).
    1.0 - tol.cos()
}

/// Fraction of isotropic orientations within `tol` of the ecliptic-symmetric
/// configuration: `i_* ≈ 90°` (a band in `cos i` of half-width `sin tol`) AND
/// `ω_* ≈ 0` (mod π), a band in `ω_*` of fractional width `2·(2 tol)/(2π)`.
pub fn symmetric_fraction(tol: f64) -> f64 {
    // |cos i| < sin(tol) → measure 2 sin(tol) out of 2 → fraction sin(tol).
    let i_band = tol.sin();
    // ω within ±tol of 0 or π → 2 windows of width 2·tol out of 2π.
    let omega_band = 2.0 * (2.0 * tol) / (2.0 * PI);
    i_band * omega_band
}

/// Total fraction of isotropic encounter orientations that satisfy the
/// inclination/alignment constraint (coplanar OR ecliptic-symmetric).
///
/// The two bands are disjoint (one near i_*=0/180°, the other near i_*=90°),
/// so the fraction is the sum.
pub fn geometry_fraction() -> f64 {
    coplanar_fraction(COPLANAR_TOL_RAD) + symmetric_fraction(SYMMETRIC_TOL_RAD)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fraction_is_small_solid_angle() {
        let f = geometry_fraction();
        assert!((0.0..1.0).contains(&f), "fraction out of range: {f}");
        // It is a *small* slice of the sphere — well under 10%.
        assert!(f < 0.10, "geometry fraction should be small, got {f:.4}");
    }

    #[test]
    fn coplanar_band_grows_with_tolerance() {
        assert!(coplanar_fraction(0.2) > coplanar_fraction(0.1));
    }

    #[test]
    fn both_bands_contribute() {
        assert!(coplanar_fraction(COPLANAR_TOL_RAD) > 0.0);
        assert!(symmetric_fraction(SYMMETRIC_TOL_RAD) > 0.0);
    }

    #[test]
    fn isotropic_full_sphere_normalisation() {
        // As tol → 90°, the coplanar caps cover the whole sphere: fraction → 1.
        assert!((coplanar_fraction(PI / 2.0) - 1.0).abs() < 1e-12);
    }
}
