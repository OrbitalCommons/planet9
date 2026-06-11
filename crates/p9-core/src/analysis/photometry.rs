//! Photometry of distant solar-system bodies: mass-radius relation, absolute
//! magnitude, and the reflected-sunlight apparent-magnitude law.
//!
//! Single source for the workspace (p9-2016-constraints, p9-2019-review,
//! p9-2021-orbit, p9-2022-des and p9-2024-panstarrs previously carried
//! inconsistent copies, one using the stellar 5·log10(d/10) distance law).
//!
//! Reflected sunlight falls off as 1/(r²Δ²) — once for the Sun→body leg, once
//! for the body→Earth leg — so
//!
//!   m = H + 5 log10(r·Δ) + phase terms (≈ 0 near opposition at large r),
//!
//! never the stellar 5 log10(d/10 pc).

/// Neptune's volumetric mean radius in Earth radii (24 622 km / 6 371 km).
pub const NEPTUNE_RADIUS_EARTH: f64 = 3.865;

/// Neptune's mass in Earth masses.
pub const NEPTUNE_MASS_EARTH: f64 = 17.147;

/// Neptune's geometric albedo (V band).
pub const ALBEDO_NEPTUNE: f64 = 0.41;

/// Mass-radius relation for Neptunian (volatile-envelope) planets, anchored
/// at Neptune:
///
///   R/R⊕ = 3.865 · (M / 17.147 M⊕)^0.27
///
/// giving R ≈ 3.4 R⊕ at 10 M⊕ and R ≈ 2.6 R⊕ at 5 M⊕, consistent with the
/// Fortney et al. (2007) ice-giant models (~3.5 R⊕ at 10 M⊕) and the
/// Chen & Kipping (2017) Neptunian-branch slope. Returns Earth radii.
///
/// (The previous inline relations — 3.0·M^0.27, larger than Neptune at
/// 10 M⊕, and 1.0·M^0.27, a 3.5× IR-flux underestimate — bracketed this and
/// flipped detectability conclusions in both directions.)
pub fn mass_radius_neptunian(mass_earth: f64) -> f64 {
    NEPTUNE_RADIUS_EARTH * (mass_earth / NEPTUNE_MASS_EARTH).powf(0.27)
}

/// Absolute magnitude H from radius (km) and geometric albedo:
///
///   H = 5 log10( 1329 km / (√p · D_km) ),  D = 2R.
pub fn absolute_magnitude(radius_km: f64, albedo: f64) -> f64 {
    5.0 * (1329.0 / (albedo.sqrt() * 2.0 * radius_km)).log10()
}

/// Apparent magnitude of reflected sunlight:
///
///   m = H + 5 log10(r·Δ)
///
/// with `r_au` the heliocentric and `delta_au` the geocentric distance
/// (phase function ≈ 1 near opposition for distant bodies).
pub fn apparent_magnitude(h: f64, r_au: f64, delta_au: f64) -> f64 {
    h + 5.0 * (r_au * delta_au).log10()
}

/// Geocentric distance at opposition: Δ = r − 1 AU.
pub fn opposition_delta(r_au: f64) -> f64 {
    r_au - 1.0
}

/// Apparent V magnitude of a planet of mass `mass_earth` and geometric albedo
/// `albedo` at heliocentric distance `r_au`, observed at opposition.
pub fn planet_apparent_magnitude(mass_earth: f64, albedo: f64, r_au: f64) -> f64 {
    let radius_km = mass_radius_neptunian(mass_earth) * crate::constants::EARTH_RADIUS_KM;
    let h = absolute_magnitude(radius_km, albedo);
    apparent_magnitude(h, r_au, opposition_delta(r_au))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_neptune_absolute_magnitude() {
        // Neptune: R = 24 622 km, p_V ≈ 0.41 → H ≈ −6.9 (IAU value −6.87).
        let h = absolute_magnitude(24_622.0, ALBEDO_NEPTUNE);
        assert!((h - (-6.87)).abs() < 0.15, "H_neptune = {h:.2}");
    }

    #[test]
    fn test_neptune_apparent_magnitude() {
        // Neptune at 30.1 AU near opposition: V ≈ 7.8.
        let h = absolute_magnitude(24_622.0, ALBEDO_NEPTUNE);
        let m = apparent_magnitude(h, 30.1, opposition_delta(30.1));
        assert!((m - 7.8).abs() < 0.3, "V_neptune = {m:.2}");
    }

    #[test]
    fn test_mass_radius_anchored_at_neptune() {
        let r = mass_radius_neptunian(NEPTUNE_MASS_EARTH);
        assert!((r - NEPTUNE_RADIUS_EARTH).abs() < 1e-12);
        // 10 M⊕: ~3.3–3.5 R⊕ (Fortney et al. 2007 ice giants)
        let r10 = mass_radius_neptunian(10.0);
        assert!(r10 > 3.0 && r10 < 3.6, "R(10 M⊕) = {r10:.2} R⊕");
        // Smaller than Neptune below Neptune's mass
        assert!(r10 < NEPTUNE_RADIUS_EARTH);
    }

    #[test]
    fn test_p9_faint_at_aphelion() {
        // A 6.2 M⊕ Planet Nine near aphelion (~500 AU) must be far beyond
        // naked-eye brightness — the stellar distance law gave m ≈ 12–15.
        let m = planet_apparent_magnitude(6.2, 0.4, 500.0);
        assert!(m > 20.0 && m < 28.0, "V_P9(500 AU) = {m:.1}");

        // Reflected-light scaling: doubling distance dims by ~3 mag
        // (10 log10 2), not 1.5 mag.
        let m2 = planet_apparent_magnitude(6.2, 0.4, 1000.0);
        let dm = m2 - m;
        assert!((dm - 3.01).abs() < 0.05, "Δm for 2x distance = {dm:.2}");
    }
}
