//! Photometry of distant solar-system bodies: mass-radius relation, absolute
//! magnitude, the reflected-sunlight apparent-magnitude law, and the IAU
//! H-G phase function.
//!
//! Single source for the workspace (p9-2016-constraints, p9-2019-review,
//! p9-2021-orbit, p9-2022-des and p9-2024-panstarrs previously carried
//! inconsistent copies, one using the stellar 5·log10(d/10) distance law).
//!
//! Reflected sunlight falls off as 1/(r²Δ²) — once for the Sun→body leg, once
//! for the body→Earth leg — so
//!
//!   m = H + 5 log10(r·Δ) − 2.5 log10 Φ(α),
//!
//! never the stellar 5 log10(d/10 pc). The phase function Φ(α) is the
//! two-term IAU H-G law of Bowell et al. (1989, in *Asteroids II*, 524),
//! implemented locally ([`hg_phase_factor`]) with the same coefficients as
//! starfield-mpc's `hg_apparent_magnitude` (unpublished crate; issue #33
//! tracks contributing this generic-body photometry upstream, so the APIs
//! are kept convergent).
//!
//! Note the H-G law's opposition surge: dΦ/dα → −∞ as α → 0 (the tan^0.63
//! term), so even the tiny phase angles reachable at Planet Nine distances
//! (α ≤ asin(1 AU/r) ≈ 0.1–0.2°) darken the body by a few hundredths of a
//! magnitude relative to the exact-opposition limit — not the µmag a linear
//! phase law would suggest.

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

/// Apparent magnitude of reflected sunlight in the exact-opposition
/// (zero-phase) limit:
///
///   m = H + 5 log10(r·Δ)
///
/// with `r_au` the heliocentric and `delta_au` the geocentric distance.
/// This is [`apparent_magnitude_with_phase`] at α = 0 (Φ(0) = 1); use the
/// phase-aware form whenever the observer geometry yields a phase angle.
pub fn apparent_magnitude(h: f64, r_au: f64, delta_au: f64) -> f64 {
    h + 5.0 * (r_au * delta_au).log10()
}

/// Default IAU H-G slope parameter G = 0.15, the standard assumption for
/// bodies with unmeasured phase behavior (Bowell et al. 1989; the MPC uses
/// the same default for asteroids without fitted G).
pub const DEFAULT_SLOPE_G: f64 = 0.15;

/// IAU H-G phase function Φ(α) (Bowell et al. 1989):
///
///   Φ(α) = (1 − G) Φ₁(α) + G Φ₂(α),
///   Φᵢ(α) = exp(−Aᵢ tan^{Bᵢ}(α/2)),
///   A₁ = 3.332, B₁ = 0.631, A₂ = 1.862, B₂ = 1.218,
///
/// with `phase_rad` the Sun–body–observer angle in radians (valid for
/// 0 ≤ α ≲ 120°). Φ(0) = 1 by construction; Φ decreases monotonically
/// with α. Coefficients match starfield-mpc's `hg_apparent_magnitude`.
pub fn hg_phase_factor(g: f64, phase_rad: f64) -> f64 {
    let tan_half = (phase_rad / 2.0).tan();
    let phi1 = (-3.332 * tan_half.powf(0.631)).exp();
    let phi2 = (-1.862 * tan_half.powf(1.218)).exp();
    (1.0 - g) * phi1 + g * phi2
}

/// Phase integral of the H-G system, q(G) = 0.290 + 0.684 G (Bowell et al.
/// 1989), relating Bond albedo to geometric albedo: A_B = q · p.
pub fn phase_integral(g: f64) -> f64 {
    0.290 + 0.684 * g
}

/// Apparent magnitude of reflected sunlight with the IAU H-G phase term:
///
///   m = H + 5 log10(r·Δ) − 2.5 log10 Φ(α; G)
///
/// `phase_rad` is the Sun–body–observer phase angle (see [`phase_angle`]);
/// `g` is the H-G slope parameter ([`DEFAULT_SLOPE_G`] = 0.15 when
/// unmeasured). Reduces to [`apparent_magnitude`] at α = 0.
pub fn apparent_magnitude_with_phase(
    h: f64,
    r_au: f64,
    delta_au: f64,
    phase_rad: f64,
    g: f64,
) -> f64 {
    h + 5.0 * (r_au * delta_au).log10() - 2.5 * hg_phase_factor(g, phase_rad).log10()
}

/// Phase angle α (Sun–body–observer, radians) from the triangle of
/// heliocentric distance `r_helio`, body–observer distance `delta`, and
/// Sun–observer distance `r_earth_sun` (law of cosines):
///
///   cos α = (r² + Δ² − R²) / (2 r Δ)
///
/// All three sides in the same units (AU). At Planet Nine distances the
/// phase angle is bounded by asin(R/r) ≈ 0.1–0.2°.
pub fn phase_angle(r_helio: f64, delta: f64, r_earth_sun: f64) -> f64 {
    let cos_alpha =
        (r_helio * r_helio + delta * delta - r_earth_sun * r_earth_sun) / (2.0 * r_helio * delta);
    cos_alpha.clamp(-1.0, 1.0).acos()
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
    fn test_hg_phase_factor_is_one_at_opposition() {
        // Φ(0) = 1 exactly: tan(0)^B = 0, exp(0) = 1, for any G.
        for &g in &[0.0, 0.15, 0.25, 1.0] {
            assert_eq!(hg_phase_factor(g, 0.0), 1.0, "G = {g}");
        }
        // ... so the phase-aware law reduces to the opposition form.
        let (h, r, d) = (-6.87, 30.1, 29.1);
        assert_eq!(
            apparent_magnitude_with_phase(h, r, d, 0.0, DEFAULT_SLOPE_G),
            apparent_magnitude(h, r, d)
        );
    }

    #[test]
    fn test_hg_phase_factor_monotone_decreasing() {
        for &g in &[0.0, 0.15, 0.4] {
            let mut prev = 1.0;
            for k in 1..=120 {
                let alpha = (k as f64).to_radians();
                let f = hg_phase_factor(g, alpha);
                assert!(f < prev, "Φ not decreasing at α = {k}° (G = {g})");
                assert!(f > 0.0);
                prev = f;
            }
        }
    }

    #[test]
    fn test_hg_phase_factor_pinned_values() {
        // Direct evaluation of the Bowell et al. (1989) two-term law with
        // G = 0.15 (independent reference computation, matching
        // starfield-mpc's coefficients):
        //   α = 10°: Φ = 0.55159, Δm = 0.6460 mag
        //   α = 30°: Φ = 0.30225, Δm = 1.2991 mag
        let f10 = hg_phase_factor(0.15, 10.0_f64.to_radians());
        assert!((f10 - 0.551585).abs() < 1e-5, "Φ(10°) = {f10:.6}");
        let f30 = hg_phase_factor(0.15, 30.0_f64.to_radians());
        assert!((f30 - 0.302247).abs() < 1e-5, "Φ(30°) = {f30:.6}");

        // Ceres-like opposition case from starfield-mpc's test suite:
        // H = 3.35, G = 0.15, r = 2.77, Δ = 1.77, α = 0 → V = 6.802.
        let v = apparent_magnitude_with_phase(3.35, 2.77, 1.77, 0.0, 0.15);
        assert!((v - 6.802).abs() < 0.001, "V_ceres = {v:.3}");
    }

    #[test]
    fn test_phase_integral_standard_values() {
        // q(G) = 0.290 + 0.684 G: q(0.15) = 0.3926.
        assert!((phase_integral(0.15) - 0.3926).abs() < 1e-12);
        assert!((phase_integral(0.0) - 0.290).abs() < 1e-12);
    }

    #[test]
    fn test_phase_angle_geometry() {
        // Opposition (Sun–Earth–body collinear): Δ = r − R → α = 0.
        assert!(phase_angle(500.0, 499.0, 1.0) < 1e-7);
        // Conjunction (body behind the Sun): Δ = r + R → α = 0 too.
        assert!(phase_angle(500.0, 501.0, 1.0) < 1e-7);
        // Quadrature right triangle at the *body*: r² + Δ² = R² is
        // impossible for r > R; instead check the right angle at the Sun:
        // Δ² = r² + R² → cos α = r/Δ.
        let r = 400.0;
        let delta = (r * r + 1.0_f64).sqrt();
        let alpha = phase_angle(r, delta, 1.0);
        assert!((alpha.cos() - r / delta).abs() < 1e-12);
        // The phase angle is bounded by the annual maximum asin(R/r).
        let alpha_max = (1.0_f64 / r).asin();
        assert!(alpha <= alpha_max + 1e-12, "α = {alpha} > {alpha_max}");
    }

    #[test]
    fn test_phase_darkening_small_at_p9_distances() {
        // At r = 300–700 AU the phase angle never exceeds asin(1/r)
        // (0.08–0.19°), but the H-G opposition surge still produces a few
        // hundredths of a magnitude: Δm(α_max) = 0.054 mag at 300 AU,
        // 0.039 mag at 500 AU. (A linear 0.04 mag/deg law would give
        // ~5 mmag — the surge dominates at these tiny angles.)
        for &(r, dm_expect) in &[(300.0, 0.0542), (500.0, 0.0393), (700.0, 0.0318)] {
            let alpha_max = (1.0_f64 / r).asin();
            let dm = -2.5 * hg_phase_factor(DEFAULT_SLOPE_G, alpha_max).log10();
            assert!(
                (dm - dm_expect).abs() < 5e-3,
                "Δm(α_max) at {r} AU = {dm:.4}, expected ~{dm_expect}"
            );
        }
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
