//! Headline reproduction of Benakli (2026), "A Solar-System Window for Hidden
//! Stellar Companions" (arXiv:2606.25190).
//!
//! Two claims, both recomputed (never hard-coded):
//!
//! 1. The tidal mass–distance envelope grows as M_max(d) ∝ d³, and with the
//!    Planet-Nine-like anchor (5 M⊕ at 500 AU) opens the documented window:
//!    Earth-scale at 300 AU, sub-Saturn at 1000 AU, ~1 Jupiter near 2000 AU.
//! 2. A smooth dark-matter halo holds only a sub-Pluto mass within 1000 AU, so
//!    it cannot supply a nearby planet-scale companion.
//!
//! Computed values (this build):
//!   M_max(300)  ≈ 1.08 M⊕      (Earth-scale)
//!   M_max(1000) ≈ 40 M⊕        (sub-Saturn, 30–95 M⊕)
//!   M_max(2000) ≈ 320 M⊕ ≈ 1.006 M_J  (Jovian)
//!   M_DM(<1000 AU) ≈ 9.5e21 kg < M_Pluto ≈ 1.3e22 kg (sub-Pluto, ~0.73 M_P)

use p9_2026_stellar_companions::{
    dm_halo_mass_within_kg, mass_envelope_earth, ANCHOR_DISTANCE_AU, ANCHOR_MASS_EARTH, M_PLUTO_KG,
};
use p9_core::constants::{EARTH_MASS_SOLAR, MASS_JUPITER_SOLAR};

/// One Jupiter mass expressed in Earth masses, from p9-core mass ratios.
fn jupiter_in_earth_masses() -> f64 {
    MASS_JUPITER_SOLAR / EARTH_MASS_SOLAR
}

#[test]
fn envelope_grows_as_distance_cubed() {
    // Ratio test: M(2d)/M(d) = 8 independent of the base distance.
    for &d in &[300.0_f64, 500.0, 1000.0] {
        let ratio = mass_envelope_earth(2.0 * d) / mass_envelope_earth(d);
        assert!(
            (ratio - 8.0).abs() < 1e-9,
            "M(2·{d})/M({d}) = {ratio:.6}, expected 8 (d³ scaling)"
        );
    }
    // And the anchor is reproduced exactly.
    assert!((mass_envelope_earth(ANCHOR_DISTANCE_AU) - ANCHOR_MASS_EARTH).abs() < 1e-9);
}

#[test]
fn window_is_earth_scale_at_300_au() {
    let m = mass_envelope_earth(300.0);
    assert!(
        (0.5..=3.0).contains(&m),
        "M_max(300 AU) = {m:.3} M⊕, expected Earth-scale"
    );
}

#[test]
fn window_is_sub_saturn_at_1000_au() {
    // Sub-Saturn ≈ 0.1–0.3 M_J ≈ 30–95 M⊕.
    let m = mass_envelope_earth(1000.0);
    assert!(
        (30.0..=95.0).contains(&m),
        "M_max(1000 AU) = {m:.1} M⊕, expected sub-Saturn (30–95 M⊕)"
    );
}

#[test]
fn window_is_jovian_near_2000_au() {
    let m_earth = mass_envelope_earth(2000.0);
    let m_jup = m_earth / jupiter_in_earth_masses();
    assert!(
        (0.8..=1.2).contains(&m_jup),
        "M_max(2000 AU) = {m_earth:.0} M⊕ = {m_jup:.3} M_J, expected ~1 Jovian"
    );
}

#[test]
fn dark_matter_within_1000_au_is_sub_pluto() {
    let m_dm = dm_halo_mass_within_kg(1000.0);
    assert!(
        m_dm < M_PLUTO_KG,
        "smooth DM within 1000 AU = {m_dm:e} kg, must be below Pluto = {M_PLUTO_KG:e} kg"
    );
    // The paper's "sub-Pluto" claim: the enclosed smooth-halo mass is a
    // fraction of Pluto's, not a planet-scale body.
    let fraction = m_dm / M_PLUTO_KG;
    assert!(
        (0.5..1.0).contains(&fraction),
        "DM(<1000 AU) = {fraction:.3} M_Pluto, expected sub-Pluto (< 1)"
    );
}
