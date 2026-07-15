//! Giant-planet-induced apsidal precession of a detached test particle.
//!
//! Huang & Gladman (2024) argue that the sednoids were *born* with aligned
//! longitudes of perihelion ϖ and have since been differentially de-aligned by
//! the secular apsidal precession forced by the known giant planets. The
//! engine is purely the four giant planets — no Planet Nine.
//!
//! For a distant (a ≫ a_j) low-perturber-eccentricity test particle, the
//! Laplace–Lagrange secular theory gives the *prograde* apsidal precession rate
//! (Murray & Dermott 1999, Eq. 7.13, the diagonal "A_jj" frequency for an
//! exterior test particle):
//!
//!   ϖ̇ = +(n/4) Σ_j (m_j/M_sun) α_j ᾱ_j b_{3/2}^{(1)}(α_j),   α_j = a_j/a.
//!
//! For an exterior particle ᾱ_j = α_j, and as α_j → 0 the Laplace coefficient
//! b_{3/2}^{(1)}(α) → 3α, so the leading term is
//!
//!   ϖ̇ = +(3/4) n / (1−e²)² · Σ_j (m_j/M_sun) (a_j/a)³.
//!
//! with mean motion n = √(GM_sun/a³). Same forcing bodies as the nodal rate in
//! `p9-2025-clustering::orbital_poles::nodal_precession_rate`, opposite sign
//! (nodes regress, apsides advance) and one higher power of α (the nodal
//! diagonal frequency keeps b_{3/2}^{(1)}/cos i; the apsidal one carries the
//! extra ᾱ_j from the eccentricity disturbing function).
//!
//! Because n ∝ a^{-3/2} and the interior-ring sum carries another a^{-2}, the
//! rate scales as ϖ̇ ∝ a^{-7/2}: higher-a objects precess much more slowly.
//! That *differential* precession is exactly what shears an initially aligned
//! set apart and is the mechanism the paper back-integrates. At sednoid
//! distances (e ≈ 0.85) the forced apsidal periods run from ~0.16 Gyr
//! (a ≈ 228 AU) to ~2.6 Gyr (a ≈ 506 AU): the closest detached objects wrap
//! many times over the solar system's age while Sedna completes a fraction
//! of a cycle — consistent with Huang & Gladman's figures.
//!
//! Giant-planet masses and semi-major axes are taken from the p9-core J2000
//! ephemeris states (no local planet table), matching the nodal-rate code.

use p9_core::analysis::elements::mean_motion;
use p9_core::constants::{GM_SUN, YEAR_DAYS};
use p9_core::initial_conditions::planets::giant_planets_j2000;
use p9_core::types::cartesian_to_elements;

/// The dimensionless secular forcing sum Σ_j (m_j/M_sun)(a_j/a)² from the four
/// giant planets (Jupiter…Neptune), evaluated from their J2000 ephemeris
/// elements. Decreases as a^{-2} with the test-particle semi-major axis.
///
/// For an INTERIOR perturber (planet inside the particle's orbit) the
/// Murray & Dermott coefficient is α_j·ᾱ_j·b ≈ (a_j/a)² at leading order
/// (ᾱ = 1 for internal perturbers) — equivalently the interior-ring
/// quadrupole with J2_eff = ½Σ(m_j/M)a_j². A previous version applied the
/// exterior-perturber convention (a_j/a)³, making every rate low by
/// a/a_j,eff ≈ 16–27× across the sample and the shear scale wrong.
pub fn giant_planet_forcing_sum(a: f64) -> f64 {
    giant_planets_j2000()
        .iter()
        .map(|body| {
            let a_j = cartesian_to_elements(&body.state, GM_SUN).a;
            body.mass * (a_j / a).powi(2)
        })
        .sum()
}

/// Apsidal precession rate dϖ/dt (rad/yr), prograde (positive), forced by the
/// giant planets on a particle of semi-major axis `a` (AU) and eccentricity
/// `e`. Inclination enters only weakly at this order and is carried through the
/// standard cos i factor (`i` in radians) shared with the nodal theory.
///
///   ϖ̇ = +(3/4) n cos(i)/(1−e²)² · Σ_j (m_j/M_sun)(a_j/a)²
pub fn apsidal_precession_rate(a: f64, e: f64, i: f64) -> f64 {
    let n = mean_motion(a); // rad/day
    let eta_sq = (1.0 - e * e) * (1.0 - e * e);
    let sum = giant_planet_forcing_sum(a);
    let rate_day = 0.75 * n * i.cos() / eta_sq * sum;
    rate_day * YEAR_DAYS // rad/yr
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn test_rate_is_prograde() {
        let r = apsidal_precession_rate(300.0, 0.7, 20.0 * DEG2RAD);
        assert!(r > 0.0, "apsidal precession should be prograde (advance)");
    }

    #[test]
    fn test_rate_monotonically_decreasing_in_a() {
        // The differential-precession engine: distant objects precess slower.
        let mut prev = f64::INFINITY;
        for a in [228.0, 261.0, 300.0, 350.0, 430.0, 506.0] {
            let r = apsidal_precession_rate(a, 0.85, 0.0);
            assert!(r < prev, "rate not decreasing at a = {a}: {r} >= {prev}");
            prev = r;
        }
    }

    #[test]
    fn test_rate_scales_as_a_to_minus_seven_halves() {
        // ϖ̇ ∝ n·a^{-2} ∝ a^{-7/2}: doubling a drops the rate by 2^{7/2} ≈ 11.3
        // (the lib-level a^{-7/2} claim, which the old a^{-9/2} code violated).
        let r1 = apsidal_precession_rate(250.0, 0.0, 0.0);
        let r2 = apsidal_precession_rate(500.0, 0.0, 0.0);
        let ratio = r1 / r2;
        assert!(
            (ratio - 2f64.powf(3.5)).abs() / 2f64.powf(3.5) < 1e-6,
            "ratio = {ratio}, expected {}",
            2f64.powf(3.5)
        );
    }

    #[test]
    fn test_rate_matches_interior_ring_quadrupole() {
        // At i = 0 the rate must equal the exact ring form
        // (3/2)·n·J2_eff/(a²η⁴) with J2_eff = ½Σ(m_j/M)a_j² (built from the
        // same p9-core giants the forcing sum uses).
        use p9_core::forces::j2_secular::{combined_j2_jsu, effective_j2};
        let (j2_jsu, _, _) = combined_j2_jsu();
        let a_nep = cartesian_to_elements(&giant_planets_j2000()[3].state, GM_SUN).a;
        let m_nep = giant_planets_j2000()[3].mass;
        let j2_eff = j2_jsu + effective_j2(m_nep, a_nep, 1.0);
        let (a, e) = (300.0, 0.5);
        let n = mean_motion(a);
        let eta4 = (1.0f64 - e * e).powi(2);
        let ring = 1.5 * n * j2_eff / (a * a * eta4) * YEAR_DAYS;
        let rate = apsidal_precession_rate(a, e, 0.0);
        // ~0.2% residual: the forcing sum uses J2000 osculating a_j from the
        // ephemeris states; combined_j2_jsu uses the mean-element constants.
        assert!(
            ((rate - ring) / ring).abs() < 5e-3,
            "rate {rate:.3e} vs ring {ring:.3e} rad/yr"
        );
    }

    #[test]
    fn test_sedna_like_period_is_gigayears() {
        // With the corrected interior-perturber scaling, a Sedna-like orbit
        // (a ≈ 506 AU, e = 0.85) has an apsidal period of a few Gyr — it has
        // precessed order-a-radian over the solar system's age, not the
        // ~10¹¹-yr near-frozen value the a^{-9/2} version produced.
        let rate = apsidal_precession_rate(506.0, 0.85, 0.0); // rad/yr
        let period_gyr = std::f64::consts::TAU / rate / 1.0e9;
        assert!(
            (1.0..40.0).contains(&period_gyr),
            "Sedna-like apsidal period = {period_gyr:.1} Gyr"
        );
    }

    #[test]
    fn test_forcing_sum_positive_and_decreasing() {
        let s1 = giant_planet_forcing_sum(250.0);
        let s2 = giant_planet_forcing_sum(500.0);
        assert!(s1 > 0.0 && s2 > 0.0);
        assert!(s2 < s1);
    }

    #[test]
    fn test_rate_magnitude_reasonable() {
        // A sednoid at a ≈ 500 AU precesses by ~ a few ×10⁻³ rad over a Gyr,
        // i.e. it has barely precessed; an a ≈ 230 AU object precesses an
        // order of magnitude faster. Sanity-bound the per-Gyr drift.
        let r_far = apsidal_precession_rate(500.0, 0.85, 0.0); // rad/yr
        let drift_far = r_far * 1.0e9; // rad over 1 Gyr
        assert!(
            drift_far > 0.0 && drift_far < 100.0,
            "implausible drift {drift_far} rad/Gyr"
        );
    }
}
