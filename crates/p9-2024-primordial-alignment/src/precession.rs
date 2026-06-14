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
//! Because n ∝ a^{-3/2} and the sum carries another a^{-3}, the rate scales as
//! ϖ̇ ∝ a^{-9/2}: higher-a objects precess much more slowly. That *differential*
//! precession is exactly what shears an initially aligned set apart and is the
//! mechanism the paper back-integrates. At sednoid distances the forced
//! precession periods run from ~2.6×10¹⁰ yr (a ≈ 228 AU) to ~9×10¹¹ yr
//! (a ≈ 506 AU), so over the 4.5 Gyr age of the solar system the closest
//! detached objects swing through ~1 rad while Sedna barely moves.
//!
//! Giant-planet masses and semi-major axes are taken from the p9-core J2000
//! ephemeris states (no local planet table), matching the nodal-rate code.

use p9_core::constants::{GM_SUN, YEAR_DAYS};
use p9_core::initial_conditions::planets::giant_planets_j2000;
use p9_core::types::cartesian_to_elements;

/// Mean motion n = √(GM_sun/a³) in rad/day.
pub fn mean_motion(a: f64) -> f64 {
    (GM_SUN / (a * a * a)).sqrt()
}

/// The dimensionless secular forcing sum Σ_j (m_j/M_sun)(a_j/a)³ from the four
/// giant planets (Jupiter…Neptune), evaluated from their J2000 ephemeris
/// elements. Decreases as a^{-3} with the test-particle semi-major axis.
pub fn giant_planet_forcing_sum(a: f64) -> f64 {
    giant_planets_j2000()
        .iter()
        .map(|body| {
            let a_j = cartesian_to_elements(&body.state, GM_SUN).a;
            body.mass * (a_j / a).powi(3)
        })
        .sum()
}

/// Apsidal precession rate dϖ/dt (rad/yr), prograde (positive), forced by the
/// giant planets on a particle of semi-major axis `a` (AU) and eccentricity
/// `e`. Inclination enters only weakly at this order and is carried through the
/// standard cos i factor (`i` in radians) shared with the nodal theory.
///
///   ϖ̇ = +(3/4) n cos(i)/(1−e²)² · Σ_j (m_j/M_sun)(a_j/a)³
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
    fn test_rate_scales_as_a_to_minus_nine_halves() {
        // ϖ̇ ∝ a^{-9/2}: doubling a should drop the rate by 2^{9/2} ≈ 22.6.
        let r1 = apsidal_precession_rate(250.0, 0.0, 0.0);
        let r2 = apsidal_precession_rate(500.0, 0.0, 0.0);
        let ratio = r1 / r2;
        assert!(
            (ratio - 2f64.powf(4.5)).abs() / 2f64.powf(4.5) < 1e-6,
            "ratio = {ratio}, expected {}",
            2f64.powf(4.5)
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
