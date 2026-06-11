//! Mean-motion resonance machinery: resonant angles, libration detection,
//! and the Chirikov resonance-overlap (chaos) criterion for the Neptune
//! 2:j chain.
//!
//! Single sourced implementation: the workspace previously carried two
//! inconsistent resonant-angle conventions (one with p and q swapped, which
//! circulates even for an exactly resonant orbit) and three mutually
//! inconsistent Chirikov forms differing by exp(q²/4a_N²) ≈ 1.4 at
//! q = 35 AU.
//!
//! Conventions: a particle *exterior* to the planet in the p:q resonance
//! (p > q, both positive, p orbits of the planet per q orbits of the
//! particle) sits at a_res = a_planet (p/q)^{2/3} and has the stationary
//! angle
//!
//!   φ = p λ − q λ_planet − (p − q) ϖ
//!
//! which satisfies the d'Alembert rule (coefficients sum to zero) and is
//! stationary when n/n_planet = q/p.
//!
//! Reference: Murray & Dermott (1999) Ch. 8; Batygin & Brown (2021)

use crate::constants::{A_NEPTUNE_AU, MASS_NEPTUNE_SOLAR, TWO_PI};

/// Semi-major axis of the exterior p:q resonance with a planet at
/// `a_planet` (AU): a_res = a_planet (p/q)^{2/3}.
pub fn resonance_semi_major_axis(p: u32, q: u32, a_planet: f64) -> f64 {
    a_planet * (p as f64 / q as f64).powf(2.0 / 3.0)
}

/// Resonant angle φ = p λ − q λ_planet − (p − q) ϖ for the exterior p:q
/// resonance, wrapped to [0, 2π).
///
/// `lambda`/`varpi` are the particle's mean longitude and longitude of
/// perihelion; `lambda_planet` the planet's mean longitude (all radians).
pub fn resonant_angle(p: u32, q: u32, lambda: f64, lambda_planet: f64, varpi: f64) -> f64 {
    let (pf, qf) = (p as f64, q as f64);
    (pf * lambda - qf * lambda_planet - (pf - qf) * varpi).rem_euclid(TWO_PI)
}

/// Libration amplitude of a resonant-angle time series, or `None` if the
/// angle circulates.
///
/// Detection: sort the angles on the circle and find the largest angular gap.
/// A librating angle occupies an arc, leaving an empty gap > `min_gap`; a
/// circulating angle covers the circle (max gap → 2π·ln N/N for N uniform
/// samples). The returned amplitude is the half-width of the occupied arc.
///
/// `min_gap` of ~1 rad is a robust default for ≳ 50 samples spanning several
/// libration periods.
pub fn libration_amplitude(angles: &[f64], min_gap: f64) -> Option<f64> {
    if angles.len() < 2 {
        return None;
    }
    let mut sorted: Vec<f64> = angles.iter().map(|a| a.rem_euclid(TWO_PI)).collect();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut max_gap = TWO_PI - (sorted[sorted.len() - 1] - sorted[0]);
    for w in sorted.windows(2) {
        max_gap = max_gap.max(w[1] - w[0]);
    }

    if max_gap > min_gap {
        Some(0.5 * (TWO_PI - max_gap))
    } else {
        None
    }
}

/// Whether a resonant-angle time series librates (see `libration_amplitude`).
pub fn is_librating(angles: &[f64], min_gap: f64) -> bool {
    libration_amplitude(angles, min_gap).is_some()
}

/// Chirikov overlap parameter for the Neptune 2:j resonance chain:
///
///   K = Δa/δa = (24/√5) (a/a_N)^{5/4} √(m_N/M_sun) · exp(−q²/(2 a_N²))
///
/// Adjacent resonances overlap (K > 1) ⇒ chaotic transport. The exponent
/// carries q²/(2 a_N²) — this is the form consistent with the pendulum
/// half-width of the resonance potential and with `critical_perihelion`
/// below (a previous copy used q²/4a_N², a factor ≈ 1.4 at q = 35 AU).
pub fn chirikov_overlap_parameter(a_au: f64, q_au: f64) -> f64 {
    let alpha = a_au / A_NEPTUNE_AU;
    (24.0 / 5.0_f64.sqrt())
        * alpha.powf(1.25)
        * MASS_NEPTUNE_SOLAR.sqrt()
        * (-q_au * q_au / (2.0 * A_NEPTUNE_AU * A_NEPTUNE_AU)).exp()
}

/// Whether orbits at (a, q) are in the chaotic resonance-overlap regime.
pub fn is_chaotic(a_au: f64, q_au: f64) -> bool {
    chirikov_overlap_parameter(a_au, q_au) > 1.0
}

/// Critical perihelion below which the 2:j chain overlaps (the K = 1 root of
/// `chirikov_overlap_parameter`):
///
///   q_crit = a_N √( ln[ (576/5) (m_N/M_sun) (a/a_N)^{5/2} ] )
///
/// Returns 0 when the argument of the log is ≤ 1 (no chaos at any q).
pub fn critical_perihelion(a_au: f64) -> f64 {
    let alpha = a_au / A_NEPTUNE_AU;
    let argument = (576.0 / 5.0) * MASS_NEPTUNE_SOLAR * alpha.powf(2.5);
    if argument <= 1.0 {
        return 0.0;
    }
    A_NEPTUNE_AU * argument.ln().sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_resonant_angle_stationary_at_exact_resonance() {
        // A perfectly resonant orbit: n/n_N = q/p. The angle must be
        // constant in time — the p/q-swapped convention circulates at
        // n_N (p² − q²)/q and fails this test.
        let (p, q) = (5u32, 2u32);
        let n_planet = 1.0e-2; // rad/day
        let n = n_planet * (q as f64) / (p as f64);
        let varpi = 1.3;

        let phi0 = resonant_angle(p, q, varpi, 0.0, varpi);
        for step in 1..200 {
            let t = step as f64 * 50.0;
            let phi = resonant_angle(p, q, varpi + n * t, n_planet * t, varpi);
            assert_relative_eq!(phi, phi0, epsilon = 1e-9);
        }
    }

    #[test]
    fn test_resonant_angle_circulates_off_resonance() {
        let (p, q) = (3u32, 1u32);
        let n_planet = 1.0e-2;
        let n = n_planet * (q as f64) / (p as f64) * 1.05; // 5% off
        let angles: Vec<f64> = (0..500)
            .map(|s| {
                let t = s as f64 * 200.0;
                resonant_angle(p, q, n * t, n_planet * t, 0.0)
            })
            .collect();
        assert!(!is_librating(&angles, 1.0));
    }

    #[test]
    fn test_libration_detection() {
        // Librating: φ oscillates about π with amplitude 0.8.
        let librating: Vec<f64> = (0..200)
            .map(|k| std::f64::consts::PI + 0.8 * (0.1 * k as f64).sin())
            .collect();
        let amp = libration_amplitude(&librating, 1.0).expect("should librate");
        assert!((amp - 0.8).abs() < 0.05, "amplitude {amp}");

        // Circulating: uniform coverage.
        let circulating: Vec<f64> = (0..200).map(|k| k as f64 * TWO_PI / 200.0).collect();
        assert!(libration_amplitude(&circulating, 1.0).is_none());
    }

    #[test]
    fn test_resonance_location() {
        // Neptune 3:2 (Plutinos): a ≈ 39.4 AU
        let a = resonance_semi_major_axis(3, 2, A_NEPTUNE_AU);
        assert!((a - 39.4).abs() < 0.1, "a = {a}");
    }

    #[test]
    fn test_chirikov_critical_perihelion_consistent_with_overlap() {
        // The two public APIs must agree: K(a, q_crit(a)) = 1 exactly.
        for &a in &[300.0, 500.0, 1000.0] {
            let q_crit = critical_perihelion(a);
            assert!(q_crit > 0.0);
            let k = chirikov_overlap_parameter(a, q_crit);
            assert_relative_eq!(k, 1.0, epsilon = 1e-10);
            // And the classification flips across the boundary.
            assert!(is_chaotic(a, q_crit - 1.0));
            assert!(!is_chaotic(a, q_crit + 1.0));
        }
    }

    #[test]
    fn test_chirikov_trends() {
        assert!(chirikov_overlap_parameter(500.0, 35.0) > chirikov_overlap_parameter(100.0, 35.0));
        assert!(chirikov_overlap_parameter(200.0, 30.0) > chirikov_overlap_parameter(200.0, 45.0));
        // q_crit at a = 500 AU is in the trans-Neptunian range
        let q = critical_perihelion(500.0);
        assert!(q > 25.0 && q < 60.0, "q_crit = {q}");
    }
}
