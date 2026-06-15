//! Mean-motion resonance detection.
//!
//! Identifies particles trapped in mean-motion resonances with Planet Nine.
//! The paper identifies resonances: 2:1, 3:1, 5:3, 7:4, 9:4, 11:4, 13:4,
//! 23:6, 27:17, 29:17, 33:19.
//!
//! Conventions: the p:q labels here are particle:P9 orbit counts, i.e. the
//! particle is *interior* to Planet Nine (n/n₉ = p/q with p > q) and sits at
//! a_res = a₉ (q/p)^{2/3}. The stationary resonant angle is therefore
//!
//!   φ = q λ − p λ₉ + (p − q) ϖ
//!
//! which is p9-core's exterior convention with the roles of particle and
//! planet exchanged. (A previous local implementation had p and q swapped,
//! so the angle circulated at n₉(p² − q²)/q even for an exactly resonant
//! orbit and `is_librating` could never fire.)

use p9_core::analysis::circular::wrap_to_pi;
use p9_core::analysis::resonance as core_resonance;
use p9_core::constants::*;
use p9_core::types::OrbitalElements;

/// Known mean-motion resonances from Batygin & Brown (2016).
/// Each tuple is (p, q) where the resonance is p:q (particle:P9).
pub const KNOWN_RESONANCES: &[(u32, u32)] = &[
    (2, 1),
    (3, 1),
    (5, 3),
    (7, 4),
    (9, 4),
    (11, 4),
    (13, 4),
    (23, 6),
    (27, 17),
    (29, 17),
    (33, 19),
];

/// Minimum empty arc (radians) used by the gap-based libration detector.
/// ~1 rad is robust for ≳50 samples spanning several libration periods.
pub const LIBRATION_MIN_GAP: f64 = 1.0;

/// Compute the expected semi-major axis for a p:q (particle:P9) resonance
/// with a perturber at semi-major axis `a_p`: a_res = a_p (q/p)^{2/3}.
///
/// This is p9-core's exterior `resonance_semi_major_axis` with the particle
/// and planet roles exchanged (the particle here is interior to P9).
pub fn resonance_semimajor_axis(a_p: f64, p: u32, q: u32) -> f64 {
    core_resonance::resonance_semi_major_axis(q, p, a_p)
}

/// Identify which known resonance (if any) a particle with semi-major axis `a`
/// might be near, given perturber semi-major axis `a_p`.
///
/// Returns Some((p, q, delta_a)) if within tolerance, None otherwise.
pub fn identify_resonance(a: f64, a_p: f64, tolerance_frac: f64) -> Option<(u32, u32, f64)> {
    let mut best: Option<(u32, u32, f64)> = None;

    for &(p, q) in KNOWN_RESONANCES {
        let a_res = resonance_semimajor_axis(a_p, p, q);
        let delta = (a - a_res).abs() / a_res;

        if delta < tolerance_frac {
            if best.is_none() || delta < best.unwrap().2 {
                best = Some((p, q, delta));
            }
        }
    }

    best
}

/// Compute the resonant angle for a p:q (particle:P9, interior) resonance:
///
///   φ = q λ − p λ₉ + (p − q) ϖ
///
/// where λ = M + ω + Ω and ϖ = ω + Ω. Satisfies the d'Alembert rule and is
/// stationary when n/n₉ = p/q. Returned wrapped to (−π, π].
pub fn resonant_angle(
    elem_particle: &OrbitalElements,
    elem_p9: &OrbitalElements,
    p: u32,
    q: u32,
) -> f64 {
    let lambda_part = elem_particle.mean_anomaly + elem_particle.omega + elem_particle.omega_big;
    let lambda_p9 = elem_p9.mean_anomaly + elem_p9.omega + elem_p9.omega_big;
    let varpi_part = elem_particle.omega + elem_particle.omega_big;

    // Interior particle: exchange roles in p9-core's exterior convention,
    // i.e. core's (p, q, λ, λ_planet, ϖ) ← (q, p, λ, λ₉, ϖ).
    let phi = core_resonance::resonant_angle(q, p, lambda_part, lambda_p9, varpi_part);
    wrap_to_pi(phi)
}

/// Check if a series of resonant angles indicates libration (resonance
/// trapping), using p9-core's largest-empty-gap detector.
///
/// Returns (is_librating, amplitude) where amplitude is the libration
/// half-width (2π for a circulating angle).
pub fn is_librating(angles: &[f64]) -> (bool, f64) {
    match core_resonance::libration_amplitude(angles, LIBRATION_MIN_GAP) {
        Some(amp) => (true, amp),
        None => (false, TWO_PI),
    }
}

/// Scan all known resonances and report which ones have particles.
pub fn resonance_census(
    elements: &[OrbitalElements],
    a_p: f64,
    tolerance: f64,
) -> Vec<(u32, u32, usize, f64)> {
    let mut census = Vec::new();

    for &(p, q) in KNOWN_RESONANCES {
        let a_res = resonance_semimajor_axis(a_p, p, q);
        let count = elements
            .iter()
            .filter(|e| ((e.a - a_res) / a_res).abs() < tolerance)
            .count();

        if count > 0 {
            census.push((p, q, count, a_res));
        }
    }

    census
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build particle/P9 element pairs along an exactly 2:1 resonant
    /// trajectory (n = 2 n₉), with an optional libration oscillation in λ.
    fn two_to_one_series(
        lib_amp: f64,
        n_samples: usize,
    ) -> Vec<(OrbitalElements, OrbitalElements)> {
        let a_p9 = 700.0;
        let (p, q) = (2u32, 1u32);
        let a_res = resonance_semimajor_axis(a_p9, p, q);

        let n_p9 = TWO_PI / (a_p9.powf(1.5) * 365.25); // rad/day (Kepler, a in AU)
        let n_part = n_p9 * p as f64 / q as f64;
        let varpi = 1.3;

        (0..n_samples)
            .map(|k| {
                let t = k as f64 * 2.0e4; // days
                                          // Libration enters the angle through q·δλ; oscillate λ so the
                                          // resonant angle swings by ±lib_amp.
                let dlambda = (lib_amp / q as f64) * (1.7e-6 * t).sin();
                let elem = OrbitalElements {
                    a: a_res,
                    e: 0.5,
                    i: 0.0,
                    omega: varpi,
                    omega_big: 0.0,
                    mean_anomaly: n_part * t + dlambda,
                };
                let elem_p9 = OrbitalElements {
                    a: a_p9,
                    e: 0.6,
                    i: 0.0,
                    omega: 0.4,
                    omega_big: 0.0,
                    mean_anomaly: n_p9 * t,
                };
                (elem, elem_p9)
            })
            .collect()
    }

    #[test]
    fn test_resonance_semimajor_axes() {
        let a_p = 700.0; // P9 nominal

        // 2:1 resonance: a_res = 700 * (1/2)^(2/3) ≈ 441 AU
        let a_21 = resonance_semimajor_axis(a_p, 2, 1);
        assert!((a_21 - 441.0).abs() < 1.0, "2:1 at {:.1} AU", a_21);

        // 3:1 resonance: a_res = 700 * (1/3)^(2/3) ≈ 336 AU
        let a_31 = resonance_semimajor_axis(a_p, 3, 1);
        assert!((a_31 - 336.0).abs() < 1.0, "3:1 at {:.1} AU", a_31);
    }

    #[test]
    fn test_identify_resonance() {
        let a_p = 700.0;

        // Particle near the 2:1 resonance
        let result = identify_resonance(440.0, a_p, 0.02);
        assert!(result.is_some());
        let (p, q, _) = result.unwrap();
        assert_eq!((p, q), (2, 1));

        // Particle far from any resonance
        let result = identify_resonance(600.0, a_p, 0.01);
        assert!(result.is_none());
    }

    #[test]
    fn test_exact_two_to_one_trajectory_librates() {
        // A constructed 2:1 resonant trajectory with a ±0.5 rad libration:
        // the correct angle librates; the previous p↔q-swapped convention
        // circulated at n₉(p² − q²)/q and classified this as circulating.
        let series = two_to_one_series(0.5, 400);
        let angles: Vec<f64> = series
            .iter()
            .map(|(e, e9)| resonant_angle(e, e9, 2, 1))
            .collect();

        let (is_lib, amp) = is_librating(&angles);
        assert!(is_lib, "exactly resonant trajectory must librate");
        assert!((amp - 0.5).abs() < 0.1, "amplitude {amp} (expected ~0.5)");
    }

    #[test]
    fn test_exact_resonance_angle_stationary() {
        // With zero libration amplitude the resonant angle is constant.
        let series = two_to_one_series(0.0, 100);
        let phi0 = resonant_angle(&series[0].0, &series[0].1, 2, 1);
        for (e, e9) in &series {
            let phi = resonant_angle(e, e9, 2, 1);
            assert!((phi - phi0).abs() < 1e-8, "angle drifted: {phi} vs {phi0}");
        }
    }

    #[test]
    fn test_off_resonance_circulates() {
        // 10% off the 2:1 commensurability: the angle must circulate
        // (φ̇ = 0.2 n₉, so the series below spans ~2.4 circulations).
        let a_p9: f64 = 700.0;
        let n_p9 = TWO_PI / (a_p9.powf(1.5) * 365.25);
        let n_part = n_p9 * 2.0 * 1.10;

        let angles: Vec<f64> = (0..400)
            .map(|k| {
                let t = k as f64 * 2.0e5;
                let elem = OrbitalElements {
                    a: 420.0,
                    e: 0.5,
                    i: 0.0,
                    omega: 1.3,
                    omega_big: 0.0,
                    mean_anomaly: n_part * t,
                };
                let elem_p9 = OrbitalElements {
                    a: a_p9,
                    e: 0.6,
                    i: 0.0,
                    omega: 0.4,
                    omega_big: 0.0,
                    mean_anomaly: n_p9 * t,
                };
                resonant_angle(&elem, &elem_p9, 2, 1)
            })
            .collect();

        let (is_lib, _) = is_librating(&angles);
        assert!(!is_lib, "off-resonance trajectory must circulate");
    }

    #[test]
    fn test_libration_detection() {
        // Librating angles clustered near 0
        let librating: Vec<f64> = (0..100).map(|i| 0.3 * (i as f64 * 0.1).sin()).collect();
        let (is_lib, amp) = is_librating(&librating);
        assert!(is_lib, "Should detect libration");
        assert!(amp < std::f64::consts::PI, "Amplitude should be < π");

        // Circulating angles covering full range
        let circulating: Vec<f64> = (0..100)
            .map(|i| (i as f64 / 100.0 * TWO_PI) - std::f64::consts::PI)
            .collect();
        let (is_lib, _) = is_librating(&circulating);
        assert!(!is_lib, "Should detect circulation");
    }
}
