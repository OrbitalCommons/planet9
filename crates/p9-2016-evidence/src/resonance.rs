//! Mean-motion resonance detection.
//!
//! Identifies particles trapped in mean-motion resonances with Planet Nine.
//! The paper identifies resonances: 2:1, 3:1, 5:3, 7:4, 9:4, 11:4, 13:4,
//! 23:6, 27:17, 29:17, 33:19.
//!
//! A particle is in p:q resonance when:
//!   a_particle / a_p9 ≈ (q/p)^(2/3)
//!
//! More precisely, the resonant angle φ = p*λ_particle - q*λ_p9 - (p-q)*ϖ_particle
//! should librate (oscillate around a fixed value) rather than circulate.

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

/// Compute the expected semi-major axis for a p:q resonance with a perturber
/// at semi-major axis a_p.
///
/// a_res = a_p * (q/p)^(2/3)
pub fn resonance_semimajor_axis(a_p: f64, p: u32, q: u32) -> f64 {
    a_p * (q as f64 / p as f64).powf(2.0 / 3.0)
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

/// Compute the resonant angle for a p:q resonance.
///
/// φ = p*λ_particle - q*λ_p9 - (p-q)*ϖ_particle
///
/// where λ = M + ω + Ω (mean longitude) and ϖ = ω + Ω (longitude of perihelion).
pub fn resonant_angle(
    elem_particle: &OrbitalElements,
    elem_p9: &OrbitalElements,
    p: u32,
    q: u32,
) -> f64 {
    let lambda_part = elem_particle.mean_anomaly + elem_particle.omega + elem_particle.omega_big;
    let lambda_p9 = elem_p9.mean_anomaly + elem_p9.omega + elem_p9.omega_big;
    let varpi_part = elem_particle.omega + elem_particle.omega_big;

    let mut phi =
        p as f64 * lambda_part - q as f64 * lambda_p9 - (p as i64 - q as i64) as f64 * varpi_part;

    // Wrap to [-π, π]
    phi = phi % TWO_PI;
    if phi > std::f64::consts::PI {
        phi -= TWO_PI;
    }
    if phi < -std::f64::consts::PI {
        phi += TWO_PI;
    }

    phi
}

/// Check if a series of resonant angles indicates libration (resonance trapping).
///
/// A librating angle stays within a bounded range (< 2π).
/// A circulating angle covers the full 2π range.
///
/// Returns (is_librating, amplitude) where amplitude is the peak-to-peak range.
pub fn is_librating(angles: &[f64]) -> (bool, f64) {
    if angles.len() < 10 {
        return (false, TWO_PI);
    }

    // Compute the range of the angle, accounting for wraparound
    // Use the circular statistics approach
    let sin_sum: f64 = angles.iter().map(|&a| a.sin()).sum();
    let cos_sum: f64 = angles.iter().map(|&a| a.cos()).sum();
    let n = angles.len() as f64;

    let r_bar = (sin_sum * sin_sum + cos_sum * cos_sum).sqrt() / n;

    // r_bar close to 1 means tightly clustered (librating)
    // r_bar close to 0 means spread out (circulating)
    let is_lib = r_bar > 0.5;

    // Estimate amplitude from circular dispersion
    let amplitude = if r_bar > 0.0 {
        2.0 * (-2.0 * r_bar.ln()).sqrt()
    } else {
        TWO_PI
    };

    (is_lib, amplitude.min(TWO_PI))
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
