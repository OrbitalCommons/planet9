//! Circular (directional) statistics.
//!
//! Single source for the circular mean, mean resultant length, Rayleigh test
//! (with the small-n correction — the bare exp(−nR̄²) approximation is poor at
//! the n = 6–11 sample sizes of the clustering papers), and the Kuiper test.
//! Previously re-implemented inline at least six times across the workspace,
//! sometimes with arithmetic (non-circular) means.
//!
//! Reference: Mardia & Jupp (2000), "Directional Statistics"; Fisher (1993)

use crate::constants::TWO_PI;

/// Circular mean of a set of angles (radians). Returns a value in [0, 2π).
///
/// Returns `None` for an empty slice or when the resultant vector vanishes
/// (perfectly balanced angles have no defined mean direction).
pub fn circular_mean(angles: &[f64]) -> Option<f64> {
    if angles.is_empty() {
        return None;
    }
    let s: f64 = angles.iter().map(|a| a.sin()).sum();
    let c: f64 = angles.iter().map(|a| a.cos()).sum();
    if s.hypot(c) < 1e-12 * angles.len() as f64 {
        return None;
    }
    Some(s.atan2(c).rem_euclid(TWO_PI))
}

/// Mean resultant length R̄ ∈ [0, 1]: 1 = perfectly aligned, 0 = balanced.
pub fn mean_resultant_length(angles: &[f64]) -> f64 {
    if angles.is_empty() {
        return 0.0;
    }
    let n = angles.len() as f64;
    let s: f64 = angles.iter().map(|a| a.sin()).sum();
    let c: f64 = angles.iter().map(|a| a.cos()).sum();
    s.hypot(c) / n
}

/// Circular standard deviation: sqrt(−2 ln R̄) (radians).
pub fn circular_std(angles: &[f64]) -> f64 {
    let r = mean_resultant_length(angles);
    if r <= 0.0 {
        f64::INFINITY
    } else if r >= 1.0 {
        0.0
    } else {
        (-2.0 * r.ln()).sqrt()
    }
}

/// Rayleigh test for non-uniformity of circular data.
///
/// Returns the p-value for the null hypothesis that the angles are uniform on
/// the circle, using the series correction of Mardia & Jupp (2000), Eq. 6.3.5:
///
///   p ≈ exp(−Z)·[1 + (2Z − Z²)/(4n) − (24Z − 132Z² + 76Z³ − 9Z⁴)/(288n²)]
///
/// with Z = n R̄². The bare exp(−Z) first term is biased at small n (the
/// regime of the 6–11 object ETNO samples); the correction is accurate to
/// O(n⁻³).
pub fn rayleigh_p_value(angles: &[f64]) -> f64 {
    let n = angles.len() as f64;
    if n == 0.0 {
        return 1.0;
    }
    let r_bar = mean_resultant_length(angles);
    let z = n * r_bar * r_bar;

    let p = (-z).exp()
        * (1.0 + (2.0 * z - z * z) / (4.0 * n)
            - (24.0 * z - 132.0 * z * z + 76.0 * z.powi(3) - 9.0 * z.powi(4)) / (288.0 * n * n));
    p.clamp(0.0, 1.0)
}

/// Kuiper statistic V_n for circular uniformity.
///
/// V_n = D⁺ + D⁻ against the uniform CDF on [0, 2π); rotation-invariant,
/// unlike Kolmogorov-Smirnov.
pub fn kuiper_statistic(angles: &[f64]) -> f64 {
    let n = angles.len();
    if n == 0 {
        return 0.0;
    }
    let mut u: Vec<f64> = angles
        .iter()
        .map(|a| a.rem_euclid(TWO_PI) / TWO_PI)
        .collect();
    u.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let nf = n as f64;
    let mut d_plus: f64 = 0.0;
    let mut d_minus: f64 = 0.0;
    for (i, &ui) in u.iter().enumerate() {
        d_plus = d_plus.max((i as f64 + 1.0) / nf - ui);
        d_minus = d_minus.max(ui - i as f64 / nf);
    }
    d_plus + d_minus
}

/// Approximate p-value for the Kuiper statistic (Stephens 1970 asymptotic
/// series with the finite-n stabilization factor).
pub fn kuiper_p_value(angles: &[f64]) -> f64 {
    let n = angles.len() as f64;
    if n < 2.0 {
        return 1.0;
    }
    let v = kuiper_statistic(angles);
    // Stabilized statistic
    let lambda = (n.sqrt() + 0.155 + 0.24 / n.sqrt()) * v;
    if lambda < 0.4 {
        return 1.0;
    }
    let mut p = 0.0;
    for j in 1..=100 {
        let jf = j as f64;
        let a = 4.0 * jf * jf * lambda * lambda;
        let term = (2.0 * a - 2.0) * (-0.5 * a).exp();
        p += term;
        if term.abs() < 1e-12 {
            break;
        }
    }
    p.clamp(0.0, 1.0)
}

/// Wrap an angle difference into (−π, π].
pub fn wrap_to_pi(angle: f64) -> f64 {
    let mut a = angle.rem_euclid(TWO_PI);
    if a > std::f64::consts::PI {
        a -= TWO_PI;
    }
    a
}

/// Modified Bessel function I_0(x) approximation.
///
/// Uses polynomial approximation valid for all x >= 0.
pub fn bessel_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 3.75 {
        let t = (x / 3.75).powi(2);
        1.0 + t
            * (3.5156229
                + t * (3.0899424
                    + t * (1.2067492 + t * (0.2659732 + t * (0.0360768 + t * 0.0045813)))))
    } else {
        let t = 3.75 / ax;
        (ax.exp() / ax.sqrt())
            * (0.39894228
                + t * (0.01328592
                    + t * (0.00225319
                        + t * (-0.00157565
                            + t * (0.00916281
                                + t * (-0.02057706
                                    + t * (0.02635537 + t * (-0.01647633 + t * 0.00392377))))))))
    }
}

/// ln I₀(x): overflow-safe log of the modified Bessel function. `bessel_i0`
/// itself overflows to +inf for x ≳ 709 (exp(x) leaves f64), while von Mises
/// log-likelihood consumers only ever need the log — for large |x| this uses
/// the asymptotic ln I₀(x) = |x| − ln√(2π|x|) + ln(1 + 1/(8|x|) + …).
pub fn ln_bessel_i0(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 700.0 {
        bessel_i0(x).ln()
    } else {
        // Asymptotic series: I0(x) ~ e^x/sqrt(2 pi x) (1 + 1/(8x) + 9/(128x^2)).
        let corr = 1.0 + 1.0 / (8.0 * ax) + 9.0 / (128.0 * ax * ax);
        ax - 0.5 * (2.0 * std::f64::consts::PI * ax).ln() + corr.ln()
    }
}

/// Maximum-likelihood κ from the mean resultant length R̄
/// (Mardia & Jupp 2000 piecewise approximation of A⁻¹). Clamped at 500 —
/// far beyond any physical concentration in this workspace and safely below
/// the x ≈ 709 overflow of `bessel_i0`, so downstream `bessel_i0(κ).ln()`
/// stays finite (use `ln_bessel_i0` for arbitrary κ).
pub fn kappa_from_r_bar(r_bar: f64) -> f64 {
    let kappa = if r_bar < 0.53 {
        2.0 * r_bar + r_bar.powi(3) + 5.0 * r_bar.powi(5) / 6.0
    } else if r_bar < 0.85 {
        -0.4 + 1.39 * r_bar + 0.43 / (1.0 - r_bar)
    } else {
        let denom = r_bar.powi(3) - 4.0 * r_bar.powi(2) + 3.0 * r_bar;
        if denom.abs() < 1e-10 {
            500.0 // Very concentrated distribution
        } else {
            1.0 / denom
        }
    };
    kappa.clamp(0.0, 500.0)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_circular_mean_wraps_correctly() {
        // Angles straddling 0: arithmetic mean would give π, circular gives ~0.
        let angles = [0.1, TWO_PI - 0.1];
        let mean = circular_mean(&angles).unwrap();
        assert!(mean < 0.01 || mean > TWO_PI - 0.01, "mean = {mean}");
    }

    #[test]
    fn test_circular_mean_balanced_is_none() {
        let angles = [0.0, PI / 2.0, PI, 3.0 * PI / 2.0];
        assert!(circular_mean(&angles).is_none());
        assert!(circular_mean(&[]).is_none());
    }

    #[test]
    fn test_resultant_length_limits() {
        assert_relative_eq!(
            mean_resultant_length(&[1.0, 1.0, 1.0]),
            1.0,
            epsilon = 1e-12
        );
        assert!(mean_resultant_length(&[0.0, PI]) < 1e-12);
    }

    #[test]
    fn test_rayleigh_clustered_vs_uniform() {
        // Tightly clustered: small p.
        let clustered: Vec<f64> = (0..10).map(|k| 1.0 + 0.05 * k as f64).collect();
        assert!(rayleigh_p_value(&clustered) < 1e-3);

        // Evenly spread: p near 1.
        let uniform: Vec<f64> = (0..10).map(|k| k as f64 * TWO_PI / 10.0).collect();
        assert!(rayleigh_p_value(&uniform) > 0.9);
    }

    #[test]
    fn test_rayleigh_small_n_correction_matters() {
        // n = 6, R̄ moderate: the corrected p differs measurably from exp(-Z).
        let angles = [0.0, 0.4, 0.9, 1.5, 2.2, 3.0];
        let n = angles.len() as f64;
        let r = mean_resultant_length(&angles);
        let z = n * r * r;
        let naive = (-z).exp();
        let corrected = rayleigh_p_value(&angles);
        assert!(
            (naive - corrected).abs() / naive > 0.01,
            "correction too small: naive {naive:.4}, corrected {corrected:.4}"
        );
    }

    #[test]
    fn test_kuiper_uniform_vs_clustered() {
        let uniform: Vec<f64> = (0..20).map(|k| k as f64 * TWO_PI / 20.0).collect();
        let clustered: Vec<f64> = (0..20).map(|k| 1.0 + 0.02 * k as f64).collect();
        assert!(kuiper_statistic(&clustered) > kuiper_statistic(&uniform));
        assert!(kuiper_p_value(&clustered) < 0.01);
        assert!(kuiper_p_value(&uniform) > 0.5);
    }

    #[test]
    fn test_kuiper_rotation_invariant() {
        let angles: Vec<f64> = vec![0.2, 1.1, 2.0, 4.5, 5.9];
        let rotated: Vec<f64> = angles
            .iter()
            .map(|a| (a + 2.7).rem_euclid(TWO_PI))
            .collect();
        assert_relative_eq!(
            kuiper_statistic(&angles),
            kuiper_statistic(&rotated),
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_wrap_to_pi() {
        assert_relative_eq!(wrap_to_pi(3.0 * PI), PI, epsilon = 1e-12);
        assert_relative_eq!(wrap_to_pi(-0.5), -0.5, epsilon = 1e-12);
        assert_relative_eq!(wrap_to_pi(TWO_PI + 0.5), 0.5, epsilon = 1e-12);
    }

    #[test]
    fn test_bessel_i0_at_zero() {
        assert!((bessel_i0(0.0) - 1.0).abs() < 1e-6);
    }

    #[test]
    fn test_bessel_i0_positive() {
        for x in [0.5, 1.0, 2.0, 5.0, 10.0] {
            assert!(bessel_i0(x) > 0.0, "I_0({}) should be positive", x);
        }
    }

    #[test]
    fn test_kappa_from_r_bar_monotone() {
        assert!(kappa_from_r_bar(0.2) < kappa_from_r_bar(0.6));
        assert!(kappa_from_r_bar(0.6) < kappa_from_r_bar(0.9));
        assert!(kappa_from_r_bar(0.0) < 1e-9);
    }

    #[test]
    fn ln_bessel_stays_finite_where_bessel_overflows() {
        // Continuity across the asymptotic switch and finiteness far beyond
        // the exp overflow.
        let a = ln_bessel_i0(699.0);
        let b = ln_bessel_i0(701.0);
        assert!(
            (b - a - 2.0).abs() < 0.01,
            "slope ~1 near switch: {}",
            b - a
        );
        assert!(bessel_i0(800.0).is_infinite());
        let big = ln_bessel_i0(800.0);
        assert!(big.is_finite() && (big - 796.0).abs() < 1.0, "{big}");
        // Small-x agreement with the direct log.
        for &x in &[0.1, 1.0, 10.0, 100.0] {
            assert!((ln_bessel_i0(x) - bessel_i0(x).ln()).abs() < 1e-12);
        }
    }

    #[test]
    fn kappa_clamp_is_below_bessel_overflow() {
        let k = kappa_from_r_bar(0.999_999);
        assert!(k <= 500.0);
        assert!(bessel_i0(k).is_finite());
    }
}
