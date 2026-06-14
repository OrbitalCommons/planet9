//! Conversion of the log-likelihood ratio to a Gaussian-equivalent
//! significance in σ.
//!
//! Under the isotropic null H0, Wilks' theorem states that 2Λ is
//! asymptotically distributed as χ² with k = 2 degrees of freedom, where k is
//! the number of free parameters of the von Mises model (µ and κ). For k = 2
//! the χ² survival function is closed form:
//!
//!   p = P(χ²₂ ≥ 2Λ) = exp(−Λ).
//!
//! We treat p as a two-sided Gaussian tail probability and map it to a
//! Gaussian-equivalent significance via the inverse normal CDF:
//!
//!   σ = Φ⁻¹(1 − p/2),
//!
//! so that p = 0.05 ↔ 1.96σ and p = 0.0027 ↔ 3σ. The normal quantile uses
//! the Acklam (2003) rational approximation, accurate to ~1e-9 relative.

/// χ²-with-2-dof survival function P(χ²₂ ≥ x) = exp(−x/2).
fn chi2_2dof_survival(x: f64) -> f64 {
    (-0.5 * x).exp()
}

/// Two-sided p-value for the apsidal-clustering log-likelihood ratio via
/// Wilks' theorem with 2 degrees of freedom: p = P(χ²₂ ≥ 2Λ) = exp(−Λ).
pub fn lambda_to_p_value(lambda: f64) -> f64 {
    let lambda = lambda.max(0.0);
    chi2_2dof_survival(2.0 * lambda).clamp(f64::MIN_POSITIVE, 1.0)
}

/// Inverse of the standard normal CDF, Φ⁻¹(p), via Acklam's (2003) rational
/// approximation. Valid for p ∈ (0, 1); relative error ≲ 1e-9.
pub fn normal_quantile(p: f64) -> f64 {
    assert!(p > 0.0 && p < 1.0, "normal_quantile domain is (0, 1)");

    // Coefficients for the central and tail rational approximations.
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_690e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239e0,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838e0,
        -2.549_732_539_343_734e0,
        4.374_664_141_464_968e0,
        2.938_163_982_698_783e0,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996e0,
        3.754_408_661_907_416e0,
    ];

    // Break-points of the central region.
    const P_LOW: f64 = 0.024_25;
    const P_HIGH: f64 = 1.0 - P_LOW;

    if p < P_LOW {
        // Lower tail.
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= P_HIGH {
        // Central region.
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        // Upper tail.
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

/// Convert a two-sided p-value to a Gaussian-equivalent significance in σ:
///
///   σ = Φ⁻¹(1 − p/2).
pub fn p_value_to_sigma(p: f64) -> f64 {
    let p = p.clamp(f64::MIN_POSITIVE, 1.0 - 1e-15);
    normal_quantile(1.0 - 0.5 * p)
}

/// Full pipeline: log-likelihood ratio Λ → Wilks p-value → Gaussian σ.
pub fn sigma(lambda: f64) -> f64 {
    p_value_to_sigma(lambda_to_p_value(lambda))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_known_points() {
        // Median.
        assert!((normal_quantile(0.5)).abs() < 1e-9);
        // 97.5th percentile ≈ 1.95996.
        assert!((normal_quantile(0.975) - 1.959_963_98).abs() < 1e-6);
        // 99.865th percentile ≈ 3.0 (one-sided 3σ).
        assert!((normal_quantile(0.998_650_1) - 3.0).abs() < 1e-3);
    }

    #[test]
    fn quantile_symmetry() {
        for p in [0.01, 0.1, 0.3, 0.45] {
            assert!((normal_quantile(p) + normal_quantile(1.0 - p)).abs() < 1e-8);
        }
    }

    #[test]
    fn p_to_sigma_known_values() {
        // p = 0.05 (two-sided) ↔ 1.96σ.
        let s = p_value_to_sigma(0.05);
        assert!((s - 1.959_963_98).abs() < 1e-4, "σ(0.05) = {s}");
        // p = 0.0027 (two-sided) ↔ 3σ.
        let s = p_value_to_sigma(0.0027);
        assert!((s - 3.0).abs() < 2e-3, "σ(0.0027) = {s}");
    }

    #[test]
    fn sigma_monotone_in_lambda() {
        assert!(sigma(1.0) < sigma(2.0));
        assert!(sigma(2.0) < sigma(5.0));
        // Λ = 0 → p = 1 → σ = 0.
        assert!(sigma(0.0) < 1e-6);
    }
}
