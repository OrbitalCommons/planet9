//! Shared scalar statistics: the standard-normal quantile and p-value ↔ σ
//! conversions used by the significance machinery across crates.

/// Standard-normal quantile Φ⁻¹(p), via Acklam's (2003) rational
/// approximation; |relative error| ≲ 1e-9 across the tails. Panics outside
/// the open interval (0, 1).
pub fn normal_quantile(p: f64) -> f64 {
    assert!(p > 0.0 && p < 1.0, "normal_quantile domain is (0, 1): {p}");

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
    const P_LOW: f64 = 0.024_25;

    if p < P_LOW {
        let q = (-2.0 * p.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p <= 1.0 - P_LOW {
        let q = p - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        -normal_quantile(1.0 - p)
    }
}

/// One-sided Gaussian-equivalent significance of a p-value:
/// σ = Φ⁻¹(1 − p). p ≥ 0.5 maps to non-positive σ; p is clamped away from
/// the exact endpoints so extreme MC p-values (0 or 1) stay finite.
pub fn p_value_to_sigma(p: f64) -> f64 {
    let p = p.clamp(1e-300, 1.0 - 1e-16);
    -normal_quantile(p)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantile_hits_reference_values() {
        // Φ⁻¹ at textbook points.
        assert!(normal_quantile(0.5).abs() < 1e-12);
        assert!((normal_quantile(0.975) - 1.959_963_984_540_054).abs() < 1e-8);
        assert!((normal_quantile(0.841_344_746_068_543) - 1.0).abs() < 1e-8);
        // Deep tail: Φ(−5) ≈ 2.866516e-7.
        assert!((normal_quantile(2.866_515_719e-7) + 5.0).abs() < 1e-6);
        // Antisymmetry.
        for &p in &[0.001, 0.1, 0.3] {
            assert!((normal_quantile(p) + normal_quantile(1.0 - p)).abs() < 1e-9);
        }
    }

    #[test]
    fn p_to_sigma_matches_convention() {
        assert!((p_value_to_sigma(0.0034) - 2.706).abs() < 1e-2);
        assert!((p_value_to_sigma(2.866_515_719e-7) - 5.0).abs() < 1e-5);
        assert!(p_value_to_sigma(0.5).abs() < 1e-12);
        // Extreme inputs stay finite.
        assert!(p_value_to_sigma(0.0).is_finite());
        assert!(p_value_to_sigma(1.0).is_finite());
    }
}
