//! Incompatibility metric between the Siraj 2024 and Brown & Batygin 2021
//! perturber solutions.
//!
//! Two independent inferences disagree if their inferred quantity differs by
//! more than the combined uncertainty. We quote the "tension in σ" for an
//! inferred quantity x with central values x₁, x₂ and 1σ scatters σ₁, σ₂:
//!
//!   T = |x₁ − x₂| / √(σ₁² + σ₂²).
//!
//! Siraj et al. (2024) (≈4.4 M⊕, a ≈ 290 AU, i ≈ 6.8°) sit well away from
//! Brown & Batygin (2021) (≈6.2 M⊕, a ≈ 380 AU, i ≈ 16°). A tension above ~2σ
//! marks the solutions as incompatible.

/// Tension (in units of the combined 1σ) between two independent estimates of
/// the same scalar quantity.
pub fn tension_sigma(x1: f64, sigma1: f64, x2: f64, sigma2: f64) -> f64 {
    let comb = (sigma1 * sigma1 + sigma2 * sigma2).sqrt();
    if comb <= 0.0 {
        return f64::INFINITY;
    }
    (x1 - x2).abs() / comb
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_tension_basic() {
        // 10 ± 3 vs 4 ± 4: combined σ = 5, separation 6 ⇒ 1.2σ.
        assert_relative_eq!(tension_sigma(10.0, 3.0, 4.0, 4.0), 1.2, epsilon = 1e-12);
    }

    #[test]
    fn test_tension_zero_scatter_is_infinite() {
        assert!(tension_sigma(5.0, 0.0, 4.0, 0.0).is_infinite());
        assert_eq!(tension_sigma(5.0, 0.0, 5.0, 0.0), f64::INFINITY);
    }
}
