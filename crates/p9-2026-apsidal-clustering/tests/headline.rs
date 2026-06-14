//! Headline reproduction of Siraj, Chyba & Tremaine (2026), "Measuring
//! Apsidal Clustering" (arXiv:2604.25990).
//!
//! The conditional-likelihood estimator applied to the seeded stable samples
//! reproduces the paper's central qualitative claim: the tightly clustered
//! n = 21 stable-ETNO sample yields ≈2.7σ, and adding 4 less-aligned stable
//! objects (n = 25) DILUTES the signal to ≈1.9σ — the SAME estimator on the
//! larger sample returns a lower significance.
//!
//! Computed values (seed 2026, this build):
//!   n = 21: κ ≈ 1.07, Λ ≈ 4.99, σ ≈ 2.71  (published 2.7σ; residual ≈ 0.01σ)
//!   n = 25: κ ≈ 0.70, Λ ≈ 2.81, σ ≈ 1.88  (published 1.9σ; residual ≈ 0.02σ)
//! The residuals are well inside the pinned bands; the von Mises fit and Λ are
//! recomputed from the samples, not hard-coded.

use p9_2026_apsidal_clustering::{
    estimator::fit, sigma, stable_21, stable_25, PUBLISHED_SIGMA_21, PUBLISHED_SIGMA_25,
};

#[test]
fn stable_21_significance_is_about_2_7_sigma() {
    let s = stable_21();
    let f = fit(&s);
    let sig = sigma(f.lambda);
    assert_eq!(f.n, 21);
    assert!(
        (2.3..=3.1).contains(&sig),
        "σ(n=21) = {sig:.3}, want [2.3, 3.1] (published {PUBLISHED_SIGMA_21})"
    );
}

#[test]
fn stable_25_significance_is_about_1_9_sigma() {
    let s = stable_25();
    let f = fit(&s);
    let sig = sigma(f.lambda);
    assert_eq!(f.n, 25);
    assert!(
        (1.5..=2.3).contains(&sig),
        "σ(n=25) = {sig:.3}, want [1.5, 2.3] (published {PUBLISHED_SIGMA_25})"
    );
}

#[test]
fn adding_less_aligned_objects_reduces_significance() {
    let sig21 = sigma(fit(&stable_21()).lambda);
    let sig25 = sigma(fit(&stable_25()).lambda);
    assert!(
        sig25 < sig21,
        "diluting the sample should strictly reduce σ: σ(25) = {sig25:.3} !< σ(21) = {sig21:.3}"
    );
}

#[test]
fn dilution_drops_concentration() {
    // The von Mises concentration κ also falls when the less-aligned objects
    // are added — the mechanism behind the significance drop.
    let k21 = fit(&stable_21()).kappa;
    let k25 = fit(&stable_25()).kappa;
    assert!(k25 < k21, "κ(25) = {k25:.3} !< κ(21) = {k21:.3}");
}
