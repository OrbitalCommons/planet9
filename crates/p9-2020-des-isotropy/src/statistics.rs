//! A battery of circular-isotropy test statistics.
//!
//! Bernardinelli et al. (2020) apply Kuiper's test to three orbital angles in
//! four samples (12 tests). This module reproduces a *suite* of isotropy
//! statistics as a stand-in for the paper's battery:
//!
//!   * Rayleigh test — first-harmonic non-uniformity (p9_core primitive).
//!   * Kuiper test V_n — the paper's actual statistic (p9_core primitive).
//!   * Mean resultant length R̄ — clustering magnitude (p9_core primitive).
//!   * Watson U² — a quadratic EDF statistic, complementary to Kuiper.
//!   * Beran / Ajne A_n — a smooth-alternative non-parametric statistic.
//!   * Two-angle circular correlation (Fisher–Lee) — tests joint structure
//!     between two angles (e.g. ω vs Ω).
//!
//! Analytic p-values are used where standard ones exist (Rayleigh, Kuiper);
//! Watson U², Beran and the correlation use a seeded Monte Carlo null so the
//! small-n (n = 3–7) calibration is honest rather than asymptotic.

use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use p9_core::analysis::circular::{
    kuiper_p_value, kuiper_statistic, mean_resultant_length, rayleigh_p_value, wrap_to_pi,
};
use p9_core::constants::TWO_PI;

/// Watson U² statistic against the uniform distribution on the circle.
///
/// U² = Σ (u_i − ū − (2i−1)/(2n) + 1/2)²·... — using the standard EDF form
/// (Watson 1961). Rotation-invariant, like Kuiper, but weights the whole
/// distribution rather than its extreme deviation.
pub fn watson_u2(angles: &[f64]) -> f64 {
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
    let u_bar: f64 = u.iter().sum::<f64>() / nf;
    let mut sum_sq = 0.0;
    for (i, &ui) in u.iter().enumerate() {
        let ci = (2.0 * (i as f64 + 1.0) - 1.0) / (2.0 * nf);
        let d = ui - ci;
        sum_sq += d * d;
    }
    sum_sq + 1.0 / (12.0 * nf) - nf * (u_bar - 0.5).powi(2)
}

/// Beran / Ajne A_n statistic for circular uniformity.
///
/// A_n = n/4 − (1/(nπ)) Σ_{i<j} d_ij where d_ij ∈ [0, π] is the angular
/// separation between θ_i and θ_j (Ajne 1968, a special Beran case with a
/// step weight). Under uniformity E[d_ij] = π/2 so A_n ≈ 1/4; as the angles
/// cluster the separations shrink toward 0 and A_n grows toward n/4.
/// Sensitive to a single concentrated mode.
pub fn beran_an(angles: &[f64]) -> f64 {
    let n = angles.len();
    if n < 2 {
        return 0.0;
    }
    let nf = n as f64;
    let pi = std::f64::consts::PI;
    let mut sum = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            // Angular separation in [0, π].
            let sep = wrap_to_pi(angles[i] - angles[j]).abs();
            sum += sep;
        }
    }
    nf / 4.0 - sum / (nf * pi)
}

/// Fisher–Lee circular–circular correlation coefficient between two angle
/// samples (paired). r_T ∈ [−1, 1]; |r_T| near 0 → no association.
pub fn circular_correlation(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len();
    assert_eq!(n, b.len(), "paired samples must be equal length");
    if n < 2 {
        return 0.0;
    }
    let mut num = 0.0;
    let mut den_a = 0.0;
    let mut den_b = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let sa = (a[i] - a[j]).sin();
            let sb = (b[i] - b[j]).sin();
            num += sa * sb;
            den_a += sa * sa;
            den_b += sb * sb;
        }
    }
    let den = (den_a * den_b).sqrt();
    if den < 1e-15 {
        0.0
    } else {
        num / den
    }
}

/// Draw n uniform angles on [0, 2π) using `rng`.
fn uniform_angles(rng: &mut StdRng, n: usize) -> Vec<f64> {
    (0..n).map(|_| rng.gen::<f64>() * TWO_PI).collect()
}

/// Monte Carlo upper-tail p-value: fraction of uniform-null draws whose
/// statistic is >= the observed value. `seed` fixes the null for
/// reproducibility; `iters` controls calibration.
pub fn mc_p_value_upper<F>(observed: f64, n: usize, seed: u64, iters: usize, stat: F) -> f64
where
    F: Fn(&[f64]) -> f64,
{
    let mut rng = StdRng::seed_from_u64(seed);
    let mut ge = 0usize;
    for _ in 0..iters {
        let sample = uniform_angles(&mut rng, n);
        if stat(&sample) >= observed {
            ge += 1;
        }
    }
    // (ge + 1) / (iters + 1): never returns exactly 0, standard MC convention.
    (ge as f64 + 1.0) / (iters as f64 + 1.0)
}

/// Monte Carlo two-sided p-value for the circular correlation under the null
/// of independence (one angle's pairing is randomly permuted). Uses |r_T|.
pub fn mc_correlation_p_value(a: &[f64], b: &[f64], seed: u64, iters: usize) -> f64 {
    let observed = circular_correlation(a, b).abs();
    let n = a.len();
    let mut rng = StdRng::seed_from_u64(seed);
    let mut ge = 0usize;
    let mut bp = b.to_vec();
    for _ in 0..iters {
        // Fisher–Yates shuffle of b's pairing.
        for k in (1..n).rev() {
            let j = rng.gen_range(0..=k);
            bp.swap(k, j);
        }
        if circular_correlation(a, &bp).abs() >= observed {
            ge += 1;
        }
    }
    (ge as f64 + 1.0) / (iters as f64 + 1.0)
}

/// The full battery of single-angle isotropy statistics with p-values.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct IsotropySuite {
    pub n: usize,
    pub r_bar: f64,
    pub rayleigh_p: f64,
    pub kuiper_v: f64,
    pub kuiper_p: f64,
    pub watson_u2: f64,
    pub watson_p: f64,
    pub beran_an: f64,
    pub beran_p: f64,
}

impl IsotropySuite {
    /// Compute all single-angle statistics for a set of angles (radians).
    ///
    /// The Watson U² and Beran p-values use a seeded uniform Monte Carlo null
    /// (the small-n regime makes asymptotic tables unreliable).
    pub fn compute(angles: &[f64], seed: u64, mc_iters: usize) -> Self {
        let n = angles.len();
        let r_bar = mean_resultant_length(angles);
        let rayleigh_p = rayleigh_p_value(angles);
        let kuiper_v = kuiper_statistic(angles);
        let kuiper_p = kuiper_p_value(angles);
        let u2 = watson_u2(angles);
        let an = beran_an(angles);
        let watson_p = mc_p_value_upper(u2, n, seed, mc_iters, watson_u2);
        let beran_p = mc_p_value_upper(an, n, seed.wrapping_add(1), mc_iters, beran_an);
        IsotropySuite {
            n,
            r_bar,
            rayleigh_p,
            kuiper_v,
            kuiper_p,
            watson_u2: u2,
            watson_p,
            beran_an: an,
            beran_p,
        }
    }

    /// The p-values of the four hypothesis tests in the suite.
    pub fn p_values(&self) -> [f64; 4] {
        [self.rayleigh_p, self.kuiper_p, self.watson_p, self.beran_p]
    }

    /// Number of tests in the suite returning p > alpha (consistent with
    /// isotropy at level alpha).
    pub fn n_consistent(&self, alpha: f64) -> usize {
        self.p_values().iter().filter(|&&p| p > alpha).count()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn clustered(n: usize, spread: f64) -> Vec<f64> {
        (0..n)
            .map(|k| 1.0 + spread * (k as f64 - n as f64 / 2.0))
            .collect()
    }
    fn even(n: usize) -> Vec<f64> {
        (0..n).map(|k| k as f64 * TWO_PI / n as f64).collect()
    }

    #[test]
    fn test_watson_u2_clustered_gt_uniform() {
        assert!(watson_u2(&clustered(12, 0.02)) > watson_u2(&even(12)));
    }

    #[test]
    fn test_beran_clustered_gt_uniform() {
        assert!(beran_an(&clustered(12, 0.02)) > beran_an(&even(12)));
    }

    #[test]
    fn test_correlation_self_is_one() {
        let a = vec![0.2, 1.1, 2.3, 4.0, 5.5];
        assert!((circular_correlation(&a, &a) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn test_correlation_independent_near_zero() {
        let a = vec![0.1, 1.0, 2.0, 3.0, 4.0, 5.0];
        let b = vec![3.0, 0.2, 5.0, 1.5, 4.4, 2.1];
        assert!(circular_correlation(&a, &b).abs() < 0.6);
    }

    #[test]
    fn test_mc_uniform_p_is_calibrated() {
        // A genuinely uniform sample should not look significant on average.
        let p = mc_p_value_upper(watson_u2(&even(10)), 10, 42, 4000, watson_u2);
        assert!(p > 0.2, "uniform Watson p = {p}");
    }

    #[test]
    fn test_mc_clustered_p_is_small() {
        let c = clustered(10, 0.03);
        let p = mc_p_value_upper(watson_u2(&c), 10, 42, 4000, watson_u2);
        assert!(p < 0.05, "clustered Watson p = {p}");
    }
}
