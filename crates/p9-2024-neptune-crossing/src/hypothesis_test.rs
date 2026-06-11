//! Statistical hypothesis testing for Neptune-crossing TNO perihelion
//! distributions (Batygin, Morbidelli, Brown & Nesvorny 2024, Sec. 3.2).
//!
//! Each observed object j, discovered at heliocentric distance r_j, is an
//! unbiased probe of the model's perihelion distribution interior to r_j:
//! xi_j = CDF_{r_j}(q_j), where the CDF is built from *simulated* orbital
//! footprints, geometrically weighted by the fraction of each orbit's period
//! spent inside r_j (Brown 2001; Morbidelli et al. 2004). Under a correct
//! model the xi_j are U(0,1).
//!
//! Two statistics quantify departure from uniformity:
//! - Kolmogorov-Smirnov on the xi values: paper p = 0.41 (P9) vs p = 0.0034
//!   (P9-free). A one-sided Gaussian equivalent of p = 0.0034 is ~2.7 sigma.
//! - zeta = sum_j log10(xi_j), compared against its Monte-Carlo null
//!   (sum of 17 log10 U(0,1) draws: mean -17/ln10 = -7.38, sigma
//!   sqrt(17)/ln10 = 1.79; the paper quotes -7.2 and 1.8 from a smoothed
//!   million-draw histogram). zeta_P9 = -7.9 sits ~0.4 sigma from the mean;
//!   zeta_free = -16.5 sits ~5 sigma below it — this is the paper's
//!   "~5 sigma" statement. It is a deviation of zeta from the null mean,
//!   *not* a p = 0.0034 <-> 5 sigma conversion (p = 0.0034 is ~2.7 sigma).

use rand::Rng;
use serde::{Deserialize, Serialize};

use crate::observed_tnos::NeptuneCrossingTno;
use crate::simulation::Footprint;

/// Result of the hypothesis tests for a given model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HypothesisResult {
    /// Zeta statistic for the P9-inclusive model
    pub zeta_p9: f64,
    /// KS p-value for the P9-inclusive model
    pub p_value_p9: f64,
    /// Zeta statistic for the P9-free null model
    pub zeta_null: f64,
    /// KS p-value for the P9-free null model
    pub p_value_null: f64,
}

impl HypothesisResult {
    /// Returns the paper's reported values.
    pub fn paper_values() -> Self {
        Self {
            zeta_p9: -7.9,
            p_value_p9: 0.41,
            zeta_null: -16.5,
            p_value_null: 0.0034,
        }
    }
}

/// Fraction of the orbital period an orbit (a, e) spends at heliocentric
/// distance < r. From Kepler's equation with r = a(1 - e cos E):
/// t/T = (E_r - e sin E_r)/pi, E_r = acos((1 - r/a)/e).
pub fn time_fraction_inside(a: f64, e: f64, r: f64) -> f64 {
    let q = a * (1.0 - e);
    let big_q = a * (1.0 + e);
    if r <= q {
        return 0.0;
    }
    if r >= big_q {
        return 1.0;
    }
    let cos_e = ((1.0 - r / a) / e).clamp(-1.0, 1.0);
    let e_anom = cos_e.acos();
    (e_anom - e * e_anom.sin()) / std::f64::consts::PI
}

/// Empirical CDF of perihelion interior to discovery distance r, built from
/// simulated footprints with the geometric (time-inside-r) weighting:
///
///   CDF_r(q) = sum_k w_k(r) * [q_k <= q] / sum_k w_k(r),
///   w_k(r) = fraction of orbit k's period spent at distance < r.
///
/// Returns 0 if no footprint can appear inside r (the model assigns zero
/// probability to such a discovery; callers clamp before taking logs).
pub fn empirical_cdf_q(footprints: &[Footprint], q: f64, r: f64) -> f64 {
    let mut w_total = 0.0;
    let mut w_below = 0.0;
    for f in footprints {
        let e = 1.0 - f.q / f.a;
        if !(0.0..1.0).contains(&e) {
            continue;
        }
        let w = time_fraction_inside(f.a, e, r);
        w_total += w;
        if f.q <= q {
            w_below += w;
        }
    }
    if w_total <= 0.0 {
        return 0.0;
    }
    w_below / w_total
}

/// Heuristic discovery heliocentric distances (AU) — documented stand-in.
///
/// True discovery circumstances are not part of the SBDB orbit record. The
/// paper (footnote 6) notes that this sample's discovery magnitudes span
/// ~21.5-24.5 and that the steep r^-4 reflected-flux scaling biases
/// discovery toward perihelion. We therefore approximate the discovery
/// distance as r_disc = 1.15 * q. This is an *approximation with a stated
/// rationale*, not data; replace with real circumstances if they are
/// compiled.
pub fn approximate_discovery_distances(sample: &[NeptuneCrossingTno]) -> Vec<f64> {
    sample.iter().map(|tno| 1.15 * tno.q).collect()
}

/// xi_j = CDF_{r_j}(q_j) for each observed object, using a simulated
/// footprint population as the model.
pub fn xi_values(
    sample: &[NeptuneCrossingTno],
    discovery_distances: &[f64],
    footprints: &[Footprint],
) -> Vec<f64> {
    sample
        .iter()
        .zip(discovery_distances)
        .map(|(tno, &r)| empirical_cdf_q(footprints, tno.q, r))
        .collect()
}

/// zeta = sum_j log10(xi_j) (the paper's logarithmic statistic; xi clamped
/// away from 0 to keep the sum finite when a model assigns zero probability).
pub fn compute_zeta(xi: &[f64]) -> f64 {
    xi.iter().map(|&x| x.clamp(1e-15, 1.0).log10()).sum()
}

/// Monte-Carlo null distribution of zeta: sum of `n_obj` log10 U(0,1) draws
/// (the paper draws 10^6 of these). Seeded and reproducible.
pub fn zeta_null_samples(n_obj: usize, n_mc: usize, seed: u64) -> Vec<f64> {
    use rand::SeedableRng;
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    (0..n_mc)
        .map(|_| {
            (0..n_obj)
                .map(|_| rng.gen_range(1e-300_f64..1.0).log10())
                .sum()
        })
        .collect()
}

/// One-sided p-value of an observed zeta against the MC null (fraction of
/// null draws at least as negative).
pub fn zeta_p_value(zeta_obs: f64, null_samples: &[f64]) -> f64 {
    if null_samples.is_empty() {
        return f64::NAN;
    }
    let n_below = null_samples.iter().filter(|&&z| z <= zeta_obs).count();
    n_below as f64 / null_samples.len() as f64
}

/// Deviation of an observed zeta below the null mean, in null standard
/// deviations. The paper's "~5 sigma": (mean - (-16.5))/sigma ~ 5.
pub fn zeta_sigma_deviation(zeta_obs: f64, null_samples: &[f64]) -> f64 {
    let n = null_samples.len() as f64;
    let mean = null_samples.iter().sum::<f64>() / n;
    let var = null_samples.iter().map(|z| (z - mean).powi(2)).sum::<f64>() / (n - 1.0);
    (mean - zeta_obs) / var.sqrt()
}

/// Kolmogorov-Smirnov statistic of xi values against uniform [0,1]:
/// D_n = max(|F_n(x) - x|).
pub fn ks_test(xi_values: &[f64]) -> f64 {
    let n = xi_values.len() as f64;
    let mut sorted = xi_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let mut d_max = 0.0_f64;
    for (i, &xi) in sorted.iter().enumerate() {
        let f_n = (i + 1) as f64 / n;
        let d_plus = (f_n - xi).abs();
        let d_minus = (xi - i as f64 / n).abs();
        d_max = d_max.max(d_plus).max(d_minus);
    }
    d_max
}

/// KS p-value via the Kolmogorov distribution with the Stephens/Marsaglia
/// finite-n correction (Numerical Recipes eq. 14.3.9; Marsaglia, Tsang &
/// Wang 2003): lambda = (sqrt(n) + 0.12 + 0.11/sqrt(n)) D,
/// Q(lambda) = 2 sum_{k>=1} (-1)^{k-1} exp(-2 k^2 lambda^2).
/// Good to a few percent down to n ~ 5 (we use n = 17).
pub fn ks_p_value(d: f64, n: usize) -> f64 {
    if d <= 0.0 {
        return 1.0;
    }
    let sqrt_n = (n as f64).sqrt();
    let lambda = (sqrt_n + 0.12 + 0.11 / sqrt_n) * d;
    let mut sum = 0.0;
    let mut sign = 1.0;
    for k in 1..=100 {
        let term = (-2.0 * (k as f64).powi(2) * lambda * lambda).exp();
        sum += sign * term;
        sign = -sign;
        if term < 1e-12 {
            break;
        }
    }
    (2.0 * sum).clamp(0.0, 1.0)
}

/// Anderson-Darling test statistic for uniformity of xi values.
///
/// A2 = -n - (1/n) * sum_{j=1}^{n} [(2j-1)(ln(xi_j) + ln(1 - xi_{n+1-j}))]
pub fn anderson_darling_test(xi_values: &[f64]) -> f64 {
    let n = xi_values.len();
    let nf = n as f64;

    let mut sorted = xi_values.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    // Clamp to avoid log(0)
    let sorted: Vec<f64> = sorted
        .iter()
        .map(|&x| x.clamp(1e-15, 1.0 - 1e-15))
        .collect();

    let mut sum = 0.0;
    for j in 0..n {
        let weight = (2 * j + 1) as f64;
        sum += weight * (sorted[j].ln() + (1.0 - sorted[n - 1 - j]).ln());
    }

    -nf - sum / nf
}

/// Anderson-Darling p-value for the fully-specified (case 0) uniform null,
/// from the piecewise approximation of D'Agostino & Stephens (1986),
/// Table 4.9 / Sec. 4.10:
pub fn anderson_darling_p_value(a2: f64) -> f64 {
    let p = if a2 < 0.2 {
        1.0 - (-13.436 + 101.14 * a2 - 223.73 * a2 * a2).exp()
    } else if a2 < 0.34 {
        1.0 - (-8.318 + 42.796 * a2 - 59.938 * a2 * a2).exp()
    } else if a2 < 0.6 {
        (0.9177 - 4.279 * a2 - 1.38 * a2 * a2).exp()
    } else {
        (1.2937 - 5.709 * a2 + 0.0186 * a2 * a2).exp()
    };
    p.clamp(0.0, 1.0)
}

/// One-sided Gaussian sigma equivalent of a p-value (inverse survival
/// function). p = 0.0034 -> ~2.7 sigma. (The paper's "~5 sigma" statement is
/// the zeta deviation from the null mean — see `zeta_sigma_deviation` — not
/// this conversion.)
pub fn sigma_rejection(p_value: f64) -> f64 {
    // Abramowitz & Stegun 26.2.23 rational approximation for the inverse
    // normal CDF.
    if p_value <= 0.0 || p_value >= 1.0 {
        return 0.0;
    }

    let p = if p_value > 0.5 {
        1.0 - p_value
    } else {
        p_value
    };

    let t = (-2.0 * p.ln()).sqrt();
    let c0 = 2.515517;
    let c1 = 0.802853;
    let c2 = 0.010328;
    let d1 = 1.432788;
    let d2 = 0.189269;
    let d3 = 0.001308;

    let z = t - (c0 + c1 * t + c2 * t * t) / (1.0 + d1 * t + d2 * t * t + d3 * t * t * t);

    if p_value > 0.5 {
        -z
    } else {
        z
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::observed_tnos::observed_sample;
    use crate::simulation::quick_test_simulation;

    #[test]
    fn test_paper_values() {
        let result = HypothesisResult::paper_values();
        assert!((result.zeta_p9 - (-7.9)).abs() < 0.01);
        assert!((result.p_value_p9 - 0.41).abs() < 0.01);
        assert!((result.zeta_null - (-16.5)).abs() < 0.01);
        assert!((result.p_value_null - 0.0034).abs() < 0.001);
    }

    #[test]
    fn test_time_fraction_inside_limits() {
        assert_eq!(time_fraction_inside(150.0, 0.8, 20.0), 0.0); // r < q
        assert_eq!(time_fraction_inside(150.0, 0.8, 300.0), 1.0); // r > Q
        let f = time_fraction_inside(150.0, 0.8, 40.0);
        assert!(f > 0.0 && f < 0.5, "f = {f}");
        // Monotone in r
        let f2 = time_fraction_inside(150.0, 0.8, 80.0);
        assert!(f2 > f);
    }

    #[test]
    fn test_empirical_cdf_monotone_and_bounded() {
        let fp = quick_test_simulation(true).selected_footprints();
        assert!(!fp.is_empty(), "P9 run must produce crossing footprints");
        let r = 33.0;
        let mut prev = 0.0;
        for k in 0..=30 {
            let q = k as f64;
            let c = empirical_cdf_q(&fp, q, r);
            assert!((0.0..=1.0).contains(&c));
            assert!(c >= prev - 1e-12, "CDF must be monotone");
            prev = c;
        }
    }

    #[test]
    fn test_zeta_null_distribution_matches_theory() {
        // zeta_null = sum of 17 log10 U(0,1): mean = -17/ln10 = -7.383,
        // sigma = sqrt(17)/ln10 = 1.7906 (paper: -7.2, 1.8).
        let samples = zeta_null_samples(17, 40_000, 7);
        let n = samples.len() as f64;
        let mean = samples.iter().sum::<f64>() / n;
        let std = (samples.iter().map(|z| (z - mean).powi(2)).sum::<f64>() / (n - 1.0)).sqrt();
        assert!((mean - (-7.383)).abs() < 0.05, "MC mean = {mean}");
        assert!((std - 1.7906).abs() < 0.05, "MC sigma = {std}");
    }

    #[test]
    fn test_paper_zeta_free_is_about_5_sigma() {
        // Regression against the paper's "~5 sigma": zeta = -16.5 lies
        // ~(16.5 - 7.38)/1.79 = 5.1 null standard deviations below the mean.
        let samples = zeta_null_samples(17, 40_000, 7);
        let dev = zeta_sigma_deviation(-16.5, &samples);
        assert!((4.6..=5.6).contains(&dev), "deviation = {dev} sigma");
        // And the P9 model's zeta is well within the null bulk.
        let dev_p9 = zeta_sigma_deviation(-7.9, &samples);
        assert!(dev_p9.abs() < 1.0, "P9 deviation = {dev_p9} sigma");
    }

    #[test]
    fn test_zeta_direction_p9_beats_null_model() {
        // X1 regression (reduced scale): the observed 17 objects must be
        // *more* consistent (less negative zeta) with the P9-inclusive
        // simulation's perihelion CDF than with the P9-free control, which
        // at this scale produces essentially no low-q footprints.
        let sample = observed_sample();
        let r_disc = approximate_discovery_distances(&sample);

        let fp_p9 = quick_test_simulation(true).selected_footprints();
        let fp_free = quick_test_simulation(false).selected_footprints();

        let zeta_p9 = compute_zeta(&xi_values(&sample, &r_disc, &fp_p9));
        let zeta_free = compute_zeta(&xi_values(&sample, &r_disc, &fp_free));

        assert!(
            zeta_p9 > zeta_free,
            "zeta_P9 = {zeta_p9} should exceed zeta_free = {zeta_free}"
        );
        assert!(zeta_p9.is_finite() && zeta_free.is_finite());
    }

    #[test]
    fn test_ks_statistic_range() {
        let xi_values = vec![0.1, 0.25, 0.4, 0.55, 0.7, 0.85, 0.95];
        let d = ks_test(&xi_values);
        assert!((0.0..=1.0).contains(&d), "KS statistic = {}", d);
    }

    #[test]
    fn test_ks_uniform_sample() {
        let n = 10;
        let xi: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / n as f64).collect();
        let d = ks_test(&xi);
        assert!(d < 0.15, "KS D = {} for near-uniform sample", d);
    }

    #[test]
    fn test_ks_p_value_mapping() {
        // p(D -> 0) -> 1, p decreasing in D.
        assert!(ks_p_value(0.01, 17) > 0.99);
        assert!(ks_p_value(0.2, 17) > ks_p_value(0.3, 17));
        // The paper's p = 0.0034 at n = 17 corresponds to D ~ 0.42.
        let p = ks_p_value(0.42, 17);
        assert!(
            (0.001..0.01).contains(&p),
            "p(D=0.42, n=17) = {p}, expected ~0.0034"
        );
    }

    #[test]
    fn test_anderson_darling_uniform() {
        let n = 20;
        let xi: Vec<f64> = (0..n).map(|i| (i as f64 + 0.5) / n as f64).collect();
        let a2 = anderson_darling_test(&xi);
        assert!(a2 < 2.5, "A2 = {} for near-uniform sample", a2);
        // Near-uniform sample: large p-value; skewed sample: small p-value.
        assert!(anderson_darling_p_value(a2) > 0.05);
        let skewed: Vec<f64> = (0..n).map(|i| 0.02 + 0.1 * (i as f64) / n as f64).collect();
        let a2_skewed = anderson_darling_test(&skewed);
        assert!(anderson_darling_p_value(a2_skewed) < 0.01);
    }

    #[test]
    fn test_sigma_rejection_known_values() {
        // p = 0.05 => ~1.645 sigma (one-tail)
        let sigma = sigma_rejection(0.05);
        assert!(
            (sigma - 1.645).abs() < 0.05,
            "sigma(0.05) = {}, expected ~1.645",
            sigma
        );

        // p = 0.0034 => ~2.7 sigma (one-tail) — the paper's KS rejection.
        let sigma = sigma_rejection(0.0034);
        assert!(
            (sigma - 2.7).abs() < 0.1,
            "sigma(0.0034) = {sigma}, expected ~2.7"
        );
    }
}
