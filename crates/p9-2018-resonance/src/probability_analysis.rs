//! Resonance probability analysis and a₉ distribution computation.
//!
//! Implements the statistical analysis from Bailey+ (2018) Section 4:
//! - Compute P(N/1 or N/2) for each e₉ value
//! - Show that P(all 6 KBOs in N/1 or N/2) < 5%
//! - Compare a₉ distributions using Farey F5 vs full resonance spectrum
//!
//! Key result: the prominent peak at a₉ ≈ 660 AU (Millholland & Laughlin
//! 2017) dissolves into a broad plateau when high-order resonances are
//! included.
//!
//! Observational input: the vetted ETNO table in `p9_core::data::etno`
//! (the local `OBSERVED_KBO_AXES` copy this replaces had scrambled
//! semi-major axes — e.g. 2013 RF98 at 780 AU instead of ~325 — displacing
//! every implied-a₉ peak by tens of AU).

use crate::resonance_catalog::{extended_catalog, farey_f5, Resonance};

/// Names of the 6 clustered ETNOs of Batygin & Brown (2016) — the sample
/// Millholland & Laughlin (2017) and Bailey+ (2018) analyze ("all 6 KBOs").
const BB16_SAMPLE: [&str; 6] = [
    "Sedna",
    "2012 VP113",
    "2004 VN112",
    "2010 GB174",
    "2007 TG422",
    "2013 RF98",
];

/// Semi-major axes (AU) of the observed 6-KBO sample, drawn from the single
/// vetted workspace table in `p9_core::data::etno`.
pub fn observed_kbo_axes() -> Vec<f64> {
    p9_core::data::etno::BROWN_2017_SAMPLE
        .iter()
        .filter(|k| BB16_SAMPLE.contains(&k.name))
        .map(|k| k.a)
        .collect()
}

/// Compute the probability that a randomly chosen resonant particle
/// is in an N/1 or N/2 resonance, given the resonance catalog.
///
/// Uses the occupation weights from a simulation census.
pub fn p_simple_resonance(census: &[(Resonance, usize)]) -> f64 {
    let total: usize = census.iter().map(|(_, c)| c).sum();
    if total == 0 {
        return 0.0;
    }

    let simple: usize = census
        .iter()
        .filter(|(r, _)| r.is_simple())
        .map(|(_, c)| c)
        .sum();

    simple as f64 / total as f64
}

/// Compute P(all n KBOs in N/1 or N/2) = p_simple^n.
pub fn p_all_simple(p_simple: f64, n_kbos: usize) -> f64 {
    p_simple.powi(n_kbos as i32)
}

/// Compute the implied a₉ given a KBO semi-major axis and a period ratio.
///
/// If the KBO is in a p:q resonance with P9:
///   a_kbo = a₉ * (q/p)^(2/3)
///   a₉ = a_kbo / (q/p)^(2/3) = a_kbo * (p/q)^(2/3)
pub fn implied_a9(a_kbo: f64, res: &Resonance) -> f64 {
    a_kbo * (res.p as f64 / res.q as f64).powf(2.0 / 3.0)
}

/// Explicit, honest prior on a₉ and the kernel model for the
/// commensurability scan. Every parameter is a documented assumption (the
/// previous implementation hid a hard a₉ > 500 AU step inside an expression
/// dressed up as a mass marginalization).
#[derive(Debug, Clone)]
pub struct A9Prior {
    /// Lower edge of the uniform a₉ prior (AU). Dynamical-stability and
    /// survey arguments exclude a₉ ≲ 300 AU for a body shepherding ETNOs
    /// with a up to ~500 AU.
    pub a9_min: f64,
    /// Upper edge of the uniform a₉ prior (AU): beyond ~1000 AU the
    /// perturber is too distant to maintain the observed clustering.
    pub a9_max: f64,
    /// Fractional 1σ uncertainty on the observed KBO semi-major axes,
    /// setting the kernel width σ = frac · a_kbo of the commensurability
    /// score. ~2% is conservative for the multi-opposition orbits in the
    /// vetted table.
    pub a_kbo_fractional_sigma: f64,
    /// When true, weight each commensurability by 1/order (a proxy for
    /// resonance width/occupancy) instead of equally.
    pub order_weighting: bool,
}

impl Default for A9Prior {
    fn default() -> Self {
        Self {
            a9_min: 300.0,
            a9_max: 1000.0,
            a_kbo_fractional_sigma: 0.02,
            order_weighting: false,
        }
    }
}

/// Result of a₉ distribution computation.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct A9Distribution {
    /// Bin centers (AU)
    pub bins: Vec<f64>,
    /// Probability density
    pub density: Vec<f64>,
    /// Label for this distribution
    pub label: String,
}

/// Compute the a₉ commensurability distribution for a given resonance
/// catalog (the Millholland & Laughlin 2017 scan, generalized to an
/// arbitrary catalog as in Bailey+ 2018).
///
/// For each trial a₉ on the explicit prior range, each KBO is scored by its
/// *best-fitting* commensurability: max over interior resonances of
/// w · exp(−(a_kbo − a₉ (q/p)^{2/3})² / 2σ²) with σ = frac · a_kbo from the
/// KBO's a-uncertainty. Scoring the best match per KBO (rather than summing
/// a kernel per resonance) is what makes the scan a consensus statistic:
/// the peak appears where *all* KBOs are simultaneously near low-order
/// commensurabilities, and adding high-order resonances to the catalog
/// dilutes it toward a plateau — the paper's headline comparison.
///
/// Only interior resonances (period ratio < 1; the perturber exterior to
/// the clustered ETNOs, as the Planet Nine hypothesis requires) are scored.
pub fn compute_a9_distribution(
    kbo_axes: &[f64],
    catalog: &[Resonance],
    prior: &A9Prior,
    n_bins: usize,
    label: &str,
) -> A9Distribution {
    let da = (prior.a9_max - prior.a9_min) / n_bins as f64;
    let bins: Vec<f64> = (0..n_bins)
        .map(|i| prior.a9_min + (i as f64 + 0.5) * da)
        .collect();
    let mut density = vec![0.0_f64; n_bins];

    for (j, &a9) in bins.iter().enumerate() {
        let mut score = 0.0;
        for &a_kbo in kbo_axes {
            let sigma = prior.a_kbo_fractional_sigma * a_kbo;
            let mut best = 0.0_f64;
            for res in catalog {
                if res.period_ratio() >= 1.0 {
                    continue; // exterior: P9 inside the KBO orbit
                }
                let a_res = res.semimajor_axis(a9);
                let z = (a_kbo - a_res) / sigma;
                let weight = if prior.order_weighting {
                    1.0 / res.order().max(1) as f64
                } else {
                    1.0
                };
                best = best.max(weight * (-0.5 * z * z).exp());
            }
            score += best;
        }
        density[j] = score;
    }

    // Normalize to a probability density over the prior range.
    let total: f64 = density.iter().sum::<f64>() * da;
    if total > 0.0 {
        for d in &mut density {
            *d /= total;
        }
    }

    A9Distribution {
        bins,
        density,
        label: label.to_string(),
    }
}

/// Compute a₉ distribution using the Farey F5 sequence (Millholland & Laughlin 2017).
pub fn a9_distribution_farey_f5(kbo_axes: &[f64]) -> A9Distribution {
    compute_a9_distribution(kbo_axes, &farey_f5(), &A9Prior::default(), 100, "Farey F5")
}

/// Compute a₉ distribution using the full extended catalog (Bailey+ 2018).
pub fn a9_distribution_extended(kbo_axes: &[f64]) -> A9Distribution {
    compute_a9_distribution(
        kbo_axes,
        &extended_catalog(),
        &A9Prior::default(),
        100,
        "Extended (Bailey+ 2018)",
    )
}

/// Compare Farey F5 vs extended resonance distributions.
///
/// Returns (f5_peak, f5_max, extended_peak, extended_max, ratio).
/// If the extended distribution has a much lower peak-to-mean ratio,
/// the constraint is dissolved.
pub fn compare_distributions(kbo_axes: &[f64]) -> DistributionComparison {
    let f5 = a9_distribution_farey_f5(kbo_axes);
    let ext = a9_distribution_extended(kbo_axes);

    let f5_max = f5.density.iter().copied().fold(0.0_f64, f64::max);
    let f5_mean = f5.density.iter().sum::<f64>() / f5.density.len() as f64;
    let f5_peak_idx = f5
        .density
        .iter()
        .position(|&d| (d - f5_max).abs() < 1e-15)
        .unwrap_or(0);

    let ext_max = ext.density.iter().copied().fold(0.0_f64, f64::max);
    let ext_mean = ext.density.iter().sum::<f64>() / ext.density.len() as f64;
    let ext_peak_idx = ext
        .density
        .iter()
        .position(|&d| (d - ext_max).abs() < 1e-15)
        .unwrap_or(0);

    DistributionComparison {
        f5_peak_a9: f5.bins[f5_peak_idx],
        f5_peak_to_mean: if f5_mean > 0.0 { f5_max / f5_mean } else { 0.0 },
        ext_peak_a9: ext.bins[ext_peak_idx],
        ext_peak_to_mean: if ext_mean > 0.0 {
            ext_max / ext_mean
        } else {
            0.0
        },
        f5_distribution: f5,
        ext_distribution: ext,
    }
}

/// Comparison between F5 and extended resonance distributions.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DistributionComparison {
    /// Peak a₉ from F5 distribution
    pub f5_peak_a9: f64,
    /// Peak-to-mean ratio for F5 (higher = more constraining)
    pub f5_peak_to_mean: f64,
    /// Peak a₉ from extended distribution
    pub ext_peak_a9: f64,
    /// Peak-to-mean ratio for extended (lower = less constraining)
    pub ext_peak_to_mean: f64,
    /// Full F5 distribution
    pub f5_distribution: A9Distribution,
    /// Full extended distribution
    pub ext_distribution: A9Distribution,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_implied_a9() {
        // If KBO at a=378 is in 2:1 with P9:
        // a₉ = 378 * (2/1)^(2/3) = 378 * 1.587 ≈ 600
        let a9 = implied_a9(378.0, &Resonance::new(2, 1));
        assert!((a9 - 600.0).abs() < 5.0, "Implied a₉ = {:.1}", a9);
    }

    #[test]
    fn test_observed_axes_vetted() {
        // Single source: the p9-core ETNO table (Sedna 506, 2013 RF98 ~325 —
        // not the scrambled 472/780 of the deleted local copy), restricted
        // to the 6-object Batygin & Brown (2016) sample the papers analyze.
        let axes = observed_kbo_axes();
        assert_eq!(axes.len(), 6);
        assert!(axes.contains(&506.0));
        assert!(axes.contains(&325.0));
        assert!(!axes.contains(&780.0));
        assert!(!axes.contains(&472.0));
    }

    #[test]
    fn test_p_simple_resonance() {
        let census = vec![
            (Resonance::new(1, 3), 10), // N/1 (p=1, period ratio 3)
            (Resonance::new(5, 3), 5),  // High-order (p=5)
            (Resonance::new(2, 5), 5),  // N/2 (p=2, period ratio 5/2)
        ];

        let p = p_simple_resonance(&census);
        assert!((p - 0.75).abs() < 1e-10, "P(simple) = {}", p);
    }

    #[test]
    fn test_p_all_simple() {
        let p = p_all_simple(0.3, 6);
        assert!(p < 0.001, "P(all 6) = {:.4}", p);
    }

    #[test]
    fn test_a9_distribution_f5() {
        let dist = a9_distribution_farey_f5(&observed_kbo_axes());
        assert!(!dist.bins.is_empty());
        assert!(!dist.density.is_empty());

        let total: f64 = dist.density.iter().sum::<f64>();
        assert!(total > 0.0, "Distribution should have positive mass");
    }

    #[test]
    fn test_f5_peak_near_paper_value() {
        // Paper-number regression (Millholland & Laughlin 2017): the F5
        // implied-a₉ distribution peaks near 660 AU.
        let comparison = compare_distributions(&observed_kbo_axes());
        assert!(
            (comparison.f5_peak_a9 - 660.0).abs() < 60.0,
            "F5 peak at {:.0} AU, expected ~660 AU",
            comparison.f5_peak_a9
        );
    }

    #[test]
    fn test_extended_distribution_flatter() {
        // The main result of Bailey+ (2018): including the full resonance
        // spectrum makes the implied-a₉ distribution strictly *less* peaked
        // than the F5 baseline. (The previous assertion allowed the extended
        // peak to be 1.5x larger — vacuous.)
        let comparison = compare_distributions(&observed_kbo_axes());
        assert!(
            comparison.ext_peak_to_mean < 0.8 * comparison.f5_peak_to_mean,
            "Extended peak/mean ({:.2}) should be well below F5 ({:.2})",
            comparison.ext_peak_to_mean,
            comparison.f5_peak_to_mean,
        );
    }
}
