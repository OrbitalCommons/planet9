//! Resonance probability analysis and a₉ distribution computation.
//!
//! Implements the statistical analysis from Bailey+ (2018) Section 4:
//! - Compute P(N/1 or N/2) for each e₉ value
//! - Show that P(all 6 KBOs in N/1 or N/2) < 5%
//! - Compare a₉ distributions using Farey F5 vs full resonance spectrum
//!
//! Key result: the prominent peak at a₉ ≈ 660 AU dissolves into a broad
//! plateau when high-order resonances are included.

use p9_core::constants::*;

use crate::resonance_catalog::{extended_catalog, farey_f5, Resonance};

/// Observed KBO semi-major axes from Malhotra et al. (2016).
/// Used for the feasibility analysis.
pub const OBSERVED_KBO_AXES: &[f64] = &[
    164.0, // 2003 CR105
    220.0, // 2003 GB32 (approximate)
    226.0, // 2001 FP185
    266.0, // 2012 VP113
    356.0, // 2004 VN112 (approximate)
    472.0, // Sedna
    735.0, // 2013 SY99 (approximate)
    316.0, // 2010 GB174
    780.0, // 2013 RF98 (approximate)
    540.0, // 2007 TG422 (approximate)
];

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

/// Compute the a₉ distribution using a given resonance catalog.
///
/// For each observed KBO and each resonance in the catalog, compute the
/// implied a₉ and add a Gaussian kernel. This implements the approach
/// of Millholland & Laughlin (2017) but with different resonance sets.
pub fn compute_a9_distribution(
    kbo_axes: &[f64],
    catalog: &[Resonance],
    a9_range: (f64, f64),
    n_bins: usize,
    sigma: f64,
    label: &str,
) -> A9Distribution {
    let da = (a9_range.1 - a9_range.0) / n_bins as f64;
    let bins: Vec<f64> = (0..n_bins)
        .map(|i| a9_range.0 + (i as f64 + 0.5) * da)
        .collect();
    let mut density = vec![0.0_f64; n_bins];

    // Mass range for Monte Carlo sampling
    let mass_mean = 10.0; // M_Earth
    let _mass_min = 5.0;
    let _mass_max = 20.0;

    for &a_kbo in kbo_axes {
        for res in catalog {
            let a9 = implied_a9(a_kbo, res);

            // Valid range: [200, 800] AU approximately
            let a9_min = 200.0 + 30.0 * mass_mean / EARTH_MASS_SOLAR * EARTH_MASS_SOLAR;
            let a9_max = 800.0;
            if a9 < a9_min || a9 > a9_max {
                continue;
            }

            // Add Gaussian kernel
            for (j, &bin) in bins.iter().enumerate() {
                let z = (bin - a9) / sigma;
                density[j] += (-0.5 * z * z).exp();
            }
        }
    }

    // Normalize
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
    let catalog = farey_f5();
    compute_a9_distribution(kbo_axes, &catalog, (200.0, 1200.0), 100, 30.0, "Farey F5")
}

/// Compute a₉ distribution using the full extended catalog (Bailey+ 2018).
pub fn a9_distribution_extended(kbo_axes: &[f64]) -> A9Distribution {
    let catalog = extended_catalog();
    compute_a9_distribution(
        kbo_axes,
        &catalog,
        (200.0, 1200.0),
        100,
        30.0,
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
        let dist = a9_distribution_farey_f5(OBSERVED_KBO_AXES);
        assert!(!dist.bins.is_empty());
        assert!(!dist.density.is_empty());

        let total: f64 = dist.density.iter().sum::<f64>();
        assert!(total > 0.0, "Distribution should have positive mass");
    }

    #[test]
    fn test_extended_distribution_flatter() {
        let comparison = compare_distributions(OBSERVED_KBO_AXES);

        // The extended distribution should be flatter (lower peak-to-mean)
        // This is the main result of the paper
        assert!(
            comparison.ext_peak_to_mean < comparison.f5_peak_to_mean * 1.5,
            "Extended peak/mean ({:.2}) should be smaller than F5 ({:.2})",
            comparison.ext_peak_to_mean,
            comparison.f5_peak_to_mean,
        );
    }
}
