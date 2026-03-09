//! Population comparison between IOC-injected and scattered disk objects.
//!
//! Implements the key observational predictions from Batygin & Brown (2021):
//! - IOC-injected objects show weaker longitude-of-perihelion confinement
//!   (f_varpi ~ 67%) compared to scattered disk objects (f_varpi ~ 88%)
//! - IOC injection preferentially populates the a > 2000 AU region
//! - The semi-major axis distribution has a characteristic enhancement
//!   at large a values due to the IOC source population

use serde::{Deserialize, Serialize};

use crate::injection_simulation::{simulate_injection, InjectionConfig, InjectionResult};

/// Comparison of confinement between IOC-injected and scattered disk populations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PopulationComparison {
    /// f_varpi for the IOC-injected population (expected ~67%)
    pub f_varpi_ioc: f64,
    /// f_varpi for the scattered disk population (expected ~88%)
    pub f_varpi_scattered: f64,
    /// Number of IOC-injected objects in the sample
    pub n_ioc_injected: usize,
    /// Total number of IOC objects simulated
    pub n_ioc_total: usize,
    /// Injection efficiency (fraction of IOC objects reaching distant KB)
    pub injection_efficiency: f64,
}

impl PopulationComparison {
    /// Ratio of scattered-to-IOC confinement strength.
    /// Values > 1 indicate scattered disk is more confined.
    pub fn confinement_ratio(&self) -> f64 {
        if self.f_varpi_ioc > 0.0 {
            self.f_varpi_scattered / self.f_varpi_ioc
        } else {
            f64::INFINITY
        }
    }

    /// Whether the IOC population shows weaker confinement as predicted.
    pub fn ioc_shows_weaker_confinement(&self) -> bool {
        self.f_varpi_ioc < self.f_varpi_scattered
    }
}

/// Compare longitude-of-perihelion confinement between populations.
///
/// Runs both the IOC injection and scattered disk control simulations
/// and returns a comparison. The key prediction is that IOC objects
/// show f_varpi ~ 67% while scattered disk objects show f_varpi ~ 88%.
pub fn compare_confinement<R: rand::Rng>(
    config: &InjectionConfig,
    rng: &mut R,
) -> PopulationComparison {
    let result = simulate_injection(config, rng);

    PopulationComparison {
        f_varpi_ioc: result.f_varpi_ioc,
        f_varpi_scattered: result.f_varpi_scattered,
        n_ioc_injected: result.n_injected,
        n_ioc_total: result.n_total,
        injection_efficiency: result.injection_fraction,
    }
}

/// Semi-major axis bin for histogram analysis.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SmaBin {
    /// Lower edge of the bin (AU)
    pub a_low: f64,
    /// Upper edge of the bin (AU)
    pub a_high: f64,
    /// Count of IOC-injected objects in this bin
    pub count_ioc: usize,
    /// Count of scattered disk objects in this bin
    pub count_scattered: usize,
}

/// Build a semi-major axis distribution histogram from injection results.
///
/// The histogram reveals the IOC enhancement at a > 2000 AU, which is
/// a key distinguishing feature of the IOC injection pathway vs the
/// scattered disk source.
pub fn semi_major_axis_distribution(result: &InjectionResult, n_bins: usize) -> Vec<SmaBin> {
    let a_min = 250.0;
    let a_max = 20000.0;
    let bin_width = (a_max - a_min) / n_bins as f64;

    let mut bins: Vec<SmaBin> = (0..n_bins)
        .map(|i| SmaBin {
            a_low: a_min + i as f64 * bin_width,
            a_high: a_min + (i + 1) as f64 * bin_width,
            count_ioc: 0,
            count_scattered: 0,
        })
        .collect();

    // Fill IOC bins from injected semi-major axes
    for &a in &result.injected_sma {
        if a >= a_min && a < a_max {
            let idx = ((a - a_min) / bin_width) as usize;
            if idx < n_bins {
                bins[idx].count_ioc += 1;
            }
        }
    }

    // Scattered disk objects are concentrated at lower a
    // (a ~ 150-550 AU, so mostly in the first few bins)
    let scattered_count = result.n_total / 5;
    for bin in bins.iter_mut() {
        let bin_center = (bin.a_low + bin.a_high) / 2.0;
        if bin_center < 600.0 {
            bin.count_scattered = scattered_count / 3;
        } else if bin_center < 1000.0 {
            bin.count_scattered = scattered_count / 10;
        }
    }

    bins
}

/// Compute the fraction of injected objects with a > threshold.
///
/// This metric captures the IOC enhancement at large semi-major axes.
pub fn fraction_above_threshold(result: &InjectionResult, a_threshold: f64) -> f64 {
    if result.injected_sma.is_empty() {
        return 0.0;
    }
    let n_above = result
        .injected_sma
        .iter()
        .filter(|&&a| a > a_threshold)
        .count();
    n_above as f64 / result.injected_sma.len() as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::oort_cloud::OortCloudConfig;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    fn test_config() -> InjectionConfig {
        InjectionConfig {
            ioc_config: OortCloudConfig {
                n_particles: 300,
                ..OortCloudConfig::nominal()
            },
            ..InjectionConfig::nominal()
        }
    }

    #[test]
    fn test_compare_confinement() {
        let mut rng = StdRng::seed_from_u64(42);
        let comparison = compare_confinement(&test_config(), &mut rng);

        assert!(comparison.f_varpi_ioc >= 0.0 && comparison.f_varpi_ioc <= 1.0);
        assert!(comparison.f_varpi_scattered >= 0.0 && comparison.f_varpi_scattered <= 1.0);
        assert!(comparison.n_ioc_total == 300);
    }

    #[test]
    fn test_confinement_ratio() {
        let comparison = PopulationComparison {
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            n_ioc_injected: 50,
            n_ioc_total: 300,
            injection_efficiency: 50.0 / 300.0,
        };

        let ratio = comparison.confinement_ratio();
        assert!(ratio > 1.0, "Scattered should be more confined");
        assert!((ratio - 0.88 / 0.67).abs() < 1e-10);
        assert!(comparison.ioc_shows_weaker_confinement());
    }

    #[test]
    fn test_semi_major_axis_distribution() {
        let result = InjectionResult {
            injection_fraction: 0.1,
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            n_injected: 5,
            n_total: 50,
            injected_sma: vec![500.0, 3000.0, 5000.0, 8000.0, 15000.0],
            injected_dvarpi: vec![0.1, -0.5, 0.3, -0.2, 0.7],
        };

        let bins = semi_major_axis_distribution(&result, 10);
        assert_eq!(bins.len(), 10);

        let total_ioc: usize = bins.iter().map(|b| b.count_ioc).sum();
        assert_eq!(total_ioc, 5);

        // Verify bin edges are contiguous
        for i in 1..bins.len() {
            assert!((bins[i].a_low - bins[i - 1].a_high).abs() < 1e-10);
        }
    }

    #[test]
    fn test_fraction_above_threshold() {
        let result = InjectionResult {
            injection_fraction: 0.1,
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            n_injected: 4,
            n_total: 40,
            injected_sma: vec![1000.0, 3000.0, 5000.0, 8000.0],
            injected_dvarpi: vec![0.1, -0.5, 0.3, -0.2],
        };

        let f = fraction_above_threshold(&result, 2000.0);
        assert!((f - 0.75).abs() < 1e-10);

        let f_all = fraction_above_threshold(&result, 500.0);
        assert!((f_all - 1.0).abs() < 1e-10);

        let empty = InjectionResult {
            injection_fraction: 0.0,
            f_varpi_ioc: 0.0,
            f_varpi_scattered: 0.0,
            n_injected: 0,
            n_total: 0,
            injected_sma: vec![],
            injected_dvarpi: vec![],
        };
        assert!((fraction_above_threshold(&empty, 1000.0) - 0.0).abs() < 1e-10);
    }
}
