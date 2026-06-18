//! Population comparison between IOC-injected and scattered disk objects.
//!
//! Implements the key observational predictions from Batygin & Brown (2021):
//! - IOC-injected objects show weaker longitude-of-perihelion confinement
//!   than scattered disk objects (the paper's full-scale ensemble gives
//!   ~67% vs ~88%; this crate measures both from reduced-scale runs and a
//!   genuine P9-free control)
//! - IOC injection preferentially populates the large-a region
//!
//! The semi-major-axis comparison histogram is filled from the *simulated*
//! populations carried in `InjectionResult` (the previous version hard-coded
//! the scattered-disk counts).

use serde::{Deserialize, Serialize};

use p9_core::units::{au, Length};

use crate::injection_simulation::{simulate_injection, InjectionConfig, InjectionResult};

/// Comparison of confinement between IOC-injected and scattered disk populations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PopulationComparison {
    /// f_varpi for the IOC-injected population (end-state, measured)
    pub f_varpi_ioc: f64,
    /// f_varpi for the scattered disk population (end-state, measured)
    pub f_varpi_scattered: f64,
    /// f_varpi for the P9-free control run (~0.5 expected)
    pub f_varpi_control: f64,
    /// Number of IOC-injected objects in the sample
    pub n_ioc_injected: usize,
    /// Total number of IOC objects simulated
    pub n_ioc_total: usize,
    /// Injection efficiency (fraction of injectable IOC objects injected)
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
pub fn compare_confinement<R: rand::Rng>(
    config: &InjectionConfig,
    rng: &mut R,
) -> PopulationComparison {
    let result = simulate_injection(config, rng);

    PopulationComparison {
        f_varpi_ioc: result.f_varpi_ioc,
        f_varpi_scattered: result.f_varpi_scattered,
        f_varpi_control: result.f_varpi_control,
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

impl SmaBin {
    /// Lower bin edge as a [`Length`].
    pub fn lower_edge(&self) -> Length {
        au(self.a_low)
    }
    /// Upper bin edge as a [`Length`].
    pub fn upper_edge(&self) -> Length {
        au(self.a_high)
    }
}

/// Build a semi-major axis distribution histogram from the simulated
/// populations in an injection result (injected IOC end-state a vs the
/// simulated scattered-disk end-state a).
pub fn semi_major_axis_distribution(result: &InjectionResult, n_bins: usize) -> Vec<SmaBin> {
    let a_min = 150.0;
    let a_max = result
        .injected_sma
        .iter()
        .chain(result.scattered_sma.iter())
        .cloned()
        .fold(1000.0_f64, f64::max)
        * 1.05;
    let bin_width = (a_max - a_min) / n_bins as f64;

    let mut bins: Vec<SmaBin> = (0..n_bins)
        .map(|i| SmaBin {
            a_low: a_min + i as f64 * bin_width,
            a_high: a_min + (i + 1) as f64 * bin_width,
            count_ioc: 0,
            count_scattered: 0,
        })
        .collect();

    let mut fill = |a: f64, scattered: bool| {
        if a >= a_min && a < a_max {
            let idx = ((a - a_min) / bin_width) as usize;
            if idx < n_bins {
                if scattered {
                    bins[idx].count_scattered += 1;
                } else {
                    bins[idx].count_ioc += 1;
                }
            }
        }
    };
    for &a in &result.injected_sma {
        fill(a, false);
    }
    for &a in &result.scattered_sma {
        fill(a, true);
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
    use approx::assert_relative_eq;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_bin_edges_match_f64() {
        let bin = SmaBin {
            a_low: 150.0,
            a_high: 220.0,
            count_ioc: 3,
            count_scattered: 5,
        };
        assert_relative_eq!(
            bin.lower_edge().get::<astronomical_unit>(),
            bin.a_low,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            bin.upper_edge().get::<astronomical_unit>(),
            bin.a_high,
            epsilon = 1e-9
        );
    }

    fn sample_result() -> InjectionResult {
        InjectionResult {
            injection_fraction: 0.1,
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            f_varpi_control: 0.5,
            n_injected: 5,
            n_total: 50,
            injected_sma: vec![500.0, 3000.0, 5000.0, 8000.0, 15000.0],
            injected_dvarpi: vec![0.1, -0.5, 0.3, -0.2, 0.7],
            scattered_sma: vec![200.0, 280.0, 350.0, 420.0, 500.0, 540.0],
        }
    }

    #[test]
    fn test_confinement_ratio() {
        let comparison = PopulationComparison {
            f_varpi_ioc: 0.67,
            f_varpi_scattered: 0.88,
            f_varpi_control: 0.5,
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
    fn test_semi_major_axis_distribution_from_simulated_populations() {
        let result = sample_result();
        let bins = semi_major_axis_distribution(&result, 10);
        assert_eq!(bins.len(), 10);

        let total_ioc: usize = bins.iter().map(|b| b.count_ioc).sum();
        let total_sd: usize = bins.iter().map(|b| b.count_scattered).sum();
        assert_eq!(total_ioc, result.injected_sma.len());
        assert_eq!(total_sd, result.scattered_sma.len());

        // The scattered population sits at small a, the IOC reaches large a.
        let last_half_sd: usize = bins[5..].iter().map(|b| b.count_scattered).sum();
        assert_eq!(last_half_sd, 0);
        let last_half_ioc: usize = bins[5..].iter().map(|b| b.count_ioc).sum();
        assert!(last_half_ioc > 0);

        // Verify bin edges are contiguous
        for i in 1..bins.len() {
            assert!((bins[i].a_low - bins[i - 1].a_high).abs() < 1e-10);
        }
    }

    #[test]
    fn test_fraction_above_threshold() {
        let result = sample_result();
        let f = fraction_above_threshold(&result, 2000.0);
        assert!((f - 0.8).abs() < 1e-10);

        let f_all = fraction_above_threshold(&result, 100.0);
        assert!((f_all - 1.0).abs() < 1e-10);

        let empty = InjectionResult {
            injection_fraction: 0.0,
            f_varpi_ioc: f64::NAN,
            f_varpi_scattered: f64::NAN,
            f_varpi_control: f64::NAN,
            n_injected: 0,
            n_total: 0,
            injected_sma: vec![],
            injected_dvarpi: vec![],
            scattered_sma: vec![],
        };
        assert!((fraction_above_threshold(&empty, 1000.0) - 0.0).abs() < 1e-10);
    }
}
