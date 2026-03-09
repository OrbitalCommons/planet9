//! Kozai-Lidov oscillation detection and pathway analysis.
//!
//! The paper identifies a two-step pathway for generating high-i TNOs:
//! 1. Planet Nine drives Kozai-Lidov oscillations (coupled e-i oscillations)
//! 2. Neptune scattering reduces semi-major axis at high inclination
//!
//! This module provides tools to detect and characterize these oscillations
//! in simulation time-series data.

use p9_core::constants::*;
use p9_core::types::OrbitalElements;

/// Time series record for a single particle.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ParticleHistory {
    pub id: usize,
    pub times: Vec<f64>,
    pub elements: Vec<OrbitalElements>,
}

/// Classification of a particle's dynamical pathway.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Pathway {
    /// Still in the scattered disk (a > 100, i < 50°)
    ScatteredDisk,
    /// Undergoing Kozai-Lidov oscillations (a > 100, i oscillating)
    KozaiLidov,
    /// Decoupled high-i object (a < 100, i > 50°)
    HighInclination,
    /// Ejected or accreted
    Removed,
}

/// Detect Kozai-Lidov oscillations in a particle's history.
///
/// Kozai-Lidov signature: anti-correlated oscillations in e and i
/// with period much longer than orbital period.
pub fn detect_kozai_lidov(history: &ParticleHistory) -> bool {
    if history.elements.len() < 10 {
        return false;
    }

    // Look for inclination oscillations with amplitude > 20°
    let i_values: Vec<f64> = history.elements.iter().map(|e| e.i * RAD2DEG).collect();
    let i_min = i_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let i_max = i_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if i_max - i_min < 20.0 {
        return false;
    }

    // Check for anti-correlation between e and cos(i)
    // In Kozai-Lidov: sqrt(1-e²)*cos(i) ≈ const (Kozai integral)
    let kozai_values: Vec<f64> = history
        .elements
        .iter()
        .map(|e| (1.0 - e.e * e.e).sqrt() * e.i.cos())
        .collect();

    let k_mean = kozai_values.iter().sum::<f64>() / kozai_values.len() as f64;
    let k_var = kozai_values
        .iter()
        .map(|&k| (k - k_mean).powi(2))
        .sum::<f64>()
        / kozai_values.len() as f64;

    // If Kozai integral is approximately conserved, std should be small
    let k_std = k_var.sqrt();

    // Kozai-Lidov: large i oscillation + approximately conserved Kozai integral
    // k_std < 0.2 means the integral is roughly constant (allowing for Neptune perturbations)
    k_std < 0.2
}

/// Classify a particle based on its final state.
pub fn classify_particle(elem: &OrbitalElements) -> Pathway {
    let i_deg = elem.i * RAD2DEG;

    if elem.a > 100.0 {
        if i_deg > 50.0 {
            Pathway::KozaiLidov
        } else {
            Pathway::ScatteredDisk
        }
    } else if i_deg > 50.0 {
        Pathway::HighInclination
    } else {
        Pathway::ScatteredDisk
    }
}

/// Count particles in each pathway classification.
pub fn pathway_census(elements: &[OrbitalElements]) -> (usize, usize, usize) {
    let mut scattered = 0;
    let mut kozai = 0;
    let mut high_i = 0;

    for elem in elements {
        match classify_particle(elem) {
            Pathway::ScatteredDisk => scattered += 1,
            Pathway::KozaiLidov => kozai += 1,
            Pathway::HighInclination => high_i += 1,
            Pathway::Removed => {}
        }
    }

    (scattered, kozai, high_i)
}

/// Check if a particle matches the orbital region of known high-i TNOs.
///
/// Returns true if a < 100 AU, i > 60°, and q < 40 AU.
pub fn matches_high_i_region(elem: &OrbitalElements) -> bool {
    let q = elem.a * (1.0 - elem.e);
    elem.a < 100.0 && elem.i > 60.0 * DEG2RAD && q < 40.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_classify_scattered_disk() {
        let elem = OrbitalElements {
            a: 300.0,
            e: 0.7,
            i: 20.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert_eq!(classify_particle(&elem), Pathway::ScatteredDisk);
    }

    #[test]
    fn test_classify_high_inclination() {
        let elem = OrbitalElements {
            a: 50.0,
            e: 0.5,
            i: 110.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert_eq!(classify_particle(&elem), Pathway::HighInclination);
    }

    #[test]
    fn test_classify_kozai() {
        let elem = OrbitalElements {
            a: 200.0,
            e: 0.8,
            i: 80.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert_eq!(classify_particle(&elem), Pathway::KozaiLidov);
    }

    #[test]
    fn test_pathway_census() {
        let elements = vec![
            OrbitalElements {
                a: 300.0,
                e: 0.7,
                i: 20.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 50.0,
                e: 0.5,
                i: 110.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 200.0,
                e: 0.8,
                i: 80.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ];

        let (scattered, kozai, high_i) = pathway_census(&elements);
        assert_eq!(scattered, 1);
        assert_eq!(kozai, 1);
        assert_eq!(high_i, 1);
    }

    #[test]
    fn test_matches_high_i_region() {
        let drac = OrbitalElements {
            a: 41.4,
            e: 0.49,
            i: 103.4 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert!(matches_high_i_region(&drac));

        let distant = OrbitalElements {
            a: 300.0,
            e: 0.7,
            i: 100.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert!(!matches_high_i_region(&distant));
    }

    #[test]
    fn test_kozai_lidov_detection() {
        // Create a synthetic Kozai-Lidov oscillation
        let n = 50;
        let mut elements = Vec::new();
        let mut times = Vec::new();

        for j in 0..n {
            let phase = 2.0 * std::f64::consts::PI * j as f64 / n as f64;
            // Kozai: e and i oscillate anti-correlated
            // sqrt(1-e²)*cos(i) ≈ const
            // Use moderate range so e stays physical: i ∈ (10°, 50°)
            let i = 30.0 * DEG2RAD + 20.0 * DEG2RAD * phase.sin();
            let kozai_const = 0.5;
            let cos_i = i.cos();
            let e = (1.0 - (kozai_const / cos_i).powi(2))
                .sqrt()
                .clamp(0.01, 0.99);

            elements.push(OrbitalElements {
                a: 300.0,
                e,
                i,
                omega: phase,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            });
            times.push(j as f64 * 1e6 * YEAR_DAYS);
        }

        let history = ParticleHistory {
            id: 0,
            times,
            elements,
        };

        assert!(
            detect_kozai_lidov(&history),
            "Should detect Kozai-Lidov in synthetic data"
        );
    }
}
