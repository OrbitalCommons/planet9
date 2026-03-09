//! Population classification and perihelion distribution analysis.
//!
//! Classifies surviving particles as aligned or anti-aligned with Planet Nine
//! based on Δϖ, and computes perihelion distance statistics for each population.
//!
//! Key result from Khain+ 2018:
//! - Narrow initial q → only anti-aligned survivors
//! - Broad initial q → bimodal: low-q anti-aligned + high-q aligned

use std::f64::consts::PI;

use crate::simulation::{delta_varpi_values, perihelion_distances, KbSnapshot};

/// Classification of a particle's alignment with Planet Nine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Alignment {
    /// |Δϖ| < π/2: perihelion longitude near P9's
    Aligned,
    /// |Δϖ| > π/2: perihelion longitude opposite P9's
    AntiAligned,
}

/// Statistics for a population of particles.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PopulationStats {
    pub count: usize,
    pub median_q: f64,
    pub mean_q: f64,
    pub min_q: f64,
    pub max_q: f64,
}

/// Result of population classification for a single snapshot.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PopulationResult {
    pub aligned: PopulationStats,
    pub anti_aligned: PopulationStats,
    pub total_survivors: usize,
    pub aligned_fraction: f64,
}

/// Classify particles as aligned or anti-aligned with Planet Nine.
///
/// `varpi_p9` is Planet Nine's longitude of perihelion (ω₉ + Ω₉).
/// Objects with |Δϖ| < π/2 are aligned; |Δϖ| > π/2 are anti-aligned.
pub fn classify_alignment(snapshot: &KbSnapshot, varpi_p9: f64) -> Vec<Alignment> {
    let dvs = delta_varpi_values(snapshot, varpi_p9);
    dvs.iter()
        .map(|&dv| {
            if dv.abs() < PI / 2.0 {
                Alignment::Aligned
            } else {
                Alignment::AntiAligned
            }
        })
        .collect()
}

/// Compute population statistics for aligned and anti-aligned groups.
pub fn population_statistics(snapshot: &KbSnapshot, varpi_p9: f64) -> PopulationResult {
    let alignments = classify_alignment(snapshot, varpi_p9);
    let qs = perihelion_distances(snapshot);

    let mut aligned_qs: Vec<f64> = Vec::new();
    let mut anti_aligned_qs: Vec<f64> = Vec::new();

    for (i, &alignment) in alignments.iter().enumerate() {
        match alignment {
            Alignment::Aligned => aligned_qs.push(qs[i]),
            Alignment::AntiAligned => anti_aligned_qs.push(qs[i]),
        }
    }

    let total = aligned_qs.len() + anti_aligned_qs.len();
    let aligned_frac = if total > 0 {
        aligned_qs.len() as f64 / total as f64
    } else {
        0.0
    };

    PopulationResult {
        aligned: compute_stats(&aligned_qs),
        anti_aligned: compute_stats(&anti_aligned_qs),
        total_survivors: total,
        aligned_fraction: aligned_frac,
    }
}

fn compute_stats(qs: &[f64]) -> PopulationStats {
    if qs.is_empty() {
        return PopulationStats {
            count: 0,
            median_q: 0.0,
            mean_q: 0.0,
            min_q: 0.0,
            max_q: 0.0,
        };
    }

    let mut sorted = qs.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());

    let median = if sorted.len() % 2 == 0 {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    } else {
        sorted[sorted.len() / 2]
    };

    PopulationStats {
        count: qs.len(),
        median_q: median,
        mean_q: qs.iter().sum::<f64>() / qs.len() as f64,
        min_q: sorted[0],
        max_q: *sorted.last().unwrap(),
    }
}

/// Compute perihelion distance histogram.
///
/// Returns (bin_centers, counts) for bins spanning q_min to q_max.
pub fn perihelion_histogram(
    snapshot: &KbSnapshot,
    q_range: (f64, f64),
    n_bins: usize,
) -> (Vec<f64>, Vec<usize>) {
    let bin_width = (q_range.1 - q_range.0) / n_bins as f64;
    let centers: Vec<f64> = (0..n_bins)
        .map(|i| q_range.0 + (i as f64 + 0.5) * bin_width)
        .collect();
    let mut counts = vec![0usize; n_bins];

    let qs = perihelion_distances(snapshot);
    for q in qs {
        let bin = ((q - q_range.0) / bin_width).floor() as isize;
        if bin >= 0 && bin < n_bins as isize {
            counts[bin as usize] += 1;
        }
    }

    (centers, counts)
}

/// Compute separate perihelion histograms for aligned and anti-aligned populations.
pub fn perihelion_histograms_by_alignment(
    snapshot: &KbSnapshot,
    varpi_p9: f64,
    q_range: (f64, f64),
    n_bins: usize,
) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
    let bin_width = (q_range.1 - q_range.0) / n_bins as f64;
    let centers: Vec<f64> = (0..n_bins)
        .map(|i| q_range.0 + (i as f64 + 0.5) * bin_width)
        .collect();
    let mut aligned_counts = vec![0usize; n_bins];
    let mut anti_counts = vec![0usize; n_bins];

    let alignments = classify_alignment(snapshot, varpi_p9);
    let qs = perihelion_distances(snapshot);

    for (i, q) in qs.iter().enumerate() {
        let bin = ((q - q_range.0) / bin_width).floor() as isize;
        if bin >= 0 && bin < n_bins as isize {
            match alignments[i] {
                Alignment::Aligned => aligned_counts[bin as usize] += 1,
                Alignment::AntiAligned => anti_counts[bin as usize] += 1,
            }
        }
    }

    (centers, aligned_counts, anti_counts)
}

/// Check if a population shows bimodal structure (both aligned and anti-aligned present).
pub fn is_bimodal(result: &PopulationResult) -> bool {
    result.aligned.count > 0 && result.anti_aligned.count > 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::types::OrbitalElements;

    fn make_test_snapshot(elements: Vec<OrbitalElements>) -> KbSnapshot {
        let n = elements.len();
        KbSnapshot {
            t: 0.0,
            elements,
            active_count: n,
            total_count: n,
        }
    }

    #[test]
    fn test_classify_aligned() {
        // P9 at ϖ = π (180°)
        let varpi_p9 = PI;

        let snap = make_test_snapshot(vec![
            // ϖ = π → Δϖ = 0 → aligned
            OrbitalElements {
                a: 300.0,
                e: 0.8,
                i: 0.0,
                omega: PI / 2.0,
                omega_big: PI / 2.0,
                mean_anomaly: 0.0,
            },
            // ϖ = 0 → Δϖ = -π → anti-aligned
            OrbitalElements {
                a: 200.0,
                e: 0.5,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ]);

        let alignments = classify_alignment(&snap, varpi_p9);
        assert_eq!(alignments[0], Alignment::Aligned);
        assert_eq!(alignments[1], Alignment::AntiAligned);
    }

    #[test]
    fn test_population_statistics() {
        let varpi_p9 = PI;
        let snap = make_test_snapshot(vec![
            // Aligned, q = 100 AU
            OrbitalElements {
                a: 500.0,
                e: 0.8,
                i: 0.0,
                omega: PI,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            // Anti-aligned, q = 35 AU
            OrbitalElements {
                a: 350.0,
                e: 0.9,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ]);

        let result = population_statistics(&snap, varpi_p9);
        assert_eq!(result.aligned.count, 1);
        assert_eq!(result.anti_aligned.count, 1);
        assert!(is_bimodal(&result));
        assert!((result.aligned.median_q - 100.0).abs() < 1e-10);
        assert!((result.anti_aligned.median_q - 35.0).abs() < 1e-10);
    }

    #[test]
    fn test_anti_aligned_only() {
        let varpi_p9 = PI;
        let snap = make_test_snapshot(vec![
            OrbitalElements {
                a: 300.0,
                e: 0.9,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: 250.0,
                e: 0.85,
                i: 0.0,
                omega: 0.3,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ]);

        let result = population_statistics(&snap, varpi_p9);
        assert_eq!(result.aligned.count, 0);
        assert_eq!(result.anti_aligned.count, 2);
        assert!(!is_bimodal(&result));
    }

    #[test]
    fn test_perihelion_histogram() {
        let snap = make_test_snapshot(vec![
            OrbitalElements {
                a: 300.0,
                e: 0.9,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            }, // q=30
            OrbitalElements {
                a: 200.0,
                e: 0.5,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            }, // q=100
            OrbitalElements {
                a: 500.0,
                e: 0.6,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            }, // q=200
        ]);

        let (centers, counts) = perihelion_histogram(&snap, (0.0, 300.0), 6);
        assert_eq!(centers.len(), 6);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 3);
    }

    #[test]
    fn test_empty_population_stats() {
        let stats = compute_stats(&[]);
        assert_eq!(stats.count, 0);
        assert_eq!(stats.median_q, 0.0);
    }

    #[test]
    fn test_median_computation() {
        // Odd count
        let stats = compute_stats(&[10.0, 30.0, 20.0]);
        assert!((stats.median_q - 20.0).abs() < 1e-10);

        // Even count
        let stats = compute_stats(&[10.0, 20.0, 30.0, 40.0]);
        assert!((stats.median_q - 25.0).abs() < 1e-10);
    }
}
