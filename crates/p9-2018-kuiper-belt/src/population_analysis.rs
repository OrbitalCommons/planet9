//! Population classification and perihelion distribution analysis.
//!
//! Classifies surviving particles as aligned, anti-aligned, or circulating
//! from their full Δϖ(t) time series (circular mean + concentration across
//! the collected snapshots), and computes perihelion distance statistics for
//! each population.
//!
//! A single-snapshot |Δϖ| cut cannot distinguish a librator from a
//! circulator caught at a random phase — circulators pollute both bins —
//! so classification requires the time series.
//!
//! Key result from Khain+ 2018:
//! - Narrow initial q → only anti-aligned objects
//! - Broad initial q → bimodal: low-q anti-aligned + high-q aligned

use std::f64::consts::PI;

use p9_core::analysis::circular::{circular_mean, mean_resultant_length};

use crate::simulation::{delta_varpi_series, perihelion_distances, SimulationResult};

/// Classification of a particle's apsidal behavior relative to Planet Nine.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
pub enum Alignment {
    /// Δϖ(t) concentrated about 0: librating aligned with P9
    Aligned,
    /// Δϖ(t) concentrated about ±π: librating opposite P9
    AntiAligned,
    /// Δϖ(t) uniform on the circle: apsidally circulating
    Circulating,
}

/// Minimum mean resultant length R̄ of the Δϖ(t) series for a particle to
/// count as apsidally confined rather than circulating. R̄ ≈ sin(A)/A for
/// libration amplitude A, so 0.5 corresponds to A ≈ 110°.
pub const LIBRATION_R_BAR_MIN: f64 = 0.5;

/// Minimum number of snapshots a particle must survive to be classified.
pub const MIN_SERIES_LEN: usize = 4;

/// Statistics for a population of particles.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PopulationStats {
    pub count: usize,
    pub median_q: f64,
    pub mean_q: f64,
    pub min_q: f64,
    pub max_q: f64,
}

/// Result of population classification.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct PopulationResult {
    pub aligned: PopulationStats,
    pub anti_aligned: PopulationStats,
    pub circulating: PopulationStats,
    pub total_survivors: usize,
    /// Aligned fraction among apsidally confined survivors.
    pub aligned_fraction: f64,
}

/// Classify one particle's Δϖ(t) series.
///
/// Returns None when the series is too short (removed early / never
/// observed). A series with mean resultant length R̄ < `LIBRATION_R_BAR_MIN`
/// circulates; otherwise the circular mean decides aligned (|⟨Δϖ⟩| < π/2)
/// vs anti-aligned.
pub fn classify_series(series: &[f64]) -> Option<Alignment> {
    if series.len() < MIN_SERIES_LEN {
        return None;
    }
    let r_bar = mean_resultant_length(series);
    if r_bar < LIBRATION_R_BAR_MIN {
        return Some(Alignment::Circulating);
    }
    let mean = circular_mean(series)?;
    let mean = if mean > PI { mean - 2.0 * PI } else { mean };
    if mean.abs() < PI / 2.0 {
        Some(Alignment::Aligned)
    } else {
        Some(Alignment::AntiAligned)
    }
}

/// Classify every particle of a simulation from its Δϖ(t) snapshot series.
/// Index-aligned with the particle slots; None = insufficient data.
pub fn classify_alignment_series(result: &SimulationResult) -> Vec<Option<Alignment>> {
    delta_varpi_series(result)
        .iter()
        .map(|s| classify_series(s))
        .collect()
}

/// Compute population statistics for aligned / anti-aligned / circulating
/// groups, classifying from the Δϖ(t) series and taking perihelia from the
/// final snapshot.
pub fn population_statistics_series(result: &SimulationResult) -> PopulationResult {
    let alignments = classify_alignment_series(result);
    let final_snapshot = result.snapshots.last().expect("at least one snapshot");
    let qs = perihelion_distances(final_snapshot);

    let mut aligned_qs: Vec<f64> = Vec::new();
    let mut anti_aligned_qs: Vec<f64> = Vec::new();
    let mut circulating_qs: Vec<f64> = Vec::new();

    for (i, alignment) in alignments.iter().enumerate() {
        // Only surviving particles have a final-snapshot perihelion.
        let (Some(alignment), Some(q)) = (alignment, qs[i]) else {
            continue;
        };
        match alignment {
            Alignment::Aligned => aligned_qs.push(q),
            Alignment::AntiAligned => anti_aligned_qs.push(q),
            Alignment::Circulating => circulating_qs.push(q),
        }
    }

    let confined = aligned_qs.len() + anti_aligned_qs.len();
    let total = confined + circulating_qs.len();
    let aligned_frac = if confined > 0 {
        aligned_qs.len() as f64 / confined as f64
    } else {
        0.0
    };

    PopulationResult {
        aligned: compute_stats(&aligned_qs),
        anti_aligned: compute_stats(&anti_aligned_qs),
        circulating: compute_stats(&circulating_qs),
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

/// Compute perihelion distance histogram for the final snapshot.
///
/// Returns (bin_centers, counts) for bins spanning q_min to q_max.
pub fn perihelion_histogram(
    result: &SimulationResult,
    q_range: (f64, f64),
    n_bins: usize,
) -> (Vec<f64>, Vec<usize>) {
    let bin_width = (q_range.1 - q_range.0) / n_bins as f64;
    let centers: Vec<f64> = (0..n_bins)
        .map(|i| q_range.0 + (i as f64 + 0.5) * bin_width)
        .collect();
    let mut counts = vec![0usize; n_bins];

    let final_snapshot = result.snapshots.last().expect("at least one snapshot");
    for q in perihelion_distances(final_snapshot).into_iter().flatten() {
        let bin = ((q - q_range.0) / bin_width).floor() as isize;
        if bin >= 0 && bin < n_bins as isize {
            counts[bin as usize] += 1;
        }
    }

    (centers, counts)
}

/// Compute separate perihelion histograms for the aligned and anti-aligned
/// populations (series-classified; circulators excluded).
pub fn perihelion_histograms_by_alignment(
    result: &SimulationResult,
    q_range: (f64, f64),
    n_bins: usize,
) -> (Vec<f64>, Vec<usize>, Vec<usize>) {
    let bin_width = (q_range.1 - q_range.0) / n_bins as f64;
    let centers: Vec<f64> = (0..n_bins)
        .map(|i| q_range.0 + (i as f64 + 0.5) * bin_width)
        .collect();
    let mut aligned_counts = vec![0usize; n_bins];
    let mut anti_counts = vec![0usize; n_bins];

    let alignments = classify_alignment_series(result);
    let final_snapshot = result.snapshots.last().expect("at least one snapshot");
    let qs = perihelion_distances(final_snapshot);

    for (i, q) in qs.iter().enumerate() {
        let (Some(q), Some(alignment)) = (q, alignments[i]) else {
            continue;
        };
        let bin = ((q - q_range.0) / bin_width).floor() as isize;
        if bin >= 0 && bin < n_bins as isize {
            match alignment {
                Alignment::Aligned => aligned_counts[bin as usize] += 1,
                Alignment::AntiAligned => anti_counts[bin as usize] += 1,
                Alignment::Circulating => {}
            }
        }
    }

    (centers, aligned_counts, anti_counts)
}

/// Check if a population shows bimodal structure (both aligned and
/// anti-aligned confined members present).
pub fn is_bimodal(result: &PopulationResult) -> bool {
    result.aligned.count > 0 && result.anti_aligned.count > 0
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation::{KbSnapshot, SimulationResult};
    use p9_core::types::OrbitalElements;

    fn elem(varpi: f64, a: f64, e: f64) -> Option<OrbitalElements> {
        Some(OrbitalElements {
            a,
            e,
            i: 0.0,
            omega: varpi,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        })
    }

    /// Build a synthetic result with prescribed Δϖ(t) per particle
    /// (varpi_p9 = 0 so Δϖ = ϖ).
    fn make_result(series: Vec<Vec<f64>>, a: f64, e: f64) -> SimulationResult {
        let n_snaps = series.iter().map(|s| s.len()).max().unwrap_or(0);
        let n = series.len();
        let snapshots: Vec<KbSnapshot> = (0..n_snaps)
            .map(|k| {
                let elements: Vec<Option<OrbitalElements>> = (0..n)
                    .map(|i| series[i].get(k).and_then(|&v| elem(v, a, e)))
                    .collect();
                let active_count = elements.iter().flatten().count();
                KbSnapshot {
                    t: k as f64,
                    elements,
                    active_count,
                    total_count: n,
                }
            })
            .collect();
        SimulationResult {
            snapshots,
            config_summary: String::new(),
            varpi_p9: 0.0,
        }
    }

    fn librating(center: f64, amplitude: f64, n: usize) -> Vec<f64> {
        (0..n)
            .map(|k| center + amplitude * (k as f64 * 0.7).sin())
            .collect()
    }

    fn circulating(n: usize) -> Vec<f64> {
        (0..n).map(|k| k as f64 * 2.0 * PI / n as f64).collect()
    }

    #[test]
    fn test_classify_series_aligned() {
        assert_eq!(
            classify_series(&librating(0.0, 0.5, 20)),
            Some(Alignment::Aligned)
        );
    }

    #[test]
    fn test_classify_series_anti_aligned() {
        assert_eq!(
            classify_series(&librating(PI, 0.5, 20)),
            Some(Alignment::AntiAligned)
        );
        // Libration about π that wraps the branch cut still classifies.
        assert_eq!(
            classify_series(&librating(-PI, 0.8, 20)),
            Some(Alignment::AntiAligned)
        );
    }

    #[test]
    fn test_classify_series_circulating() {
        // The single-snapshot cut would have placed every one of these
        // samples in aligned or anti-aligned; the series classifier sees the
        // full circle and rejects them.
        assert_eq!(
            classify_series(&circulating(24)),
            Some(Alignment::Circulating)
        );
    }

    #[test]
    fn test_classify_series_too_short() {
        assert_eq!(classify_series(&[0.1, 0.2]), None);
    }

    #[test]
    fn test_population_statistics_series() {
        let result = make_result(
            vec![
                librating(0.0, 0.4, 12), // aligned, q = 500*0.2 = 100
                librating(PI, 0.4, 12),  // anti-aligned
                circulating(12),         // circulating
                vec![0.0, 0.1],          // removed early -> unclassified
            ],
            500.0,
            0.8,
        );
        let stats = population_statistics_series(&result);
        assert_eq!(stats.aligned.count, 1);
        assert_eq!(stats.anti_aligned.count, 1);
        assert_eq!(stats.circulating.count, 1);
        assert!(is_bimodal(&stats));
        assert!((stats.aligned_fraction - 0.5).abs() < 1e-12);
        assert!((stats.aligned.median_q - 100.0).abs() < 1e-9);
    }

    #[test]
    fn test_anti_aligned_only() {
        let result = make_result(
            vec![librating(PI, 0.3, 10), librating(-PI, 0.5, 10)],
            300.0,
            0.9,
        );
        let stats = population_statistics_series(&result);
        assert_eq!(stats.aligned.count, 0);
        assert_eq!(stats.anti_aligned.count, 2);
        assert!(!is_bimodal(&stats));
    }

    #[test]
    fn test_perihelion_histogram() {
        let result = make_result(
            vec![
                librating(0.0, 0.4, 8),
                librating(PI, 0.4, 8),
                circulating(8),
            ],
            300.0,
            0.9, // q = 30
        );
        let (centers, counts) = perihelion_histogram(&result, (0.0, 300.0), 6);
        assert_eq!(centers.len(), 6);
        let total: usize = counts.iter().sum();
        assert_eq!(total, 3);

        let (_, aligned, anti) = perihelion_histograms_by_alignment(&result, (0.0, 300.0), 6);
        assert_eq!(aligned.iter().sum::<usize>(), 1);
        assert_eq!(anti.iter().sum::<usize>(), 1);
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
