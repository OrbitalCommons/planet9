//! Run the isotropy suite across the paper's Cases and angles, under a chosen
//! null (flat-isotropic or DES-selection-folded).

use crate::des_sample::{angle_values, case_sample, Angle, DesEtno, SampleCase};
use crate::selection::{selection_aware_p, NullModel};
use crate::statistics::{beran_an, watson_u2};
use p9_core::analysis::circular::{kuiper_statistic, mean_resultant_length};

/// Per-(angle) result: the four statistics' p-values under the chosen null,
/// plus the observed R̄.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct AnglePValues {
    pub n: usize,
    pub r_bar: f64,
    pub rayleigh_p: f64,
    pub kuiper_p: f64,
    pub watson_p: f64,
    pub beran_p: f64,
}

impl AnglePValues {
    /// The four hypothesis-test p-values.
    pub fn p_values(&self) -> [f64; 4] {
        [self.rayleigh_p, self.kuiper_p, self.watson_p, self.beran_p]
    }

    /// Number of the four tests with p > alpha (consistent with isotropy).
    pub fn n_consistent(&self, alpha: f64) -> usize {
        self.p_values().iter().filter(|&&p| p > alpha).count()
    }
}

/// First-harmonic Rayleigh statistic n·R̄² (monotone in R̄); used as the
/// statistic in the selection-aware Monte Carlo so its p-value is computed
/// against the *same* null as the others (the analytic Rayleigh assumes a flat
/// null, which is exactly what the paper argues against).
fn rayleigh_z(angles: &[f64]) -> f64 {
    let n = angles.len() as f64;
    let r = mean_resultant_length(angles);
    n * r * r
}

/// Compute all four selection-aware p-values for one angle of one sample.
pub fn angle_p_values(
    sample: &[DesEtno],
    angle: Angle,
    null: NullModel,
    seed: u64,
    iters: usize,
) -> AnglePValues {
    let vals = angle_values(sample, angle);
    AnglePValues {
        n: sample.len(),
        r_bar: mean_resultant_length(&vals),
        rayleigh_p: selection_aware_p(sample, angle, null, seed, iters, rayleigh_z),
        kuiper_p: selection_aware_p(
            sample,
            angle,
            null,
            seed.wrapping_add(1),
            iters,
            kuiper_statistic,
        ),
        watson_p: selection_aware_p(sample, angle, null, seed.wrapping_add(2), iters, watson_u2),
        beran_p: selection_aware_p(sample, angle, null, seed.wrapping_add(3), iters, beran_an),
    }
}

/// One (Case, angle) cell of the isotropy battery.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BatteryCell {
    pub case: String,
    pub angle: String,
    pub p: AnglePValues,
}

/// Run the full battery: every angle (Ω, ω, ϖ) in every Case, under `null`.
pub fn run_battery(null: NullModel, seed: u64, mc_iters: usize) -> Vec<BatteryCell> {
    let mut cells = Vec::new();
    for (ci, case) in SampleCase::all().iter().enumerate() {
        let sample = case_sample(*case);
        for (ai, angle) in Angle::all().iter().enumerate() {
            let cell_seed = seed
                .wrapping_add((ci as u64) * 1000)
                .wrapping_add((ai as u64) * 10);
            let p = angle_p_values(&sample, *angle, null, cell_seed, mc_iters);
            cells.push(BatteryCell {
                case: format!("{case:?}"),
                angle: angle.label().to_string(),
                p,
            });
        }
    }
    cells
}

/// Summary of the battery: how many of all (cells × tests) p-values fall below
/// alpha, and the smallest p-value seen.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BatterySummary {
    pub n_cells: usize,
    pub n_tests: usize,
    pub n_significant: usize,
    pub min_p: f64,
    pub fraction_consistent: f64,
}

pub fn summarize(cells: &[BatteryCell], alpha: f64) -> BatterySummary {
    let mut n_tests = 0usize;
    let mut n_significant = 0usize;
    let mut min_p = 1.0f64;
    for cell in cells {
        for p in cell.p.p_values() {
            n_tests += 1;
            if p < alpha {
                n_significant += 1;
            }
            min_p = min_p.min(p);
        }
    }
    BatterySummary {
        n_cells: cells.len(),
        n_tests,
        n_significant,
        min_p,
        fraction_consistent: 1.0 - n_significant as f64 / n_tests as f64,
    }
}

/// The paper's qualitative reference conclusion, kept as a labelled string.
pub const PAPER_CONCLUSION: &str =
    "The DES data taken on their own are consistent with azimuthal isotropy and \
     do not require a 'Planet 9' hypothesis (Bernardinelli et al. 2020).";

/// The paper applied Kuiper's test to 3 angles in 4 Cases = 12 tests; exactly
/// 2 of them reached p = 0.03 (both for the Ω distribution, a > 250 AU). The
/// remaining 10 of 12 were not significant — the basis of the "consistent with
/// isotropy" conclusion.
pub const PAPER_KUIPER_TOTAL_TESTS: usize = 12;
pub const PAPER_KUIPER_SIGNIFICANT_TESTS: usize = 2;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_battery_has_12_cells() {
        let cells = run_battery(NullModel::DesSelection, 2020, 1000);
        assert_eq!(cells.len(), 12, "4 Cases x 3 angles");
    }

    #[test]
    fn test_paper_reference_constants() {
        assert_eq!(PAPER_KUIPER_TOTAL_TESTS, 12);
        assert_eq!(PAPER_KUIPER_SIGNIFICANT_TESTS, 2);
        assert!(PAPER_CONCLUSION.contains("isotropy"));
    }
}
