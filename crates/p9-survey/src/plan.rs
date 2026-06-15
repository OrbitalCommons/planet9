//! Observing-plan optimizer: given a fixed time budget on JBT/SPENCER, how much
//! of the Planet Nine probability can you actually cover, and how?
//!
//! The key asymmetry for this instrument: **depth is cheap, area is dear.** A
//! single 300 s visit already reaches V≈23.8 while Planet Nine is V≈19–23, so
//! the binding constraint is the 0.32 deg² field — you cannot tile thousands of
//! square degrees in days. The optimizer therefore (1) converts a time budget
//! into a coverable *area* and a stacked *depth*, and (2) greedily covers the
//! highest-probability sky cells up to that area, returning the captured
//! probability — the detection chance, assuming the stacked depth reaches the
//! source and the multi-epoch cadence links the motion.

use crate::schema::{SkyGrid, StudyResult};
use p9_core::constants::DEG2RAD;

/// A time budget and how it is spent per field.
#[derive(Debug, Clone)]
pub struct TimeBudget {
    /// Total wall-clock time on the telescope (hours).
    pub hours_total: f64,
    /// Fraction of wall-clock actually integrating on target (Earth
    /// occultation from LEO, SAA, slew/settle between fields, downlink).
    pub duty_cycle: f64,
    /// Open-shutter seconds per single visit.
    pub exposure_s: f64,
    /// Dead time per visit (readout + slew + settle).
    pub overhead_s: f64,
    /// Revisits of each field across the run, for moving-object linking.
    pub epochs: u32,
    /// Instantaneous field of view (deg²).
    pub fov_deg2: f64,
}

impl TimeBudget {
    pub fn useful_seconds(&self) -> f64 {
        self.hours_total * 3600.0 * self.duty_cycle
    }
    /// Total time invested per field (all epochs).
    pub fn seconds_per_field(&self) -> f64 {
        self.epochs as f64 * (self.exposure_s + self.overhead_s)
    }
    pub fn n_fields(&self) -> f64 {
        self.useful_seconds() / self.seconds_per_field()
    }
    /// Sky area coverable with this budget (deg²).
    pub fn area_deg2(&self) -> f64 {
        self.n_fields() * self.fov_deg2
    }
    /// Total open-shutter integration co-added per field (all epochs).
    pub fn integration_per_field_s(&self) -> f64 {
        self.epochs as f64 * self.exposure_s
    }
}

/// Background-limited limiting magnitude scaling from a reference depth: a
/// √(t) sensitivity gain → Δm = 1.25·log10(t/t_ref). Anchored at the JBT
/// single-visit depth (V≈23.8 in 300 s; see telescope.rs).
pub const JBT_REF_DEPTH_300S: f64 = 23.8;
pub fn stacked_depth(total_integration_s: f64) -> f64 {
    JBT_REF_DEPTH_300S + 1.25 * (total_integration_s / 300.0).log10()
}

/// Result of greedily covering the most-probable cells up to an area budget.
#[derive(Debug, Clone)]
pub struct Coverage {
    /// Position probability enclosed by the covered cells (0..1).
    pub captured_prob: f64,
    pub area_used_deg2: f64,
    pub n_cells: usize,
    /// Bounding box and probability-weighted centroid of the covered region.
    pub ra_lo_deg: f64,
    pub ra_hi_deg: f64,
    pub dec_lo_deg: f64,
    pub dec_hi_deg: f64,
    pub ra_centroid_deg: f64,
    pub dec_centroid_deg: f64,
}

/// Greedily select the highest-probability grid cells until `area_budget_deg2`
/// of (cos Dec-weighted) solid angle is used, returning the enclosed
/// probability and the region geometry.
pub fn capture(prob: &[f64], grid: &SkyGrid, area_budget_deg2: f64) -> Coverage {
    let cell_base = grid.dra() * grid.ddec();
    let mut cells: Vec<(usize, f64, f64)> = (0..prob.len())
        .map(|idx| {
            let dec = grid.dec_center(idx / grid.n_ra);
            let solid = cell_base * (dec * DEG2RAD).cos();
            (idx, prob[idx], solid)
        })
        .collect();
    cells.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap());

    let (mut captured, mut area) = (0.0, 0.0);
    let (mut n, mut sx, mut sy, mut sw) = (0usize, 0.0, 0.0, 0.0);
    let (mut ra_lo, mut ra_hi) = (f64::MAX, f64::MIN);
    let (mut dec_lo, mut dec_hi) = (f64::MAX, f64::MIN);
    for (idx, p, solid) in cells {
        if area + solid > area_budget_deg2 && n > 0 {
            break;
        }
        let ra = grid.ra_min_deg + (idx % grid.n_ra) as f64 * grid.dra() + grid.dra() / 2.0;
        let dec = grid.dec_center(idx / grid.n_ra);
        captured += p;
        area += solid;
        n += 1;
        sx += p * ra;
        sy += p * dec;
        sw += p;
        ra_lo = ra_lo.min(ra);
        ra_hi = ra_hi.max(ra);
        dec_lo = dec_lo.min(dec);
        dec_hi = dec_hi.max(dec);
    }
    Coverage {
        captured_prob: captured,
        area_used_deg2: area,
        n_cells: n,
        ra_lo_deg: ra_lo,
        ra_hi_deg: ra_hi,
        dec_lo_deg: dec_lo,
        dec_hi_deg: dec_hi,
        ra_centroid_deg: if sw > 0.0 { sx / sw } else { 0.0 },
        dec_centroid_deg: if sw > 0.0 { sy / sw } else { 0.0 },
    }
}

/// Equal-weight blend of the studies' probability maps (a prior that does not
/// bet on one orbit solution), renormalized to sum 1.
pub fn blended_prob(studies: &[StudyResult]) -> Vec<f64> {
    let n = studies.first().map(|s| s.prob.len()).unwrap_or(0);
    let mut blend = vec![0.0; n];
    for s in studies {
        for (b, p) in blend.iter_mut().zip(&s.prob) {
            *b += p;
        }
    }
    let total: f64 = blend.iter().sum();
    if total > 0.0 {
        for b in blend.iter_mut() {
            *b /= total;
        }
    }
    blend
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{default_grid, run_survey};

    fn budget() -> TimeBudget {
        TimeBudget {
            hours_total: 96.0,
            duty_cycle: 0.5,
            exposure_s: 120.0,
            overhead_s: 60.0,
            epochs: 3,
            fov_deg2: 0.32,
        }
    }

    #[test]
    fn budget_area_is_order_100_deg2() {
        let b = budget();
        // 96 h × 0.5 = 172,800 s; per field 3×180 = 540 s → 320 fields → ~102 deg².
        assert!(
            (b.area_deg2() - 102.0).abs() < 5.0,
            "area = {}",
            b.area_deg2()
        );
    }

    #[test]
    fn stacked_depth_increases_with_time() {
        assert!((stacked_depth(300.0) - JBT_REF_DEPTH_300S).abs() < 1e-9);
        assert!(stacked_depth(3600.0) > stacked_depth(300.0));
        // An hour co-added reaches ~V25.1.
        assert!((stacked_depth(3600.0) - 25.1).abs() < 0.2);
    }

    #[test]
    fn capture_is_monotone_and_bounded() {
        let d = run_survey(40_000, 7);
        let grid = default_grid();
        let s = &d.studies[3]; // 2021 MCMC
        let small = capture(&s.prob, &grid, 50.0).captured_prob;
        let big = capture(&s.prob, &grid, 500.0).captured_prob;
        assert!(small <= big);
        assert!((0.0..=1.0).contains(&big));
        // The peaked aphelion lobe means even ~100 deg² captures a real chunk.
        let c = capture(&s.prob, &grid, 100.0);
        assert!(c.captured_prob > 0.02, "captured = {}", c.captured_prob);
        // Region lands in the expected RA/Dec lobe.
        assert!((40.0..130.0).contains(&c.ra_centroid_deg));
    }

    #[test]
    fn concentrated_solution_captures_more_per_area() {
        let d = run_survey(60_000, 11);
        let grid = default_grid();
        let siraj = d
            .studies
            .iter()
            .find(|s| s.solution.name.contains("Siraj"))
            .unwrap();
        let b2016 = d
            .studies
            .iter()
            .find(|s| s.solution.name.contains("2016 Batygin & Brown (nominal)"))
            .unwrap();
        // The tighter 2024 solution yields more captured probability per 100 deg²
        // than the diffuse 2016 one.
        let cs = capture(&siraj.prob, &grid, 100.0).captured_prob;
        let cb = capture(&b2016.prob, &grid, 100.0).captured_prob;
        assert!(cs > cb, "Siraj {cs} should beat 2016 {cb}");
    }
}
