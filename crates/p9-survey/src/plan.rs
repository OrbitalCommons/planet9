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
/// Smallest arc [lo, hi] (degrees; hi may exceed 360 when the arc wraps
/// through 0) containing all the given RAs: sort, find the largest gap
/// between consecutive values on the circle, and take its complement.
fn smallest_ra_arc(ras: &mut [f64]) -> (f64, f64) {
    if ras.is_empty() {
        return (f64::MAX, f64::MIN);
    }
    ras.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = ras.len();
    let (mut gap_max, mut gap_at) = (360.0 - (ras[n - 1] - ras[0]), n - 1);
    for k in 0..n - 1 {
        let gap = ras[k + 1] - ras[k];
        if gap > gap_max {
            gap_max = gap;
            gap_at = k;
        }
    }
    if gap_at == n - 1 {
        (ras[0], ras[n - 1])
    } else {
        // Arc wraps through 0: starts after the gap, ends at gap start +360.
        (ras[gap_at + 1], ras[gap_at] + 360.0)
    }
}

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
    // RA is circular: centroid via the probability-weighted vector mean,
    // bounds via the smallest arc containing all covered cells (a lobe
    // straddling RA = 0 must not report lo = 0 / hi = 360 with an antipodal
    // centroid).
    let (mut n, mut sc, mut ss, mut sy, mut sw) = (0usize, 0.0, 0.0, 0.0, 0.0);
    let mut ras: Vec<f64> = Vec::new();
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
        let ra_rad = ra.to_radians();
        sc += p * ra_rad.cos();
        ss += p * ra_rad.sin();
        sy += p * dec;
        sw += p;
        ras.push(ra);
        dec_lo = dec_lo.min(dec);
        dec_hi = dec_hi.max(dec);
    }
    let (ra_lo, ra_hi) = smallest_ra_arc(&mut ras);
    Coverage {
        captured_prob: captured,
        area_used_deg2: area,
        n_cells: n,
        ra_lo_deg: ra_lo,
        ra_hi_deg: ra_hi,
        dec_lo_deg: dec_lo,
        dec_hi_deg: dec_hi,
        ra_centroid_deg: if sw > 0.0 {
            ss.atan2(sc).to_degrees().rem_euclid(360.0)
        } else {
            0.0
        },
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

    #[test]
    fn wrapped_lobe_reports_circular_ra_stats() {
        // A lobe straddling RA = 0 (cells near 355° and 5°) must report a
        // tight wrapped arc and a centroid near 0° — not lo≈0/hi≈360 with an
        // antipodal centroid.
        let grid = SkyGrid {
            ra_min_deg: 0.0,
            ra_max_deg: 360.0,
            n_ra: 72,
            dec_min_deg: -60.0,
            dec_max_deg: 60.0,
            n_dec: 24,
        };
        let mut prob = vec![0.0; grid.len()];
        for ra in [352.5, 357.5, 2.5, 7.5] {
            let idx = grid.index(ra, 0.0).unwrap();
            prob[idx] = 0.25;
        }
        // Budget sized to the four hot cells (5 deg cells, ~25 deg^2 each).
        let cov = capture(&prob, &grid, 110.0);
        let centroid = cov.ra_centroid_deg;
        let d0 = (centroid - 0.0)
            .rem_euclid(360.0)
            .min((0.0 - centroid).rem_euclid(360.0));
        assert!(d0 < 5.0, "centroid {centroid:.1} should sit near RA 0");
        let span = cov.ra_hi_deg - cov.ra_lo_deg;
        assert!(
            (10.0..30.0).contains(&span),
            "wrapped arc span {span:.1} (lo {:.1}, hi {:.1})",
            cov.ra_lo_deg,
            cov.ra_hi_deg
        );
    }
}
