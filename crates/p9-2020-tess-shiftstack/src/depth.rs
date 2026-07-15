//! Shift-stack depth gain for TESS, and the Planet Nine / TNO reach comparison.
//!
//! The depth machinery is *not* re-derived here: it is the matched-filter √N
//! law already implemented and pinned in `p9-2025-stacking`. Co-adding N
//! background-limited frames along the correct track deepens the limiting
//! magnitude by Δm = 1.25·log₁₀N, so TESS's shallow single-frame depth
//! (T ≈ 15) becomes a deep effective depth once thousands of FFIs are stacked.
//! This module just feeds the TESS frame counts ([`crate::tess`]) into that law
//! and compares the result against bodies of interest using `p9-core`
//! photometry.
//!
//! Rice & Laughlin recover Sedna, 2015 BP519 and 2007 TG422 — all in the
//! V ≈ 20.6–22.3 range — using ONE- to TWO-sector stacks, and report a search
//! sensitivity of V < 21 for distances d ≲ 150 au. The depth budget therefore
//! starts from the search's EFFECTIVE per-frame depth
//! (`published::TESS_SINGLE_DEPTH_EFFECTIVE` ≈ 17.15, back-derived from that
//! V < 21 baseline sensitivity), not the naive 5σ FFI catalog depth T ≈ 15 —
//! which would misrepresent the paper's own recoveries as 20–600-sector
//! co-adds (~3 mag of missing budget). The faintest recovery (TG422,
//! V = 22.3) additionally rides its red V−T color (~1 mag for these TNOs):
//! in the detection band it sits within the two-sector reach.

use p9_core::analysis::photometry::planet_apparent_magnitude;
use p9_core::analysis::stacking::matched_filter::{
    frames_to_reach_depth, stack_depth_gain_mag, stacked_limiting_mag,
};
use p9_core::types::P9Params;

use crate::tess;

/// Stacked TESS limiting magnitude after co-adding `n_frames` FFIs, starting
/// from a single-frame depth (T magnitudes). Thin wrapper over the
/// p9-2025-stacking √N law so the depth model is shared, not duplicated.
pub fn stacked_depth(single_frame_depth: f64, n_frames: f64) -> f64 {
    stacked_limiting_mag(single_frame_depth, n_frames.round() as usize)
}

/// Stacked TESS limiting magnitude reached over `n_sectors` observing sectors,
/// from a single-frame depth.
pub fn stacked_depth_over_sectors(single_frame_depth: f64, n_sectors: f64) -> f64 {
    stacked_depth(single_frame_depth, tess::frames_in_sectors(n_sectors))
}

/// The pure stacking gain (magnitudes deeper than one frame) for `n_frames`.
pub fn depth_gain(n_frames: f64) -> f64 {
    stack_depth_gain_mag(n_frames.round() as usize)
}

/// How many sectors of TESS FFIs must be stacked to reach `target_depth` from
/// `single_frame_depth`. Inverts the √N law and divides by the per-sector
/// frame count.
pub fn sectors_to_reach_depth(single_frame_depth: f64, target_depth: f64) -> f64 {
    frames_to_reach_depth(single_frame_depth, target_depth) / tess::frames_per_sector()
}

/// Apparent V magnitude at opposition of a body of P9's mass/albedo sitting at
/// heliocentric distance `r_au`, from p9-core's reflected-sunlight photometry
/// (mass→radius→H→m). Albedo defaults to a Neptune-like 0.4 if non-positive.
pub fn p9_apparent_magnitude_at(p9: &P9Params, albedo: f64, r_au: f64) -> f64 {
    let alb = if albedo > 0.0 { albedo } else { 0.4 };
    planet_apparent_magnitude(p9.mass_earth, alb, r_au)
}

/// Nominal Planet Nine apparent V magnitude at its semi-major-axis distance —
/// the "typical" brightness stand-in. Rice & Laughlin's distance grid spans
/// 70–800 au; this evaluates at `a`, with [`p9_apparent_magnitude_at`] giving
/// the perihelion-to-aphelion brightness range.
pub fn p9_apparent_magnitude(p9: &P9Params, albedo: f64) -> f64 {
    p9_apparent_magnitude_at(p9, albedo, p9.a)
}

/// Is a body of apparent magnitude `m_body` within the stacked TESS reach
/// `m_stack`? (Brighter — smaller magnitude — than the limit is detectable.)
pub fn within_reach(m_body: f64, m_stack: f64) -> bool {
    m_body <= m_stack
}

/// Signed depth margin (mag): positive when the stacked depth is *deeper* than
/// the body (detectable headroom), negative when the body is too faint.
pub fn reach_margin(m_body: f64, m_stack: f64) -> f64 {
    m_stack - m_body
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::published;
    use approx::assert_relative_eq;

    #[test]
    fn single_sector_stack_reproduces_baseline_sensitivity() {
        // ~1.18e3 FFIs over one sector buys 1.25·log10(N) ≈ 3.85 mag: from
        // the effective per-frame depth this reproduces the paper's V < 21
        // baseline sensitivity (the calibration's defining property).
        let depth = stacked_depth_over_sectors(published::TESS_SINGLE_DEPTH_EFFECTIVE, 1.0);
        assert!(
            (depth - published::SENSITIVITY_V_LIMIT).abs() < 0.1,
            "single-sector stacked depth = {depth:.2}"
        );
        // The naive 5σ catalog depth lands 2+ mag short of the paper's own
        // demonstrated sensitivity — the ~3 mag contradiction this
        // calibration removes.
        let naive = stacked_depth_over_sectors(published::TESS_SINGLE_DEPTH, 1.0);
        assert!(
            naive < published::SENSITIVITY_V_LIMIT - 1.5,
            "naive = {naive:.2}"
        );
    }

    #[test]
    fn year_long_stack_deepens_by_t_squared_growth() {
        // A full year (13 sectors, ~1.5e4 FFIs) adds 1.25·log10(13) ≈ 1.4 mag
        // over the single-sector baseline.
        let depth = stacked_depth_over_sectors(
            published::TESS_SINGLE_DEPTH_EFFECTIVE,
            tess::SECTORS_PER_YEAR,
        );
        assert!(
            (22.0..22.8).contains(&depth),
            "one-year stacked depth = {depth:.2}"
        );
    }

    #[test]
    fn stacking_gain_follows_sqrt_n_law() {
        // 4× the frames → exactly +1.25·log10(4) ≈ 0.753 mag deeper, and the
        // gain is the shared p9-2025-stacking quantity (no local re-derivation).
        let g1 = depth_gain(1000.0);
        let g4 = depth_gain(4000.0);
        assert_relative_eq!(g4 - g1, 1.25 * 4.0_f64.log10(), epsilon = 1e-9);
    }

    #[test]
    fn recovered_tnos_within_one_to_few_sector_stacks() {
        // With the effective per-frame depth, the paper's recoveries are
        // consistent with their actual 1-2-sector stacks: Sedna (V 20.64)
        // needs about half a sector-equivalent, BP519 ~4, and TG422 ~12 on
        // the V scale — or ~1-2 sectors in the detection band once its red
        // V-T ~ 1 mag color is credited (documented in the module header).
        // The old naive-depth model demanded 20-600 sector-equivalents,
        // contradicting the paper it reproduces.
        let sedna =
            sectors_to_reach_depth(published::TESS_SINGLE_DEPTH_EFFECTIVE, published::V_SEDNA);
        assert!(
            (0.2..2.0).contains(&sedna),
            "Sedna sector-equiv = {sedna:.2}"
        );
        // Monotonically more integration for fainter targets, all bounded.
        let mut last = 0.0;
        for &v in &[published::V_SEDNA, published::V_BP519, published::V_TG422] {
            let n = sectors_to_reach_depth(published::TESS_SINGLE_DEPTH_EFFECTIVE, v);
            assert!(n > last, "non-monotone integration at V={v}: {n:.1}");
            last = n;
        }
        assert!(last < 15.0, "faintest (V-scale) sector-equiv = {last:.1}");
    }

    #[test]
    fn nominal_planet_nine_within_depth_but_outside_searched_distances() {
        // With the paper-calibrated effective depth, a year-long stack
        // (~22.4) reaches the nominal P9 across its orbit (V ~ 19.5 at a,
        // ~20.7 at 500 au) — DEPTH is not the limiter. What excludes P9 from
        // the actual search is the covered distance range: the paper's
        // V < 21 sensitivity applies to trial tracks at d ≲ 150 au, far
        // inside P9's 300-500 au. (The old naive-depth version called P9
        // "marginal" at aphelion, an artifact of the ~3 mag depth deficit.)
        let p9 = P9Params::mcmc_2021();
        let m_at_a = p9_apparent_magnitude_at(&p9, 0.4, p9.a);
        let m_aphelion = p9_apparent_magnitude_at(&p9, 0.4, 500.0);
        assert!((19.0..20.0).contains(&m_at_a), "V_P9(a) = {m_at_a:.2}");
        assert!(m_aphelion > m_at_a, "aphelion must be fainter");

        let m_stack = stacked_depth_over_sectors(
            published::TESS_SINGLE_DEPTH_EFFECTIVE,
            tess::SECTORS_PER_YEAR,
        );
        assert!((22.0..22.8).contains(&m_stack), "year depth = {m_stack:.2}");
        assert!(within_reach(m_at_a, m_stack));
        assert!(within_reach(m_aphelion, m_stack));
        assert!(reach_margin(m_aphelion, m_stack) > 1.0);

        // The searched track grid stops far inside P9's orbit.
        let q9 = p9.a * (1.0 - p9.e);
        assert!(
            q9 > 1.5 * published::SENSITIVITY_DISTANCE_AU,
            "P9 perihelion {q9:.0} au vs searched d <= {} au",
            published::SENSITIVITY_DISTANCE_AU
        );
    }
}
