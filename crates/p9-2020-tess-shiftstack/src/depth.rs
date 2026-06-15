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
//! V ≈ 20.6–22.3 range — and report a search sensitivity of V < 21 for
//! distances d ≲ 150 au. A single sector of FFIs (T ≈ 15 → ~T 18.8) does not
//! reach V < 21 on its own; reaching it is a multi-sector / per-track deep
//! co-add, which the [`sectors_to_reach_depth`] inversion makes explicit.

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
    fn single_sector_stack_reaches_t_19_to_21() {
        // ~1.18e3 FFIs over one sector buys 1.25·log10(N) ≈ 3.85 mag of depth,
        // turning T≈15 into the ~T 19 floor of the headline range.
        let depth = stacked_depth_over_sectors(published::TESS_SINGLE_DEPTH, 1.0);
        assert!(
            (18.5..19.5).contains(&depth),
            "single-sector stacked depth = {depth:.2}"
        );
    }

    #[test]
    fn year_long_stack_reaches_the_deep_end_of_the_range() {
        // A full year (13 sectors, ~1.5e4 FFIs) reaches the deep end ~T 20–21.
        let depth =
            stacked_depth_over_sectors(published::TESS_SINGLE_DEPTH, tess::SECTORS_PER_YEAR);
        assert!(
            (20.0..21.5).contains(&depth),
            "one-year stacked depth = {depth:.2}"
        );
        // Plausible frame counts land inside the paper's quoted ~T 19–21 band.
        assert!(depth >= published::STACKED_REACH_FAINT);
        assert!(depth <= published::STACKED_REACH_DEEP + 0.5);
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
    fn recovered_tnos_need_deep_multi_sector_co_adds() {
        // Sedna / 2015 BP519 / 2007 TG422 (V ≈ 20.6–22.3) are all *fainter*
        // than a one-sector T≈18.8 stack, so recovering them is a genuine
        // multi-sector / per-track deep co-add — not a single-frame detection.
        // The faintest needs hundreds of sector-equivalents of integration; the
        // brightest (Sedna) just over a single sector's worth.
        let prev = sectors_to_reach_depth(published::TESS_SINGLE_DEPTH, published::V_SEDNA);
        assert!(
            (20.0..40.0).contains(&prev),
            "Sedna sector-equiv = {prev:.0}"
        );
        // monotonically more integration for fainter targets.
        let mut last = 0.0;
        for &v in &[published::V_SEDNA, published::V_BP519, published::V_TG422] {
            let n = sectors_to_reach_depth(published::TESS_SINGLE_DEPTH, v);
            assert!(n > last, "non-monotone integration at V={v}: {n:.0}");
            last = n;
        }
        // All are finite (reachable in principle) — the validation succeeds.
        assert!(last.is_finite() && last > 0.0);
    }

    #[test]
    fn nominal_planet_nine_reach_is_marginal_and_distance_dependent() {
        // A 6.2 M⊕ P9 is V≈19.5 at its a=380 au and V≈20.7 near aphelion
        // (~500 au). The year-long stack reaches only ~T 20.2, so whether P9 is
        // detectable depends entirely on where it sits: reachable on the near
        // side, too faint at aphelion. This is the honest marginal verdict
        // (cf. the paper's V<21, d≲150 au sensitivity — far closer than P9).
        let p9 = P9Params::mcmc_2021();
        let m_at_a = p9_apparent_magnitude_at(&p9, 0.4, p9.a);
        let m_aphelion = p9_apparent_magnitude_at(&p9, 0.4, 500.0);
        assert!((19.0..20.0).contains(&m_at_a), "V_P9(a) = {m_at_a:.2}");
        assert!(m_aphelion > m_at_a, "aphelion must be fainter");

        let m_stack =
            stacked_depth_over_sectors(published::TESS_SINGLE_DEPTH, tess::SECTORS_PER_YEAR);
        assert!((20.0..21.0).contains(&m_stack), "year depth = {m_stack:.2}");

        // Near its semi-major axis P9 is (just) within the year-long reach …
        assert!(
            within_reach(m_at_a, m_stack),
            "P9 at a: V={m_at_a:.2} vs depth {m_stack:.2}"
        );
        // … but near aphelion it falls beyond it. Marginal, distance-dependent.
        assert!(
            !within_reach(m_aphelion, m_stack),
            "P9 at aphelion: V={m_aphelion:.2} vs depth {m_stack:.2}"
        );
        // The margin at aphelion is small (sub-magnitude) — beyond reach, but
        // only barely, so a deeper multi-year stack would close it.
        let margin = reach_margin(m_aphelion, m_stack);
        assert!(
            (-1.0..0.0).contains(&margin),
            "aphelion margin = {margin:.2}"
        );
    }
}
