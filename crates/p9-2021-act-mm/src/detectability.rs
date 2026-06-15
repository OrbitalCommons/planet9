//! Maximum mm detectable distance, the mm-vs-WISE reach comparison, and an
//! estimate of the number of spurious candidates over the ACT footprint.
//!
//! Because the mm flux density falls monotonically as 1/d² (radius and the
//! cold internal temperature being essentially distance-independent beyond a
//! few hundred AU), the maximum detectable distance has a closed form: the
//! flux equals the survey limit when
//!
//!   F_ν(d) = π B_ν(T) (R/d)² = S_lim   ⟹   d_max = R √(π B_ν(T) / S_lim).
//!
//! We solve it directly (and cross-check by bisection in the tests). The
//! headline of Naess et al. (2021) is that this d_max is a *large* number —
//! hundreds of AU — whereas the same cold body is essentially invisible in
//! WISE W1 (where the thermal channel is buried under the Wien cliff and only
//! reflected sunlight survives).

use p9_2018_wise_search::detectability::max_detectable_distance as wise_max_distance;
use p9_2018_wise_search::survey_model::WiseSurvey;
use p9_core::analysis::thermal::planck_bnu;
use p9_core::constants::AU_M;

use crate::survey_model::ActSurvey;
use crate::thermal_model::{P9Millimeter, MJY};

/// Whether a P9 of `mass_earth` at `distance_au` is at or above the ACT flux
/// limit at the survey frequency.
pub fn is_detectable(survey: &ActSurvey, mass_earth: f64, distance_au: f64) -> bool {
    P9Millimeter::new(mass_earth, distance_au).flux_density_mjy(survey.frequency_hz)
        >= survey.flux_limit_mjy
}

/// Maximum heliocentric distance (AU) at which a P9 of `mass_earth` reaches the
/// ACT flux limit, from the closed-form 1/d² inversion. The temperature/radius
/// are taken at a fiducial 500 AU (cold internal floor, distance-independent).
pub fn max_detectable_distance(survey: &ActSurvey, mass_earth: f64) -> f64 {
    let p = P9Millimeter::new(mass_earth, 500.0);
    let bnu = planck_bnu(p.effective_temp(), survey.frequency_hz);
    let r = p.radius_m();
    let s_lim = survey.flux_limit_mjy * MJY; // W/m²/Hz
    let d_m = r * (std::f64::consts::PI * bnu / s_lim).sqrt();
    d_m / AU_M
}

/// The mm-band reach for a cold P9 at the ACT 150 GHz limit, paired with the
/// WISE W1 reach for the *same* body. The mm number should vastly exceed W1.
#[derive(Debug, Clone, Copy)]
pub struct ReachComparison {
    pub mass_earth: f64,
    /// ACT mm maximum detectable distance (AU).
    pub act_mm_au: f64,
    /// WISE W1 maximum detectable distance (AU), or None if undetectable.
    pub wise_w1_au: Option<f64>,
}

/// Compare ACT mm vs WISE W1 reach for a given mass.
pub fn compare_reach(act: &ActSurvey, mass_earth: f64) -> ReachComparison {
    ReachComparison {
        mass_earth,
        act_mm_au: max_detectable_distance(act, mass_earth),
        wise_w1_au: wise_max_distance(&WiseSurvey::default(), mass_earth),
    }
}

/// Expected number of spurious candidates over the search footprint, given a
/// surface density `sources_per_deg2` of noise peaks / confusing point sources
/// above the detection threshold. A pure-Poisson area scaling:
///
///   N_spurious = ρ · A.
///
/// This is the quantity an mm point-source search must beat down (by the
/// slow-motion linking that distinguishes a real outer-solar-system body from
/// a static extragalactic source) before any candidate can be believed. The
/// 10 published candidates are the strongest survivors of exactly this kind of
/// winnowing.
pub fn expected_spurious(survey: &ActSurvey, sources_per_deg2: f64) -> f64 {
    sources_per_deg2 * survey.area_deg2
}

/// Expected number of pure-Gaussian-noise upward fluctuations above a
/// `n_sigma` detection threshold, given `n_pixels_per_deg2` independent
/// resolution elements per square degree. Uses the one-sided Gaussian tail
/// probability Q(nσ) = ½ erfc(nσ/√2):
///
///   N = Q(nσ) · n_pixels_per_deg2 · A.
pub fn expected_noise_peaks(survey: &ActSurvey, n_sigma: f64, n_pixels_per_deg2: f64) -> f64 {
    let q = 0.5 * erfc(n_sigma / std::f64::consts::SQRT_2);
    q * n_pixels_per_deg2 * survey.area_deg2
}

/// Complementary error function via the Abramowitz & Stegun 7.1.26 rational
/// approximation (|error| < 1.5e-7), sufficient for tail-count estimates.
fn erfc(x: f64) -> f64 {
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.3275911 * z);
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let erf = 1.0 - poly * (-z * z).exp();
    if x >= 0.0 {
        1.0 - erf
    } else {
        1.0 + erf
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::survey_model::{CANDIDATE_COUNT, NU_150_GHZ};

    #[test]
    fn closed_form_matches_bisection() {
        let survey = ActSurvey::default();
        let mass = 10.0;
        let d_closed = max_detectable_distance(&survey, mass);

        // Independent bisection on the monotone flux–distance relation.
        let (mut lo, mut hi) = (50.0_f64, 50_000.0_f64);
        for _ in 0..200 {
            let mid = 0.5 * (lo + hi);
            if is_detectable(&survey, mass, mid) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let d_bisect = 0.5 * (lo + hi);
        assert!(
            (d_closed - d_bisect).abs() / d_closed < 1e-4,
            "closed {d_closed:.2} vs bisection {d_bisect:.2} AU"
        );
    }

    #[test]
    fn crossover_flux_equals_limit() {
        let survey = ActSurvey::default();
        let d = max_detectable_distance(&survey, 10.0);
        let f = P9Millimeter::new(10.0, d).flux_density_mjy(NU_150_GHZ);
        assert!(
            (f - survey.flux_limit_mjy).abs() < 1e-6,
            "F({d:.1} AU) = {f:.4} mJy"
        );
    }

    #[test]
    fn mm_reach_is_hundreds_of_au_and_finite() {
        let survey = ActSurvey::default();
        let d = max_detectable_distance(&survey, 10.0);
        assert!(d.is_finite());
        // A 10 M⊕ cold P9 is detectable out to several hundred AU at 8 mJy,
        // consistent with Naess et al.'s 425–775 AU (4–12 mJy) range for 10 M⊕.
        assert!((300.0..900.0).contains(&d), "d_max(10 M⊕) = {d:.0} AU");
    }

    #[test]
    fn deeper_limit_reaches_farther() {
        let m = 10.0;
        let d_deep = max_detectable_distance(&ActSurvey::deepest(), m);
        let d_shallow = max_detectable_distance(&ActSurvey::shallowest(), m);
        assert!(
            d_deep > d_shallow,
            "deep {d_deep:.0} > shallow {d_shallow:.0} AU"
        );
        // The 4 mJy limit reaches √(12/4) = √3 ≈ 1.73× farther than 12 mJy.
        assert!((d_deep / d_shallow - 3.0_f64.sqrt()).abs() < 1e-6);
    }

    #[test]
    fn reach_grows_with_mass() {
        let survey = ActSurvey::default();
        let d5 = max_detectable_distance(&survey, 5.0);
        let d10 = max_detectable_distance(&survey, 10.0);
        assert!(d10 > d5, "10 M⊕ reach {d10:.0} > 5 M⊕ reach {d5:.0} AU");
    }

    #[test]
    fn mm_vastly_outreaches_wise_w1_for_a_cold_body() {
        // The headline: at the same cold (40 K) temperature, the mm band
        // reaches a P9 at far greater distance than WISE W1, which for a cold
        // body is essentially blind (W1 reach is set by reflected light and is
        // limited to the low hundreds of AU even at high mass).
        let cmp = compare_reach(&ActSurvey::default(), 10.0);
        let w1 = cmp.wise_w1_au.expect("W1 reach finite for 10 M⊕");
        assert!(
            cmp.act_mm_au > 2.0 * w1,
            "ACT mm reach {:.0} AU should far exceed WISE W1 {:.0} AU",
            cmp.act_mm_au,
            w1
        );
    }

    #[test]
    fn spurious_count_scales_with_area_and_density() {
        let survey = ActSurvey::default();
        // A modest 1e-3 confusing sources per deg² over 18,000 deg² → 18.
        let n = expected_spurious(&survey, 1.0e-3);
        assert!((n - 18.0).abs() < 1e-9, "N_spurious = {n}");
        // Linear in density.
        assert!((expected_spurious(&survey, 2.0e-3) - 2.0 * n).abs() < 1e-9);
    }

    #[test]
    fn noise_peak_count_is_plausible_for_a_candidate_list() {
        // A high (5σ) threshold over a finely sampled map still yields a
        // handful of pure-noise peaks across 18,000 deg² — the population from
        // which the 10 published candidates were winnowed. With ~1e4 beams per
        // deg² and 5σ (Q ≈ 2.87e-7), N ≈ 5e2; the slow-motion + multi-frequency
        // linking knocks that down to the ~10 strongest survivors.
        let survey = ActSurvey::default();
        let n = expected_noise_peaks(&survey, 5.0, 1.0e4);
        assert!(n > 1.0, "expected some 5σ noise peaks, got {n:.1}");
        // Sanity: the published candidate list is short.
        assert_eq!(CANDIDATE_COUNT, 10);
    }

    #[test]
    fn erfc_reference_values() {
        // erfc(0) = 1, erfc large → 0, and the 5σ one-sided tail.
        assert!((super::erfc(0.0) - 1.0).abs() < 1e-6);
        assert!(super::erfc(5.0) < 1e-10);
        let q5 = 0.5 * super::erfc(5.0 / std::f64::consts::SQRT_2);
        assert!((q5 - 2.866e-7).abs() < 1e-8, "Q(5σ) = {q5:e}");
    }
}
