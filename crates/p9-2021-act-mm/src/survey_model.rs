//! ACT millimeter survey parameters: frequencies, point-source sensitivity,
//! footprint, and the published candidate count.
//!
//! Naess et al. (2021), "The Atacama Cosmology Telescope: A search for Planet
//! 9" (arXiv:2104.10264), searched ACT maps at 98, 150 and 229 GHz over
//! ~18,000 deg² for a slowly moving mm point source. The 150 GHz channel sets
//! the headline 95%-confidence point-source flux limit of 4–12 mJy (position
//! dependent). The blind search returned a list of candidate sources of which
//! the 10 strongest were flagged for follow-up; none were confirmed as Planet
//! Nine. All numbers below are labelled reference constants from the paper.

/// ACT 90 GHz ("f090") band center, Hz.
pub const NU_98_GHZ: f64 = 98.0e9;
/// ACT 150 GHz ("f150") band center, Hz — the headline search band.
pub const NU_150_GHZ: f64 = 150.0e9;
/// ACT 220 GHz ("f220", 229 GHz effective) band center, Hz.
pub const NU_229_GHZ: f64 = 229.0e9;

/// 95%-confidence point-source flux limit at 150 GHz (mJy). Position dependent;
/// the paper quotes 4–12 mJy. We adopt a representative value for the nominal
/// detectability calculation and pin the published range in tests.
pub const FLUX_LIMIT_150_MJY: f64 = 8.0;

/// Published lower/upper bounds of the position-dependent 150 GHz flux limit
/// (mJy): "4–12 mJy depending on position" (Naess et al. 2021).
pub const FLUX_LIMIT_150_MIN_MJY: f64 = 4.0;
pub const FLUX_LIMIT_150_MAX_MJY: f64 = 12.0;

/// Search footprint in square degrees (Naess et al. 2021, ~18,000 deg²).
pub const SEARCH_AREA_DEG2: f64 = 18_000.0;

/// Number of strongest candidate sources flagged for follow-up (Naess et al.
/// 2021). None were confirmed as Planet Nine.
pub const CANDIDATE_COUNT: usize = 10;

/// The ACT mm point-source search configuration.
#[derive(Debug, Clone, Copy)]
pub struct ActSurvey {
    /// 95%-confidence point-source flux limit at the search frequency (mJy).
    pub flux_limit_mjy: f64,
    /// Search frequency (Hz).
    pub frequency_hz: f64,
    /// Footprint area (deg²).
    pub area_deg2: f64,
}

impl Default for ActSurvey {
    /// The nominal 150 GHz ACT search at the representative flux limit over the
    /// full ~18,000 deg² footprint.
    fn default() -> Self {
        Self {
            flux_limit_mjy: FLUX_LIMIT_150_MJY,
            frequency_hz: NU_150_GHZ,
            area_deg2: SEARCH_AREA_DEG2,
        }
    }
}

impl ActSurvey {
    /// The ACT search at the deepest (best, 4 mJy) position-dependent limit.
    pub fn deepest() -> Self {
        Self {
            flux_limit_mjy: FLUX_LIMIT_150_MIN_MJY,
            ..Self::default()
        }
    }

    /// The ACT search at the shallowest (12 mJy) position-dependent limit.
    pub fn shallowest() -> Self {
        Self {
            flux_limit_mjy: FLUX_LIMIT_150_MAX_MJY,
            ..Self::default()
        }
    }

    /// Total sky area in steradians (for spurious-rate estimates).
    pub fn area_sr(&self) -> f64 {
        let deg2_per_sr = (180.0 / std::f64::consts::PI).powi(2);
        self.area_deg2 / deg2_per_sr
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn representative_limit_lies_in_published_range() {
        assert!((FLUX_LIMIT_150_MIN_MJY..=FLUX_LIMIT_150_MAX_MJY).contains(&FLUX_LIMIT_150_MJY));
    }

    #[test]
    fn frequencies_ordered_and_in_mm_band() {
        assert!(NU_98_GHZ < NU_150_GHZ && NU_150_GHZ < NU_229_GHZ);
        // mm/sub-mm band: tens to hundreds of GHz.
        assert!(NU_98_GHZ > 50.0e9 && NU_229_GHZ < 400.0e9);
    }

    #[test]
    fn footprint_is_thousands_of_square_degrees() {
        // ~18,000 deg² out of 41,253 deg² all-sky: ~44% of the sky.
        let all_sky = 41_252.96;
        let frac = SEARCH_AREA_DEG2 / all_sky;
        assert!((0.3..0.5).contains(&frac), "footprint fraction = {frac:.2}");
    }

    #[test]
    fn area_in_steradians_consistent() {
        let s = ActSurvey::default();
        // 18,000 deg² ≈ 5.48 sr.
        assert!(
            (s.area_sr() - 5.48).abs() < 0.05,
            "area = {:.3} sr",
            s.area_sr()
        );
    }

    #[test]
    fn deepest_and_shallowest_bracket_the_default() {
        assert!(ActSurvey::deepest().flux_limit_mjy < ActSurvey::default().flux_limit_mjy);
        assert!(ActSurvey::shallowest().flux_limit_mjy > ActSurvey::default().flux_limit_mjy);
    }
}
