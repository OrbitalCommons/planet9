//! DES survey characteristics for the Planet Nine search of
//! Belyakov, Bernardinelli & Brown (2022), AJ 163, 216 (arXiv:2203.07642).
//!
//! Models the Dark Energy Survey wide-field footprint (as a union of
//! declination bands whose solid angle matches the surveyed ~5000 deg^2),
//! per-band difference-imaging depths, and the paper's multi-night
//! detection/linking requirement ("more than 7 unique nights of
//! detections" across the 6-year baseline).

use serde::{Deserialize, Serialize};

use p9_core::analysis::surveys::{
    des_footprint_contains, des_footprint_solid_angle_deg2, limiting_magnitude,
    poisson_binomial_tail,
};
use p9_core::coords::observer::{EarthProvider, EarthState, Time, Timescale};
use p9_core::coords::sky::{
    apparent_position_deg, apparent_position_with_earth_deg, phase_angle_with_earth,
};
use p9_core::types::OrbitalElements;

use crate::color_models::BandMagnitudes;
use crate::sky::earth_states_at;

/// Photometric band used in DES observations.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum DesBand {
    G,
    R,
    I,
    Z,
    Y,
}

/// DES survey model parameters for Planet Nine detection.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DesSurvey {
    /// Declared effective footprint area in square degrees (the box-union
    /// footprint's computed solid angle matches this to <1%; see
    /// `solid_angle_deg2`).
    pub footprint_area: f64,
    /// Total number of distinct observing nights across the 6-year survey;
    /// per-position observing epochs are quantized onto this night grid.
    pub n_nights: u32,
    /// Distinct nights on which any given sky position is imaged *per band*
    /// (~10 tilings of the wide survey; Diehl et al. 2016, Neilsen et al.
    /// 2019).
    pub nights_per_band: u32,
    /// Detections on at least this many distinct nights are required to
    /// form a linkable orbit: the paper's criterion is NUNIQUE >= 7
    /// (inherited from Bernardinelli et al. 2022, "appear in at least 7
    /// nights of data"). Belyakov et al.'s prose "more than 7 unique nights"
    /// contradicts its own math; a previous version resolved it as >= 8 — a
    /// strictly harsher cut whose deficit the night_usability calibration
    /// silently absorbed.
    pub min_nights: u32,
    /// 50%-completeness magnitude in g: 24.1, the logit-fit m50 measured by
    /// injecting synthetic Planet Nines into the DES difference-imaging
    /// pipeline (Belyakov et al. 2022, Sec. 3).
    pub depth_g: f64,
    /// 50%-completeness magnitude in r: 23.8 (Bernardinelli et al. 2022) —
    /// the value shared workspace-wide via
    /// `p9_core::analysis::surveys::SURVEY_DEPTHS`.
    pub depth_r: f64,
    /// 50%-completeness magnitude in i: DES DR1 single-epoch depths run
    /// g/r/i/z = 23.57/23.34/22.78/22.10 (Abbott et al. 2018); scaling the
    /// i-g offset (-0.79) to the difference-imaging g anchor gives ~23.3.
    pub depth_i: f64,
    /// 50%-completeness magnitude in z: DR1 z-g offset (-1.47) scaled to
    /// the difference-imaging g anchor gives ~22.6.
    pub depth_z: f64,
    /// Calibrated per-night usability probability: the chance that a given
    /// tiling night yields a clean, linkable detection epoch for a bright
    /// source. It lumps focal-plane fill (~20% chip gaps/masks, Morganson
    /// et al. 2018), template contamination in difference imaging, and the
    /// arc-length/triplet structure cuts the pipeline applies but this
    /// model does not resolve. Calibrated (single parameter) so that the
    /// end-to-end injection-recovery of the Brown & Batygin (2021)
    /// reference population reproduces the paper's 87.0% recovery; see
    /// `recovery_analysis` for the regression test.
    pub night_usability: f64,
    /// Bands used in the detection simulation
    pub bands: Vec<DesBand>,
}

impl Default for DesSurvey {
    fn default() -> Self {
        Self {
            footprint_area: 5000.0,
            n_nights: 575,
            nights_per_band: 10,
            min_nights: 7,
            depth_g: 24.1,
            depth_r: limiting_magnitude("DES").expect("DES depth in p9-core survey table"),
            depth_i: 23.3,
            depth_z: 22.6,
            night_usability: 0.255,
            bands: vec![DesBand::G, DesBand::R, DesBand::I, DesBand::Z],
        }
    }
}

impl DesSurvey {
    /// 50%-completeness magnitude for a band.
    pub fn depth(&self, band: DesBand) -> f64 {
        match band {
            DesBand::G => self.depth_g,
            DesBand::R => self.depth_r,
            DesBand::I => self.depth_i,
            DesBand::Z => self.depth_z,
            // Y is shallow and unused by the detection model.
            DesBand::Y => 21.0,
        }
    }

    /// Single-epoch detection completeness at apparent magnitude `m` in
    /// `band`:
    ///
    ///   p(m) = c / (1 + exp(k (m - m50)))
    ///
    /// the logit form fitted by Belyakov et al. (2022, Sec. 3) to synthetic
    /// P9 injections, with asymptotic completeness c = 0.95, slope k = 3.0
    /// (transition width ~0.7 mag) and m50 the band's 50%-completeness
    /// depth.
    pub fn completeness(&self, apparent_mag: f64, band: DesBand) -> f64 {
        let c = 0.95;
        let k = 3.0;
        let m50 = self.depth(band);
        c / (1.0 + (k * (apparent_mag - m50)).exp())
    }

    /// Probability that a single tiling night yields a usable detection:
    /// the magnitude completeness degraded by the calibrated per-night
    /// usability factor.
    pub fn night_detection_probability(&self, apparent_mag: f64, band: DesBand) -> f64 {
        self.night_usability * self.completeness(apparent_mag, band)
    }

    /// Footprint membership (equatorial RA/dec in degrees).
    pub fn is_in_footprint(&self, ra_deg: f64, dec_deg: f64) -> bool {
        des_footprint_contains(ra_deg, dec_deg)
    }

    /// Exact solid angle of the box-union footprint in deg^2:
    /// sum over bands of dRA * (sin(dec_max) - sin(dec_min)) * 180/pi.
    pub fn solid_angle_deg2(&self) -> f64 {
        des_footprint_solid_angle_deg2()
    }

    /// Observing epochs for one sky position: `nights_per_band` nights in
    /// each band, staggered between bands and quantized onto the survey's
    /// `n_nights`-night grid spanning the 6-year baseline. Returns
    /// (band, time in days since survey start).
    pub fn observing_epochs(&self) -> Vec<(DesBand, f64)> {
        let span_days = 6.0 * p9_core::constants::YEAR_DAYS;
        let nights = self.nights_per_band.max(1) as f64;
        let grid = (self.n_nights.max(2) - 1) as f64;
        let mut epochs = Vec::with_capacity(self.bands.len() * self.nights_per_band as usize);
        for (b_idx, &band) in self.bands.iter().enumerate() {
            let offset = b_idx as f64 / self.bands.len() as f64;
            for j in 0..self.nights_per_band {
                let frac = (j as f64 + offset) / nights;
                // Snap to the discrete observing-night grid.
                let night_index = (frac * grid).round();
                epochs.push((band, night_index / grid * span_days));
            }
        }
        epochs
    }

    /// Survey-start anchor of the night grid (`sky::survey_start`,
    /// 2013-08-31): `observing_epochs` day offsets count from this epoch.
    pub fn survey_start(&self, ts: &Timescale) -> Time {
        crate::sky::survey_start(ts)
    }

    /// Observing epochs with precomputed per-night Earth states from a
    /// real ephemeris (or the analytic fallback provider). The night grid
    /// is orbit-independent, so compute this once per recovery run and
    /// share it across orbits.
    pub fn observing_epochs_with_earth(
        &self,
        earth: &mut impl EarthProvider,
    ) -> Vec<(DesBand, f64, EarthState)> {
        let ts = Timescale::default();
        let start = self.survey_start(&ts);
        let epochs = self.observing_epochs();
        let offsets: Vec<f64> = epochs.iter().map(|&(_, t)| t).collect();
        let states = earth_states_at(earth, &ts, &start, &offsets);
        epochs
            .into_iter()
            .zip(states)
            .map(|((band, t), state)| (band, t, state))
            .collect()
    }

    /// Per-night detection probabilities for an orbit across all observing
    /// epochs (D3): each tiling night contributes its magnitude
    /// completeness in that night's band, gated by whether the object's
    /// *apparent* position that night — including parallactic oscillation
    /// and orbital drift — falls inside the footprint.
    pub fn night_probabilities(&self, elem: &OrbitalElements, mags: &BandMagnitudes) -> Vec<f64> {
        self.observing_epochs()
            .iter()
            .map(|&(band, t_days)| {
                let (ra, dec, _) = apparent_position_deg(elem, t_days);
                if self.is_in_footprint(ra, dec) {
                    self.night_detection_probability(mags.magnitude(band), band)
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Ephemeris-path counterpart of `night_probabilities`: per-night
    /// apparent positions from real Earth states (parallax from the true
    /// orbit, light-time, aberration), and per-night IAU H-G phase
    /// darkening (G = 0.15) relative to the zero-phase `mags` — the H-G
    /// opposition surge contributes up to ~0.05 mag at P9 distances
    /// (α ≤ asin(1 AU/r)). `nights` comes from
    /// `observing_epochs_with_earth`.
    pub fn night_probabilities_with_earth(
        &self,
        elem: &OrbitalElements,
        mags: &BandMagnitudes,
        nights: &[(DesBand, f64, EarthState)],
    ) -> Vec<f64> {
        use p9_core::analysis::photometry::{hg_phase_factor, DEFAULT_SLOPE_G};
        nights
            .iter()
            .map(|(band, t_days, state)| {
                let (ra, dec, _) = apparent_position_with_earth_deg(elem, state, *t_days);
                if self.is_in_footprint(ra, dec) {
                    let alpha = phase_angle_with_earth(elem, state, *t_days);
                    let dm = -2.5 * hg_phase_factor(DEFAULT_SLOPE_G, alpha).log10();
                    self.night_detection_probability(mags.magnitude(*band) + dm, *band)
                } else {
                    0.0
                }
            })
            .collect()
    }

    /// Whether the orbit's apparent position enters the footprint on any
    /// observing epoch (the paper's "crosses the DES footprint").
    pub fn crosses_footprint(&self, elem: &OrbitalElements) -> bool {
        self.observing_epochs().iter().any(|&(_, t_days)| {
            let (ra, dec, _) = apparent_position_deg(elem, t_days);
            self.is_in_footprint(ra, dec)
        })
    }

    /// Ephemeris-path counterpart of `crosses_footprint`.
    pub fn crosses_footprint_with_earth(
        &self,
        elem: &OrbitalElements,
        nights: &[(DesBand, f64, EarthState)],
    ) -> bool {
        nights.iter().any(|(_, t_days, state)| {
            let (ra, dec, _) = apparent_position_with_earth_deg(elem, state, *t_days);
            self.is_in_footprint(ra, dec)
        })
    }

    /// End-to-end probability that the orbit is detected and linked:
    /// P(detections on >= `min_nights` distinct nights), with each night an
    /// independent Bernoulli trial at `night_probabilities` (exact
    /// Poisson-binomial tail).
    pub fn detection_probability_for_orbit(
        &self,
        elem: &OrbitalElements,
        mags: &BandMagnitudes,
    ) -> f64 {
        let probs = self.night_probabilities(elem, mags);
        poisson_binomial_tail(&probs, self.min_nights)
    }

    /// Ephemeris-path counterpart of `detection_probability_for_orbit`.
    pub fn detection_probability_for_orbit_with_earth(
        &self,
        elem: &OrbitalElements,
        mags: &BandMagnitudes,
        nights: &[(DesBand, f64, EarthState)],
    ) -> f64 {
        let probs = self.night_probabilities_with_earth(elem, mags, nights);
        poisson_binomial_tail(&probs, self.min_nights)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_matches_paper_values() {
        let survey = DesSurvey::default();
        assert!((survey.footprint_area - 5000.0).abs() < 1e-10);
        assert_eq!(survey.n_nights, 575);
        // g-band m50 = 24.1 (Belyakov et al. 2022 logit fit).
        assert!((survey.depth_g - 24.1).abs() < 1e-10);
        // r depth comes from the shared p9-core table (23.8).
        assert!((survey.depth_r - 23.8).abs() < 1e-10);
        // NUNIQUE >= 7 (the paper's stated criterion; the prose "more than
        // 7" contradicts its own math and a previous version pinned 8).
        assert_eq!(survey.min_nights, 7);
    }

    #[test]
    fn footprint_solid_angle_matches_declared_area() {
        // D2 regression: the footprint cut must integrate to the declared
        // 5000 deg^2 (the old box covered ~9300 deg^2).
        let survey = DesSurvey::default();
        let omega = survey.solid_angle_deg2();
        assert!(
            (omega - survey.footprint_area).abs() / survey.footprint_area < 0.02,
            "solid angle = {omega:.0} deg^2 vs declared {}",
            survey.footprint_area
        );
    }

    #[test]
    fn footprint_solid_angle_monte_carlo_agrees() {
        // Cross-check the analytic solid angle with uniform-sphere sampling.
        use rand::{Rng, SeedableRng};
        let survey = DesSurvey::default();
        let mut rng = rand::rngs::StdRng::seed_from_u64(2022);
        let n = 200_000;
        let mut hits = 0u32;
        for _ in 0..n {
            let ra = rng.gen_range(0.0..360.0);
            let z: f64 = rng.gen_range(-1.0..1.0);
            let dec = z.asin().to_degrees();
            if survey.is_in_footprint(ra, dec) {
                hits += 1;
            }
        }
        let full_sky = 4.0 * std::f64::consts::PI * (180.0 / std::f64::consts::PI).powi(2);
        let omega_mc = full_sky * hits as f64 / n as f64;
        assert!(
            (omega_mc - survey.solid_angle_deg2()).abs() < 150.0,
            "MC = {omega_mc:.0} vs analytic {:.0}",
            survey.solid_angle_deg2()
        );
    }

    #[test]
    fn bright_objects_high_completeness() {
        let survey = DesSurvey::default();
        let p = survey.completeness(20.0, DesBand::R);
        assert!(
            p > 0.90,
            "Bright objects should have >90% completeness, got {p}"
        );
    }

    #[test]
    fn faint_objects_low_completeness() {
        let survey = DesSurvey::default();
        let p = survey.completeness(26.0, DesBand::R);
        assert!(
            p < 0.10,
            "Faint objects should have <10% completeness, got {p}"
        );
    }

    #[test]
    fn completeness_monotonically_decreasing() {
        let survey = DesSurvey::default();
        let mut prev = survey.completeness(18.0, DesBand::G);
        for mag_tenth in 185..=270 {
            let mag = mag_tenth as f64 / 10.0;
            let p = survey.completeness(mag, DesBand::G);
            assert!(
                p <= prev + 1e-12,
                "Completeness should decrease with magnitude"
            );
            prev = p;
        }
    }

    #[test]
    fn footprint_rejects_northern_sky() {
        let survey = DesSurvey::default();
        assert!(!survey.is_in_footprint(30.0, 10.0));
        assert!(!survey.is_in_footprint(30.0, 45.0));
    }

    #[test]
    fn footprint_accepts_spt_region() {
        let survey = DesSurvey::default();
        assert!(survey.is_in_footprint(30.0, -55.0));
        assert!(survey.is_in_footprint(350.0, -50.0));
        assert!(survey.is_in_footprint(100.0, -45.0));
    }

    #[test]
    fn footprint_rejects_mid_north_away_from_stripe82() {
        let survey = DesSurvey::default();
        // dec -30 at RA 150 is outside the wide survey.
        assert!(!survey.is_in_footprint(150.0, -30.0));
        // dec 0 belongs to the footprint only near RA ~ 0.
        assert!(survey.is_in_footprint(10.0, 0.0));
        assert!(!survey.is_in_footprint(120.0, 0.0));
    }

    #[test]
    fn footprint_rejects_deep_south() {
        let survey = DesSurvey::default();
        assert!(!survey.is_in_footprint(30.0, -70.0));
    }

    #[test]
    fn band_depths_ordered() {
        // g (difference-imaging anchor) is deepest, then r, i, z.
        let survey = DesSurvey::default();
        assert!(survey.depth(DesBand::G) > survey.depth(DesBand::R));
        assert!(survey.depth(DesBand::R) > survey.depth(DesBand::I));
        assert!(survey.depth(DesBand::I) > survey.depth(DesBand::Z));
    }

    #[test]
    fn observing_epochs_span_survey_and_use_night_grid() {
        let survey = DesSurvey::default();
        let epochs = survey.observing_epochs();
        assert_eq!(
            epochs.len(),
            survey.bands.len() * survey.nights_per_band as usize
        );
        let span = 6.0 * p9_core::constants::YEAR_DAYS;
        for &(_, t) in &epochs {
            assert!((0.0..=span).contains(&t));
            // Quantized onto the n_nights grid.
            let grid = (survey.n_nights - 1) as f64;
            let idx = t / span * grid;
            assert!((idx - idx.round()).abs() < 1e-9, "t = {t} off-grid");
        }
    }
}
