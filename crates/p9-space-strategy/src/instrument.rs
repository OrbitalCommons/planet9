//! The space telescope being planned for: depth as a function of integration,
//! detection and linking of a slow mover, and the time a tile costs.
//!
//! The default is the JBT 0.5 m + SPENCER concept already catalogued in
//! `p9_survey::telescope` (0.485 m, f/12.3, four IMX455 sensors, 0.13″ pixels,
//! 0.32 deg² field; V ≈ 23.8 at 5σ in one 300 s exposure).
//!
//! **Depth.** Signal-to-noise for a point source in `n` sub-exposures totalling
//! `T` seconds is ∝ `T / sqrt((s + d)·T + n·RN²)`, with `s` the sky rate per
//! pixel, `d` the dark rate and `RN` the read noise. With 0.13″ pixels the sky
//! is so dark per pixel (≈0.014 e⁻/s) that an exposure is read-noise limited
//! below ~3 minutes: short exposures are expensive, and depth grows faster
//! than √T until the sky takes over. Depth is anchored at the catalogued
//! 300 s value.
//!
//! **Sky.** The background is zodiacal light, brightest on the ecliptic. The
//! anti-solar surface brightness is modelled as
//! μ_V(β) = 22.1 + 1.2·(1 − exp(−|β|/25°)) mag/arcsec² (22.1 on the ecliptic,
//! 23.3 at the ecliptic poles; Leinert et al. 1998, Table 16–17 range).
//!
//! **Linking.** A detection is a linked track: the object must be found in at
//! least `link_epochs` of `n_epochs` visits, each reaching the full depth,
//! spaced so it moves at least `min_motion_arcsec` between them.

use p9_core::analysis::surveys::{logistic_efficiency, poisson_binomial_tail};
use p9_survey::plan::JBT_REF_DEPTH_300S;
use p9_survey::telescope::{JBT_PLATE_SCALE_ARCSEC_PER_PX, SPENCER_FIELD_DEG2};
use serde::Serialize;

/// Reference sky surface brightness behind the catalogued depth (V mag/arcsec²).
pub const REFERENCE_SKY_MU: f64 = 22.5;

/// Photo-electrons per second from a V = 0 source through the 0.485 m
/// aperture (≈ 8.8×10⁵ photons s⁻¹ cm⁻² across the visible band × 1847 cm² ×
/// 0.5 end-to-end throughput).
pub const ZERO_POINT_E_PER_S: f64 = 8.1e8;

/// Anti-solar zodiacal surface brightness at ecliptic latitude `beta_deg`.
pub fn zodiacal_mu(beta_deg: f64) -> f64 {
    22.1 + 1.2 * (1.0 - (-beta_deg.abs() / 25.0).exp())
}

/// Apparent rate of a body at heliocentric distance `dist_au`, seen from
/// Earth `phase_deg` away from opposition (arcsec per hour). The reflex of
/// Earth's 29.78 km/s dominates at these distances.
pub fn sky_rate_arcsec_per_hr(dist_au: f64, phase_deg: f64) -> f64 {
    147.8 / dist_au * phase_deg.to_radians().cos().abs()
}

#[derive(Debug, Clone, Serialize)]
pub struct SpaceTelescope {
    pub name: String,
    /// Instantaneous field (deg²).
    pub fov_deg2: f64,
    /// 5σ point-source depth of one reference exposure at the reference sky.
    pub reference_depth: f64,
    pub reference_exposure_s: f64,
    /// Pixel scale (arcsec).
    pub pixel_arcsec: f64,
    /// Dark current (e⁻ s⁻¹ px⁻¹) and read noise (e⁻).
    pub dark_rate: f64,
    pub read_noise_e: f64,
    /// Longest single sub-exposure (cosmic rays, pointing stability).
    pub max_subframe_s: f64,
    /// Slew + settle per visit, and readout per sub-exposure (s).
    pub slew_settle_s: f64,
    pub readout_s: f64,
    /// Fraction of wall-clock time spent integrating or slewing to targets
    /// (the rest is Earth occultation, SAA, downlink, momentum management).
    pub duty_cycle: f64,
    /// Fraction of a tile's area that lands on live silicon once per visit
    /// (chip gaps, dither overlap).
    pub fill_factor: f64,
    /// Radius lost around each star (arcsec).
    pub mask_radius_arcsec: f64,
    /// Logistic steepness of the detection roll-off (mag⁻¹).
    pub efficiency_steepness: f64,
    /// Visits per field and how many must detect the object to link it.
    pub n_epochs: u32,
    pub link_epochs: u32,
    /// Displacement needed between linked visits (arcsec).
    pub min_motion_arcsec: f64,
    /// Fields within this angle of the Sun cannot be observed (deg).
    pub sun_exclusion_deg: f64,
}

impl Default for SpaceTelescope {
    fn default() -> Self {
        Self {
            name: "JBT 0.5 m + SPENCER".to_string(),
            fov_deg2: SPENCER_FIELD_DEG2,
            reference_depth: JBT_REF_DEPTH_300S,
            reference_exposure_s: 300.0,
            pixel_arcsec: JBT_PLATE_SCALE_ARCSEC_PER_PX,
            dark_rate: 0.0046,
            read_noise_e: 1.58,
            max_subframe_s: 600.0,
            slew_settle_s: 45.0,
            readout_s: 5.0,
            duty_cycle: 0.6,
            fill_factor: 0.9,
            mask_radius_arcsec: 0.6,
            efficiency_steepness: 4.0,
            n_epochs: 4,
            link_epochs: 3,
            min_motion_arcsec: 1.0,
            sun_exclusion_deg: 60.0,
        }
    }
}

impl SpaceTelescope {
    /// Sky electrons per second per pixel at surface brightness `mu`.
    pub fn sky_rate(&self, mu: f64) -> f64 {
        ZERO_POINT_E_PER_S * 10f64.powf(-0.4 * mu) * self.pixel_arcsec.powi(2)
    }

    fn subframes(&self, integration_s: f64) -> f64 {
        (integration_s / self.max_subframe_s).ceil().max(1.0)
    }

    /// Relative signal-to-noise of a fixed source: T / sqrt(noise variance).
    fn snr_shape(&self, integration_s: f64, mu: f64) -> f64 {
        let variance = (self.sky_rate(mu) + self.dark_rate) * integration_s
            + self.subframes(integration_s) * self.read_noise_e.powi(2);
        integration_s / variance.sqrt()
    }

    /// 5σ depth (V) of one visit of `integration_s` seconds at ecliptic
    /// latitude `beta_deg`.
    pub fn depth(&self, integration_s: f64, beta_deg: f64) -> f64 {
        let reference = self.snr_shape(self.reference_exposure_s, REFERENCE_SKY_MU);
        let here = self.snr_shape(integration_s, zodiacal_mu(beta_deg));
        self.reference_depth + 2.5 * (here / reference).log10()
    }

    /// Probability a source of magnitude `v` is linked in a field whose
    /// per-visit depth is `depth`.
    pub fn link_probability(&self, v: f64, depth: f64) -> f64 {
        let eps = logistic_efficiency(v, depth, self.efficiency_steepness);
        poisson_binomial_tail(&vec![eps; self.n_epochs as usize], self.link_epochs)
    }

    /// Telescope seconds one field costs at `integration_s` per visit.
    pub fn field_seconds(&self, integration_s: f64) -> f64 {
        self.n_epochs as f64
            * (integration_s + self.slew_settle_s + self.subframes(integration_s) * self.readout_s)
    }

    /// Fields needed to cover `area_deg2`.
    pub fn fields_for(&self, area_deg2: f64) -> f64 {
        area_deg2 / (self.fov_deg2 * self.fill_factor)
    }

    /// Wall-clock hours to survey `area_deg2` at `integration_s` per visit.
    pub fn wall_hours(&self, area_deg2: f64, integration_s: f64) -> f64 {
        self.fields_for(area_deg2) * self.field_seconds(integration_s) / self.duty_cycle / 3600.0
    }

    /// Shortest gap between linked visits for a body at `dist_au`, observed
    /// `phase_deg` from opposition (hours).
    pub fn min_baseline_hr(&self, dist_au: f64, phase_deg: f64) -> f64 {
        self.min_motion_arcsec / sky_rate_arcsec_per_hr(dist_au, phase_deg)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn depth_is_anchored_and_monotone() {
        let t = SpaceTelescope::default();
        // Reference exposure at the reference sky (μ = 22.5 ↔ |β| ≈ 10°).
        let beta_ref = -25.0 * (1.0 - (REFERENCE_SKY_MU - 22.1) / 1.2_f64).ln();
        assert!((t.depth(300.0, beta_ref) - t.reference_depth).abs() < 1e-9);
        let mut last = 0.0;
        for s in [30.0, 60.0, 120.0, 300.0, 600.0, 1200.0, 2400.0] {
            let d = t.depth(s, 10.0);
            assert!(d > last);
            last = d;
        }
    }

    #[test]
    fn short_exposures_are_read_noise_limited() {
        // Halving a 60 s exposure costs nearly 2.5·log10(2) = 0.75 mag (read
        // noise), while halving a 2400 s one costs nearer 0.38 mag (sky).
        let t = SpaceTelescope::default();
        let short = t.depth(60.0, 10.0) - t.depth(30.0, 10.0);
        let long = t.depth(2400.0, 10.0) - t.depth(1200.0, 10.0);
        assert!(short > 0.6, "short-exposure gain {short}");
        assert!(long < 0.45, "long-exposure gain {long}");
    }

    #[test]
    fn the_ecliptic_is_shallower_than_the_pole() {
        let t = SpaceTelescope::default();
        assert!(t.depth(600.0, 0.0) < t.depth(600.0, 60.0));
        assert!(t.depth(600.0, 60.0) - t.depth(600.0, 0.0) < 0.7);
    }

    #[test]
    fn planet_nine_needs_hours_between_visits() {
        let t = SpaceTelescope::default();
        // 0.30"/hr at 500 AU at opposition; a TNO at 40 AU moves 3.7"/hr.
        assert!((sky_rate_arcsec_per_hr(500.0, 0.0) - 0.2956).abs() < 1e-3);
        assert!((t.min_baseline_hr(500.0, 0.0) - 3.38).abs() < 0.05);
        assert!(t.min_baseline_hr(1000.0, 45.0) > 9.0);
    }

    #[test]
    fn linking_needs_three_of_four() {
        let t = SpaceTelescope::default();
        assert!(t.link_probability(22.0, 24.0) > 0.999);
        assert!(t.link_probability(26.0, 24.0) < 1e-6);
        // At the single-visit 50% depth, 3 of 4 succeeds 5/16 of the time.
        assert!((t.link_probability(24.0, 24.0) - 0.3125).abs() < 1e-9);
    }
}
