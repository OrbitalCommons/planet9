//! Detectability forecast: maximum detectable heliocentric distance, expected
//! detection significance, and the fraction of a reference (mass, distance)
//! box a millimeter survey can detect — for SO vs ACT.
//!
//! All quantities are computed from the mm thermal flux (`crate::thermal`) and
//! a documented point-source sensitivity (`crate::bands`). The headline
//! results pinned in the tests:
//!
//! 1. The maximum detectable distance grows with assumed mass (larger radius →
//!    brighter → detectable farther).
//! 2. SO reaches farther than ACT for the *same cold body*, because SO's deep
//!    map is fainter-limited.
//! 3. The detectable fraction of a reference mass–distance box lies in (0, 1)
//!    and is larger for SO than for ACT.

use p9_core::units::{au, earth_masses, Length, Mass};

use crate::bands::MmSensitivity;
use crate::thermal::P9Thermal;
use p9_core::analysis::thermal::max_detectable_distance as core_max_detectable_distance;

/// Maximum heliocentric distance at which a Planet Nine of `mass_earth`
/// emits at least `sens.flux_limit_mjy` at the sensitivity's band frequency,
/// as a dimension-checked [`Length`].
///
/// Flux falls monotonically as 1/d², so the reach is the unique root of
/// `F_ν(d) = flux_limit`. Solved by bisection over `[d_lo, d_hi]`; returns
/// `d_lo` if even the nearest distance is too faint and `d_hi` if the body is
/// still detectable at the far bracket.
pub fn max_detectable_distance(mass_earth: f64, sens: &MmSensitivity) -> Length {
    // Flux falls monotonically with distance, so negated flux plays the role
    // of a "magnitude" (dimmer = larger) for the shared core bisection.
    let (lo, hi) = (50.0_f64, 5000.0_f64);
    let reach = core_max_detectable_distance(lo, hi, -sens.flux_limit_mjy, |d| {
        -P9Thermal::new(mass_earth, d).flux_mjy(sens.nu_hz)
    });
    au(reach.unwrap_or(lo))
}

/// Expected detection significance (signal-to-noise ratio) of a Planet Nine of
/// `mass_earth` at `distance_au`, given a sensitivity whose `flux_limit_mjy`
/// represents the ~2σ (95%-confidence) point-source flux limit. Then the 1σ
/// flux noise is `flux_limit / 2` and
///
///   SNR = F_ν(d) / (flux_limit / 2).
///
/// A source exactly at the flux limit has SNR ≈ 2 (the 95%-CL detection
/// threshold).
pub fn detection_significance(mass_earth: f64, distance_au: f64, sens: &MmSensitivity) -> f64 {
    let flux = P9Thermal::new(mass_earth, distance_au).flux_mjy(sens.nu_hz);
    let sigma = sens.flux_limit_mjy / 2.0;
    flux / sigma
}

/// A rectangular reference parameter box in (mass, distance) over which the
/// detectable fraction is evaluated.
#[derive(Debug, Clone, Copy)]
pub struct ReferenceBox {
    pub mass_min_earth: f64,
    pub mass_max_earth: f64,
    pub dist_min_au: f64,
    pub dist_max_au: f64,
}

impl ReferenceBox {
    /// The plausible Planet Nine box used for the forecast: 5–10 M⊕ and
    /// 400–1100 AU, bracketing the Brown & Batygin posterior and extending out
    /// past the ACT/SO reach so that neither survey covers the full box (the
    /// detectable fraction stays strictly in (0, 1) for both).
    pub fn nominal() -> Self {
        Self {
            mass_min_earth: 5.0,
            mass_max_earth: 10.0,
            dist_min_au: 400.0,
            dist_max_au: 1100.0,
        }
    }

    /// Fraction of the box a survey detects: the measure of (mass, distance)
    /// points whose flux exceeds the sensitivity threshold, by a uniform grid
    /// quadrature. In (0, 1) for a sensitivity that detects part but not all
    /// of the box.
    /// Lower mass bound as a dimension-checked [`Mass`] (Earth-mass storage).
    pub fn mass_min(&self) -> Mass {
        earth_masses(self.mass_min_earth)
    }
    /// Upper mass bound as a dimension-checked [`Mass`] (Earth-mass storage).
    pub fn mass_max(&self) -> Mass {
        earth_masses(self.mass_max_earth)
    }
    /// Lower distance bound as a dimension-checked [`Length`] (AU storage).
    pub fn dist_min(&self) -> Length {
        au(self.dist_min_au)
    }
    /// Upper distance bound as a dimension-checked [`Length`] (AU storage).
    pub fn dist_max(&self) -> Length {
        au(self.dist_max_au)
    }

    pub fn detectable_fraction(&self, sens: &MmSensitivity) -> f64 {
        let n = 200usize;
        let mut detected = 0usize;
        for im in 0..n {
            let m = self.mass_min_earth
                + (self.mass_max_earth - self.mass_min_earth) * (im as f64 + 0.5) / n as f64;
            for id in 0..n {
                let d = self.dist_min_au
                    + (self.dist_max_au - self.dist_min_au) * (id as f64 + 0.5) / n as f64;
                if P9Thermal::new(m, d).flux_mjy(sens.nu_hz) >= sens.flux_limit_mjy {
                    detected += 1;
                }
            }
        }
        detected as f64 / (n * n) as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bands::{ACT_SENSITIVITY, SO_SENSITIVITY, SO_SENSITIVITY_SHALLOW};

    #[test]
    fn typed_boundary_matches_f64() {
        use approx::assert_relative_eq;
        use uom::si::length::astronomical_unit;
        use uom::si::mass::kilogram;

        // The bisection root recomputed inline must equal the typed reach.
        // The reach is defined by F(reach) == flux limit; check the
        // shared-core bisection lands on that root.
        let reach = max_detectable_distance(5.0, &SO_SENSITIVITY);
        let d = reach.get::<astronomical_unit>();
        assert!((50.0..5000.0).contains(&d), "reach = {d}");
        assert_relative_eq!(
            P9Thermal::new(5.0, d).flux_mjy(SO_SENSITIVITY.nu_hz),
            SO_SENSITIVITY.flux_limit_mjy,
            max_relative = 1e-6
        );
        let bx = ReferenceBox::nominal();
        assert_relative_eq!(
            bx.dist_min().get::<astronomical_unit>(),
            bx.dist_min_au,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            bx.mass_max().get::<kilogram>(),
            earth_masses(bx.mass_max_earth).get::<kilogram>(),
            max_relative = 1e-12
        );
    }

    #[test]
    fn reach_increases_with_mass() {
        // A heavier P9 has a larger radius, is brighter, and is detectable to a
        // greater distance — for both surveys.
        for sens in [&ACT_SENSITIVITY, &SO_SENSITIVITY] {
            let d5 = (max_detectable_distance(5.0, sens) / au(1.0)).value;
            let d10 = (max_detectable_distance(10.0, sens) / au(1.0)).value;
            let d15 = (max_detectable_distance(15.0, sens) / au(1.0)).value;
            assert!(
                d5 < d10 && d10 < d15,
                "{}: reach must rise with mass: {d5:.0} < {d10:.0} < {d15:.0}",
                sens.survey
            );
        }
    }

    #[test]
    fn so_reaches_farther_than_act_for_same_body() {
        // The headline comparison: for the same cold body SO's deeper map
        // detects it to a greater heliocentric distance than ACT.
        for &m in &[5.0, 10.0, 15.0] {
            let so = (max_detectable_distance(m, &SO_SENSITIVITY) / au(1.0)).value;
            let act = (max_detectable_distance(m, &ACT_SENSITIVITY) / au(1.0)).value;
            assert!(
                so > act,
                "{m} M⊕: SO reach {so:.0} AU should exceed ACT {act:.0} AU"
            );
        }
    }

    #[test]
    fn act_reach_matches_published_5me_range() {
        // Naess et al. (2021): a 5 M⊕ P9 has an ACT detection limit of
        // 325–625 AU depending on sky location; our median-depth (8 mJy)
        // reach lands inside that band.
        let act5 = (max_detectable_distance(5.0, &ACT_SENSITIVITY) / au(1.0)).value;
        assert!(
            (325.0..=625.0).contains(&act5),
            "ACT 5 M⊕ reach {act5:.0} AU should be within published 325–625"
        );
    }

    #[test]
    fn so_reach_matches_published_5me_band() {
        // SO forecast: 5 M⊕ detectable from 500 AU (shallow) to 900 AU (deep).
        let deep = (max_detectable_distance(5.0, &SO_SENSITIVITY) / au(1.0)).value;
        let shallow = (max_detectable_distance(5.0, &SO_SENSITIVITY_SHALLOW) / au(1.0)).value;
        assert!(
            (820.0..=980.0).contains(&deep),
            "SO deep 5 M⊕ reach {deep:.0} AU near published ~900"
        );
        assert!(
            (450.0..=560.0).contains(&shallow),
            "SO shallow 5 M⊕ reach {shallow:.0} AU near published ~500"
        );
    }

    #[test]
    fn significance_scales_inverse_distance_squared() {
        // SNR ∝ flux ∝ 1/d², so doubling distance drops SNR by 4×.
        let near = detection_significance(5.0, 400.0, &SO_SENSITIVITY);
        let far = detection_significance(5.0, 800.0, &SO_SENSITIVITY);
        let ratio = near / far;
        assert!((3.9..=4.1).contains(&ratio), "SNR ratio {ratio:.3} ≈ 4");
        // At the flux limit the SNR is ~2 (the 95%-CL threshold).
        let edge_d = (max_detectable_distance(5.0, &SO_SENSITIVITY) / au(1.0)).value;
        let snr_edge = detection_significance(5.0, edge_d, &SO_SENSITIVITY);
        assert!((snr_edge - 2.0).abs() < 0.05, "edge SNR {snr_edge:.3} ≈ 2");
    }

    #[test]
    fn detectable_fraction_in_unit_interval_and_so_beats_act() {
        let bx = ReferenceBox::nominal();
        let so = bx.detectable_fraction(&SO_SENSITIVITY);
        let act = bx.detectable_fraction(&ACT_SENSITIVITY);
        assert!((0.0..1.0).contains(&so), "SO fraction {so:.3} in (0,1)");
        assert!((0.0..1.0).contains(&act), "ACT fraction {act:.3} in (0,1)");
        assert!(so > 0.0 && act > 0.0, "both surveys detect part of the box");
        assert!(
            so > act,
            "SO fraction {so:.3} should exceed ACT {act:.3} over the box"
        );
    }
}
