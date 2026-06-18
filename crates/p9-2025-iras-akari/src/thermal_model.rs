//! Far-infrared thermal emission model for Planet Nine.
//!
//! At 500–700 AU, Planet Nine is too far from the Sun for reflected light
//! to dominate. Instead, its thermal emission from internal heat is the
//! primary signal in the far-infrared. Following Phan et al. (2025), we
//! model this as a blackbody with effective temperature 30–50 K, scaled by
//! the planet's physical cross-section, with an optional grey-emissivity
//! hook (H₂ collision-induced absorption makes ice giants deviate from a
//! perfect blackbody near 40 K).
//!
//! Radius law: the Neptune-anchored mass-radius relation from
//! `p9_core::analysis::photometry::mass_radius_neptunian`
//! (R = 3.865 R⊕ · (M / 17.147 M⊕)^0.27, consistent with Fortney et al.
//! 2007 ice-giant models and the Chen & Kipping 2017 Neptunian branch),
//! giving ≈ 3.34 R⊕ at 10 M⊕. The previous local law R⊕·M^0.27 gave
//! 1.86 R⊕ at 10 M⊕ — a ~3.2× flux underestimate that pushed the paper's
//! whole favorable corner below the IRAS detection limit.

use p9_core::analysis::photometry::mass_radius_neptunian;
use p9_core::analysis::thermal::{planck_bnu, solar_equilibrium_temp, thermal_flux_jy, C_LIGHT};
use p9_core::units::{au, earth_masses, meters, Length, Mass};

use crate::survey_model::{AkariFisSurvey, IrasSurvey};

/// Earth radius (m)
const R_EARTH: f64 = 6.371e6;

/// Physical parameters for a Planet Nine thermal model.
#[derive(Debug, Clone, Copy)]
pub struct P9ThermalParams {
    /// Mass in Earth masses
    pub mass_earth: f64,
    /// Heliocentric distance in AU
    pub distance_au: f64,
    /// Effective temperature in Kelvin (internal heat)
    pub t_eff: f64,
}

impl P9ThermalParams {
    /// Physical radius from the Neptune-anchored mass-radius relation in
    /// `p9_core::analysis::photometry` (≈ 3.34 R⊕ at 10 M⊕; see module
    /// docs for sources).
    pub fn radius_m(&self) -> f64 {
        mass_radius_neptunian(self.mass_earth) * R_EARTH
    }

    /// Mass as a dimension-checked [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }

    /// Heliocentric distance as a dimension-checked [`Length`].
    pub fn distance(&self) -> Length {
        au(self.distance_au)
    }

    /// Physical radius as a dimension-checked [`Length`].
    pub fn radius(&self) -> Length {
        meters(self.radius_m())
    }

    /// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
    pub fn planck_bnu(t_eff: f64, nu_hz: f64) -> f64 {
        planck_bnu(t_eff, nu_hz)
    }

    /// Wavelength in meters to frequency in Hz.
    pub fn wavelength_to_freq(wavelength_m: f64) -> f64 {
        C_LIGHT / wavelength_m
    }

    /// Predicted flux density in Janskys at a given wavelength for a grey
    /// emitter of the given `emissivity` (1.0 = blackbody):
    ///
    ///   F_ν = ε · π · B_ν(T) · (R / d)²
    ///
    /// The emissivity hook accommodates H₂ collision-induced-absorption
    /// deviations from a blackbody expected for a cold (≈40 K) ice giant.
    pub fn flux_density_jy_with_emissivity(&self, wavelength_m: f64, emissivity: f64) -> f64 {
        let nu = Self::wavelength_to_freq(wavelength_m);
        emissivity * thermal_flux_jy(self.t_eff, self.radius_m(), self.distance_au, nu)
    }

    /// Predicted blackbody (ε = 1) flux density in Janskys.
    pub fn flux_density_jy(&self, wavelength_m: f64) -> f64 {
        self.flux_density_jy_with_emissivity(wavelength_m, 1.0)
    }

    /// Equilibrium temperature from solar illumination alone (for comparison).
    ///
    /// T_eq = T_sun * sqrt(R_sun / (2 * d)) * (1 - A)^0.25
    ///
    /// At 500+ AU this gives ~10–12 K, much below the 30–50 K internal heat.
    pub fn solar_equilibrium_temp(distance_au: f64, albedo: f64) -> f64 {
        solar_equilibrium_temp(distance_au, albedo)
    }
}

/// Blackbody flux-density ratio F_60µm / F_90µm at temperature `t_eff` —
/// the radius and distance cancel, so the ratio is a pure temperature
/// diagnostic (0.23 at 30 K to 0.66 at 50 K). Used to derive the
/// candidate-selection flux-ratio window from the thermal model itself.
pub fn flux_ratio_60_90(t_eff: f64) -> f64 {
    let nu60 = P9ThermalParams::wavelength_to_freq(60.0e-6);
    let nu90 = P9ThermalParams::wavelength_to_freq(90.0e-6);
    P9ThermalParams::planck_bnu(t_eff, nu60) / P9ThermalParams::planck_bnu(t_eff, nu90)
}

/// Compute flux densities across a grid of masses, distances, and temperatures.
///
/// Returns a Vec of (mass, distance, t_eff, flux_60um_jy, flux_90um_jy).
pub fn flux_grid(
    masses: &[f64],
    distances: &[f64],
    temperatures: &[f64],
) -> Vec<(f64, f64, f64, f64, f64)> {
    let wavelength_60um = 60.0e-6;
    let wavelength_90um = 90.0e-6;

    let mut results = Vec::new();
    for &m in masses {
        for &d in distances {
            for &t in temperatures {
                let params = P9ThermalParams {
                    mass_earth: m,
                    distance_au: d,
                    t_eff: t,
                };
                let f60 = params.flux_density_jy(wavelength_60um);
                let f90 = params.flux_density_jy(wavelength_90um);
                results.push((m, d, t, f60, f90));
            }
        }
    }
    results
}

/// Whether a Planet Nine with the given parameters would be catalogued by
/// the two surveys, using each survey's own band and sensitivity (no
/// hard-coded duplicates of the survey limits).
pub fn is_detectable(
    params: &P9ThermalParams,
    iras: &IrasSurvey,
    akari: &AkariFisSurvey,
) -> (bool, bool) {
    let f_iras = params.flux_density_jy(iras.wavelength_um * 1.0e-6);
    let f_akari = params.flux_density_jy(akari.wavelength_um * 1.0e-6);
    (iras.is_detectable(f_iras), akari.is_detectable(f_akari))
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::length::{astronomical_unit, meter};
    use uom::si::mass::kilogram;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let p = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        assert_relative_eq!(
            p.distance().get::<astronomical_unit>(),
            p.distance_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p.radius().get::<meter>(),
            p.radius_m(),
            max_relative = 1e-12
        );
        let earth_kg = earth_masses(1.0).get::<kilogram>();
        assert_relative_eq!(
            p.mass().get::<kilogram>() / earth_kg,
            p.mass_earth,
            max_relative = 1e-12
        );
    }

    #[test]
    fn planck_peak_shifts_with_temperature() {
        // Wien's displacement: peak frequency ∝ T
        // For T=40K, peak is at ~2.4 THz (~125 µm)
        let bnu_low = P9ThermalParams::planck_bnu(40.0, 1.0e12);
        let bnu_peak = P9ThermalParams::planck_bnu(40.0, 2.4e12);
        let bnu_high = P9ThermalParams::planck_bnu(40.0, 1.0e13);
        assert!(bnu_peak > bnu_low);
        assert!(bnu_peak > bnu_high);
    }

    #[test]
    fn flux_decreases_with_distance() {
        let near = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 500.0,
            t_eff: 40.0,
        };
        let far = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 700.0,
            t_eff: 40.0,
        };
        let f_near = near.flux_density_jy(60.0e-6);
        let f_far = far.flux_density_jy(60.0e-6);
        assert!(f_near > f_far);
        // Should scale as 1/d²
        let ratio = f_near / f_far;
        let expected_ratio = (700.0 / 500.0_f64).powi(2);
        assert!((ratio - expected_ratio).abs() / expected_ratio < 0.01);
    }

    #[test]
    fn flux_increases_with_mass() {
        let light = P9ThermalParams {
            mass_earth: 7.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        let heavy = P9ThermalParams {
            mass_earth: 17.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        assert!(heavy.flux_density_jy(60.0e-6) > light.flux_density_jy(60.0e-6));
    }

    #[test]
    fn radius_matches_ice_giant_models() {
        // Fortney et al. (2007): ~3.3–3.5 R⊕ at 10 M⊕ (the old local law
        // gave 1.86 R⊕, a ~3.2× flux underestimate).
        let params = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        let r = params.radius_m() / R_EARTH;
        assert!(r > 3.0 && r < 3.6, "radius ratio = {r}");
    }

    #[test]
    fn solar_equilibrium_is_cold() {
        // At 600 AU, solar equilibrium should be ~10-12 K
        let t_eq = P9ThermalParams::solar_equilibrium_temp(600.0, 0.3);
        assert!(t_eq < 15.0, "T_eq = {t_eq} K, expected < 15 K");
        assert!(t_eq > 5.0, "T_eq = {t_eq} K, expected > 5 K");
    }

    #[test]
    fn favorable_corner_crosses_survey_thresholds() {
        // The paper's premise: at the favorable corner of the grid
        // (10 M⊕, 500 AU, warm end T = 50 K) the predicted fluxes cross
        // both survey limits. Pinned values: F_60 ≈ 0.39 Jy ≥ 0.2 Jy
        // (IRAS), F_90 ≈ 0.59 Jy ≥ 0.55 Jy (AKARI).
        let params = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 500.0,
            t_eff: 50.0,
        };
        let f60 = params.flux_density_jy(60.0e-6);
        let f90 = params.flux_density_jy(90.0e-6);
        assert!((0.30..0.50).contains(&f60), "f60 = {f60:.3} Jy");
        assert!((0.50..0.70).contains(&f90), "f90 = {f90:.3} Jy");

        let (iras_ok, akari_ok) =
            is_detectable(&params, &IrasSurvey::default(), &AkariFisSurvey::default());
        assert!(iras_ok, "IRAS should detect the favorable corner");
        assert!(akari_ok, "AKARI should detect the favorable corner");
    }

    #[test]
    fn unfavorable_corner_below_iras_threshold() {
        // Cold, light, distant corner: 7 M⊕ at 700 AU, 30 K → F_60 ≈
        // 0.007 Jy, far below the 0.2 Jy IRAS limit (detectability is
        // marginal across the grid, as the paper notes).
        let params = P9ThermalParams {
            mass_earth: 7.0,
            distance_au: 700.0,
            t_eff: 30.0,
        };
        let f60 = params.flux_density_jy(60.0e-6);
        assert!(f60 < 0.02, "f60 = {f60:.4} Jy");
        let (iras_ok, _) =
            is_detectable(&params, &IrasSurvey::default(), &AkariFisSurvey::default());
        assert!(!iras_ok);
    }

    #[test]
    fn nominal_p9_flux_pinned() {
        // Mid-grid (10 M⊕, 600 AU, 40 K): F_60 ≈ 0.081 Jy.
        let params = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        let f60 = params.flux_density_jy(60.0e-6);
        assert!((0.06..0.11).contains(&f60), "f60 = {f60:.4} Jy");
    }

    #[test]
    fn emissivity_scales_flux_linearly() {
        let params = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        let full = params.flux_density_jy_with_emissivity(60.0e-6, 1.0);
        let grey = params.flux_density_jy_with_emissivity(60.0e-6, 0.7);
        assert!((grey / full - 0.7).abs() < 1e-12);
        assert!((params.flux_density_jy(60.0e-6) - full).abs() < 1e-30);
    }

    #[test]
    fn flux_ratio_matches_blackbody_range() {
        // 30–50 K blackbody: F60/F90 from ~0.23 to ~0.66, monotonic in T.
        let r30 = flux_ratio_60_90(30.0);
        let r50 = flux_ratio_60_90(50.0);
        assert!((r30 - 0.233).abs() < 0.02, "ratio(30K) = {r30:.3}");
        assert!((r50 - 0.66).abs() < 0.03, "ratio(50K) = {r50:.3}");
        assert!(r50 > r30);
    }

    #[test]
    fn hotter_planet_brighter_at_60um() {
        let cool = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 30.0,
        };
        let warm = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 50.0,
        };
        assert!(warm.flux_density_jy(60.0e-6) > cool.flux_density_jy(60.0e-6));
    }
}
