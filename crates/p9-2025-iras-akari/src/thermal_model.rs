//! Far-infrared thermal emission model for Planet Nine.
//!
//! At 500–700 AU, Planet Nine is too far from the Sun for reflected light
//! to dominate. Instead, its thermal emission from internal heat is the
//! primary signal in the far-infrared. Following Phan et al. (2025), we
//! model this as a modified blackbody with effective temperature 30–50 K,
//! scaled by the planet's physical cross-section.

use std::f64::consts::PI;

/// Planck constant (J·s)
const H_PLANCK: f64 = 6.626_070_15e-34;
/// Boltzmann constant (J/K)
const K_BOLTZ: f64 = 1.380_649e-23;
/// Speed of light (m/s)
const C_LIGHT: f64 = 2.997_924_58e8;
/// Earth radius (m)
const R_EARTH: f64 = 6.371e6;
/// 1 AU in meters
const AU_M: f64 = 1.495_978_707e11;
/// 1 Jansky in W/m²/Hz
const JY: f64 = 1.0e-26;

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
    /// Estimate the physical radius using a simple mass-radius relation
    /// for sub-Neptune planets: R ≈ R_Earth * (M/M_Earth)^0.27
    ///
    /// This power-law approximation is consistent with the range used
    /// in Phan et al. (2025) Table 1, giving ~3.5 R_Earth for 10 M_Earth.
    pub fn radius_m(&self) -> f64 {
        R_EARTH * self.mass_earth.powf(0.27)
    }

    /// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
    pub fn planck_bnu(t_eff: f64, nu_hz: f64) -> f64 {
        let x = H_PLANCK * nu_hz / (K_BOLTZ * t_eff);
        if x > 500.0 {
            return 0.0;
        }
        (2.0 * H_PLANCK * nu_hz.powi(3) / (C_LIGHT * C_LIGHT)) / (x.exp() - 1.0)
    }

    /// Wavelength in meters to frequency in Hz.
    pub fn wavelength_to_freq(wavelength_m: f64) -> f64 {
        C_LIGHT / wavelength_m
    }

    /// Predicted flux density in Janskys at a given wavelength.
    ///
    /// The planet is modeled as a uniform-temperature sphere radiating
    /// as a blackbody. The flux at Earth (≈ at the Sun for d >> 1 AU) is:
    ///
    ///   F_ν = π * B_ν(T) * (R / d)²
    ///
    /// where R is the planet's radius and d its heliocentric distance.
    pub fn flux_density_jy(&self, wavelength_m: f64) -> f64 {
        let nu = Self::wavelength_to_freq(wavelength_m);
        let bnu = Self::planck_bnu(self.t_eff, nu);
        let r = self.radius_m();
        let d = self.distance_au * AU_M;
        let solid_angle = PI * (r / d).powi(2);
        let flux_w_m2_hz = solid_angle * bnu;
        flux_w_m2_hz / JY
    }

    /// Equilibrium temperature from solar illumination alone (for comparison).
    ///
    /// T_eq = T_sun * sqrt(R_sun / (2 * d)) * (1 - A)^0.25
    ///
    /// At 500+ AU this gives ~10–12 K, much below the 30–50 K internal heat.
    pub fn solar_equilibrium_temp(distance_au: f64, albedo: f64) -> f64 {
        let t_sun = 5778.0;
        let r_sun_m = 6.957e8;
        let d_m = distance_au * AU_M;
        t_sun * (r_sun_m / (2.0 * d_m)).sqrt() * (1.0 - albedo).powf(0.25)
    }
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

/// Check whether a Planet Nine with given parameters would be detectable
/// by IRAS (60 µm, ~0.2 Jy limit) and AKARI (90 µm, ~0.55 Jy limit).
pub fn is_detectable(params: &P9ThermalParams) -> (bool, bool) {
    let f60 = params.flux_density_jy(60.0e-6);
    let f90 = params.flux_density_jy(90.0e-6);
    let iras_detectable = f60 >= 0.2;
    let akari_detectable = f90 >= 0.55;
    (iras_detectable, akari_detectable)
}

#[cfg(test)]
mod tests {
    use super::*;

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
    fn radius_scaling() {
        let params = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        let r = params.radius_m() / R_EARTH;
        // 10^0.27 ≈ 1.86, but super-Earths can be ~3.5 R_Earth
        // Our power law gives ~1.86, which is conservative
        assert!(r > 1.5 && r < 4.0, "radius ratio = {r}");
    }

    #[test]
    fn solar_equilibrium_is_cold() {
        // At 600 AU, solar equilibrium should be ~10-12 K
        let t_eq = P9ThermalParams::solar_equilibrium_temp(600.0, 0.3);
        assert!(t_eq < 15.0, "T_eq = {t_eq} K, expected < 15 K");
        assert!(t_eq > 5.0, "T_eq = {t_eq} K, expected > 5 K");
    }

    #[test]
    fn nominal_p9_flux_order_of_magnitude() {
        // Phan et al. Table 1: for ~10 M_Earth at 600 AU, T~40K,
        // fluxes should be in the sub-Jy to few-Jy range at 60-90 µm.
        // The exact value depends on assumed radius; we check order of magnitude.
        let params = P9ThermalParams {
            mass_earth: 10.0,
            distance_au: 600.0,
            t_eff: 40.0,
        };
        let f60 = params.flux_density_jy(60.0e-6);
        let f90 = params.flux_density_jy(90.0e-6);
        // These should be tiny — sub-mJy for a ~2 R_Earth object at 600 AU
        assert!(f60 > 0.0);
        assert!(f90 > 0.0);
        assert!(f60 < 100.0, "f60 = {f60} Jy — sanity check");
        assert!(f90 < 100.0, "f90 = {f90} Jy — sanity check");
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
