//! WISE W1 (3.4 µm) detectability model for Planet Nine.
//!
//! Meisner et al. (2018) searched the WISE/NEOWISE data with a shift-and-stack
//! technique for a thermal-infrared point source. For a cold, distant planet
//! the relevant band is W1 (3.4 µm), with a co-added single-frame depth of
//! W1 ≈ 16.5 (Vega), from `p9_core::analysis::surveys` ("WISE W1").
//!
//! Two emission channels contribute at 3.4 µm and the total flux is their sum:
//!
//! 1. **Thermal emission.** The planet radiates as a grey body at an effective
//!    temperature set by the larger of (a) solar heating and (b) an internal
//!    luminosity floor. Flux ∝ (R/d)² · B_ν(T), with B_ν the Planck function.
//!
//! 2. **Reflected sunlight.** The planet reflects the 3.4 µm solar flux; this
//!    scales as 1/(r²Δ²) — once for Sun→planet, once for planet→observer.
//!
//! **Honest residual / physics note.** At W1 (3.4 µm) a cold (≈40 K) ice
//! giant sits *extremely* far down the Wien tail: hν/kT ≈ 107 at 40 K, so its
//! blackbody thermal flux at 3.4 µm is ~10⁻⁵⁵ W/m²/Hz/sr and utterly
//! negligible. The thermal channel only overtakes reflected sunlight above
//! ≈175 K (at 15 M⊕, 400 AU). So for the cold Planet Nine the paper targets, the
//! WISE W1 detectability is in practice set by **reflected sunlight**, and
//! this crate's "thermal-IR" exclusion is, at W1, a reflected-light exclusion
//! — the thermal term is retained (and would dominate for an implausibly warm
//! body) but contributes nothing at 40 K. The W2 (4.6 µm) band of the actual
//! Meisner et al. search is somewhat more thermally favorable; we model the
//! deeper, bluer W1 band that the shared survey table pins.
//!
//! Both channels are converted to a Vega-system W1 magnitude with the WISE W1
//! zero point F_ν0 = 309.540 Jy (Wright et al. 2010, Table 1). The radius law
//! is the shared Neptune-anchored mass-radius relation in
//! `p9_core::analysis::photometry` (R ≈ 3.34 R⊕ at 10 M⊕).

use p9_core::analysis::photometry::{mass_radius_neptunian, ALBEDO_NEPTUNE};
use p9_core::analysis::thermal::{
    effective_temp, flux_to_magnitude, reflected_flux_jy, solar_equilibrium_temp, thermal_flux_jy,
    C_LIGHT,
};
use p9_core::units::{au, earth_masses, meters, Length, Mass};

/// Earth radius (m).
const R_EARTH_M: f64 = 6.371e6;

/// WISE W1 effective wavelength in meters (3.3526 µm; Wright et al. 2010).
pub const W1_WAVELENGTH_M: f64 = 3.3526e-6;

/// WISE W1 Vega zero-point flux density in Jansky: a 0.0-mag (Vega) source
/// has F_ν = 309.540 Jy at W1 (Wright et al. 2010, Table 1). A source of
/// flux F has W1 = −2.5 log10(F / F_ν0).
pub const W1_ZERO_POINT_JY: f64 = 309.540;

/// Internal-luminosity temperature floor (K). A cold ice giant retains
/// formation/radiogenic heat; Fortney et al. (2016) and the Meisner et al.
/// (2018) discussion adopt an effective temperature of ~30–40 K for a
/// several-Earth-mass Planet Nine essentially independent of heliocentric
/// distance beyond a few hundred AU. We take 40 K as the nominal internal
/// floor and let solar heating raise it when it dominates (it does not at
/// these distances).
pub const INTERNAL_TEMP_FLOOR_K: f64 = 40.0;

/// Physical parameters of a candidate Planet Nine for the W1 model.
#[derive(Debug, Clone, Copy)]
pub struct P9Thermal {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Heliocentric distance in AU.
    pub distance_au: f64,
    /// Geometric albedo for the reflected-light channel.
    pub albedo: f64,
    /// Internal-luminosity effective-temperature floor in Kelvin.
    pub internal_temp_k: f64,
}

/// W1-band (3.35 µm) geometric albedo range for an ice-giant atmosphere.
/// The 3.3 µm region sits in the strong CH₄ ν₃ absorption band, where
/// Uranus/Neptune reflectance is a few percent — NOT the V-band 0.41 a
/// previous version applied (which overstated the W1 reflected reach by
/// ~2.5 mag for a Neptune-like atmosphere). The bright end covers CH₄-poor
/// or high cloud decks.
pub const W1_ALBEDO_DARK: f64 = 0.03;
/// Adopted default W1-band albedo (labelled model choice, geometric mean of
/// the dark CH₄-absorbed and bright cloud-deck ends).
pub const W1_ALBEDO_DEFAULT: f64 = 0.10;
/// Bright end of the W1 albedo range (the V-band value, as the optimistic
/// bound for a CH₄-poor atmosphere).
pub const W1_ALBEDO_BRIGHT: f64 = ALBEDO_NEPTUNE;

impl P9Thermal {
    /// A Planet Nine at the given mass and heliocentric distance with the
    /// default W1-band albedo and internal-temperature floor.
    pub fn new(mass_earth: f64, distance_au: f64) -> Self {
        Self {
            mass_earth,
            distance_au,
            albedo: W1_ALBEDO_DEFAULT,
            internal_temp_k: INTERNAL_TEMP_FLOOR_K,
        }
    }

    /// Same body with an explicit W1 albedo (sensitivity scans).
    pub fn with_albedo(mut self, albedo: f64) -> Self {
        self.albedo = albedo;
        self
    }

    /// Physical radius in meters (Neptune-anchored mass-radius relation).
    pub fn radius_m(&self) -> f64 {
        mass_radius_neptunian(self.mass_earth) * R_EARTH_M
    }

    /// Mass as a typed [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }

    /// Heliocentric distance as a typed [`Length`].
    pub fn distance(&self) -> Length {
        au(self.distance_au)
    }

    /// Physical radius as a typed [`Length`].
    pub fn radius(&self) -> Length {
        meters(self.radius_m())
    }

    /// Solar equilibrium temperature (K) for a fast rotator radiating from
    /// its full surface: T_eq = T_sun √(R_sun / 2d) (1 − A)^¼.
    pub fn solar_equilibrium_temp(&self) -> f64 {
        solar_equilibrium_temp(self.distance_au, self.albedo)
    }

    /// Effective temperature: the larger of the internal-heat floor and the
    /// solar-equilibrium temperature. Beyond a few hundred AU the solar term
    /// is ~10 K, so the internal floor dominates.
    pub fn effective_temp(&self) -> f64 {
        effective_temp(self.distance_au, self.albedo, self.internal_temp_k)
    }

    /// Thermal flux density at W1 (Jy), a grey-body emitter of unit
    /// emissivity: F_ν = π B_ν(T) (R/d)².
    pub fn thermal_flux_w1_jy(&self) -> f64 {
        let nu = C_LIGHT / W1_WAVELENGTH_M;
        thermal_flux_jy(self.effective_temp(), self.radius_m(), self.distance_au, nu)
    }

    /// Reflected-sunlight flux density at W1 (Jy), observed near opposition
    /// (Δ = d − 1 AU). The Sun is a blackbody at `T_SUN`; the solar W1 flux
    /// at the planet is π B_ν(T_sun) (R_sun / r)². The planet intercepts a
    /// disc πR_p², reflects a fraction `albedo` into ~2π sr (Lambert factor
    /// absorbed into the geometric albedo), and that reradiated flux falls
    /// as 1/Δ² to the observer:
    ///
    ///   F_refl = albedo · [π B_ν(T_sun) (R_sun/r)²] · (R_p² / (4 Δ²)).
    pub fn reflected_flux_w1_jy(&self) -> f64 {
        let nu = C_LIGHT / W1_WAVELENGTH_M;
        reflected_flux_jy(self.albedo, self.radius_m(), self.distance_au, nu)
    }

    /// Total W1 flux density (Jy): thermal + reflected.
    pub fn total_flux_w1_jy(&self) -> f64 {
        self.thermal_flux_w1_jy() + self.reflected_flux_w1_jy()
    }

    /// Apparent W1 (Vega) magnitude from the total flux density. Fainter
    /// (larger) numbers for smaller flux; +∞ when the flux is zero.
    pub fn w1_magnitude(&self) -> f64 {
        flux_to_magnitude(self.total_flux_w1_jy(), W1_ZERO_POINT_JY)
    }
}

/// Apparent W1 magnitude for a default-albedo P9 of `mass_earth` at
/// heliocentric distance `distance_au`.
pub fn w1_magnitude(mass_earth: f64, distance_au: f64) -> f64 {
    P9Thermal::new(mass_earth, distance_au).w1_magnitude()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn typed_accessors_match_f64() {
        use uom::si::length::{astronomical_unit, meter};
        use uom::si::mass::kilogram;
        let p = P9Thermal::new(10.0, 500.0);
        assert!((p.distance().get::<astronomical_unit>() - 500.0).abs() < 1e-9);
        assert!((p.radius().get::<meter>() - p.radius_m()).abs() < 1e-6);
        assert!((p.mass().get::<kilogram>() - 10.0 * p9_core::units::EARTH_MASS_KG).abs() < 1e12);
    }

    #[test]
    fn radius_matches_ice_giant_models() {
        let r = P9Thermal::new(10.0, 500.0).radius_m() / R_EARTH_M;
        assert!(r > 3.0 && r < 3.6, "R(10 M⊕) = {r:.2} R⊕");
    }

    #[test]
    fn solar_heating_is_cold_and_floor_dominates() {
        // At hundreds of AU solar equilibrium is ~10 K, below the internal
        // floor, so the effective temperature is the floor.
        let p = P9Thermal::new(10.0, 500.0);
        assert!(
            p.solar_equilibrium_temp() < 15.0,
            "T_eq = {:.1} K",
            p.solar_equilibrium_temp()
        );
        assert!((p.effective_temp() - INTERNAL_TEMP_FLOOR_K).abs() < 1e-9);
    }

    #[test]
    fn magnitude_brightens_with_mass() {
        // More massive -> larger radius -> more flux -> smaller (brighter) mag.
        let light = w1_magnitude(5.0, 500.0);
        let heavy = w1_magnitude(15.0, 500.0);
        assert!(
            heavy < light,
            "heavy {heavy:.2} should be brighter than light {light:.2}"
        );
    }

    #[test]
    fn magnitude_faints_with_distance() {
        let near = w1_magnitude(10.0, 400.0);
        let far = w1_magnitude(10.0, 800.0);
        assert!(
            far > near,
            "far {far:.2} should be fainter than near {near:.2}"
        );
    }

    #[test]
    fn reflected_dominates_at_cold_temperatures_but_thermal_wins_when_warm() {
        // The honest W1 physics: at the 40 K floor the thermal 3.4 µm flux is
        // negligible and reflected sunlight dominates...
        let cold = P9Thermal::new(15.0, 400.0);
        assert!(
            cold.reflected_flux_w1_jy() > 1e6 * cold.thermal_flux_w1_jy(),
            "at 40 K reflected {:.3e} should swamp thermal {:.3e} Jy",
            cold.reflected_flux_w1_jy(),
            cold.thermal_flux_w1_jy()
        );
        // ...but only a hot (≳175 K) body crosses over to thermal at 3.4 µm,
        // far above any plausible Planet Nine temperature — underscoring how
        // thermally hostile W1 is for a cold planet.
        let warm = P9Thermal {
            internal_temp_k: 200.0,
            ..cold
        };
        assert!(
            warm.thermal_flux_w1_jy() > warm.reflected_flux_w1_jy(),
            "at 120 K thermal {:.3e} should exceed reflected {:.3e} Jy",
            warm.thermal_flux_w1_jy(),
            warm.reflected_flux_w1_jy()
        );
    }

    #[test]
    fn zero_point_gives_zero_magnitude() {
        // A source emitting exactly the W1 zero-point flux is W1 = 0.
        let flux = W1_ZERO_POINT_JY;
        let mag = -2.5 * (flux / W1_ZERO_POINT_JY).log10();
        assert!(mag.abs() < 1e-12);
    }

    #[test]
    fn thermal_flux_steeply_temperature_dependent() {
        // On the Wien tail at 3.4 µm a 50 K planet is vastly brighter than a
        // 30 K one of the same size/distance.
        let cold = P9Thermal {
            internal_temp_k: 30.0,
            ..P9Thermal::new(10.0, 500.0)
        };
        let warm = P9Thermal {
            internal_temp_k: 50.0,
            ..P9Thermal::new(10.0, 500.0)
        };
        let ratio = warm.thermal_flux_w1_jy() / cold.thermal_flux_w1_jy();
        assert!(ratio > 100.0, "warm/cold thermal ratio = {ratio:.1}");
    }
}
