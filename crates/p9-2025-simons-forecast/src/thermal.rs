//! Millimeter-wave thermal emission of a cold Planet Nine.
//!
//! The Simons Observatory (SO) and ACT searches are *millimeter* searches, in
//! contrast to the near-/mid-IR WISE search (`p9-2018-wise-search`). At
//! mm wavelengths the relevant physics is qualitatively cleaner: for a Planet
//! Nine with an effective temperature `T ≈ 30–50 K`, a 280 GHz (λ ≈ 1.07 mm)
//! photon has
//!
//!   x = hν/kT ≈ (6.6e-34 · 2.8e11) / (1.38e-23 · 40) ≈ 0.34,
//!
//! so the source sits in the **Rayleigh–Jeans** regime (x ≪ 1), where the
//! Planck function reduces to
//!
//!   B_ν(T) ≈ 2 ν² k T / c²            (flux ∝ ν² T).
//!
//! This is the opposite of the WISE W1 case (x ≈ 107 at 3.4 µm, deep on the
//! Wien tail, where the thermal flux is ~10⁻⁵⁵ W/m²/Hz/sr and *reflected*
//! sunlight dominates). At mm wavelengths thermal emission completely
//! dominates reflected sunlight for a cold body, and the flux is a smooth
//! linear function of temperature rather than an exponential one — which is
//! exactly why a mm CMB survey is a powerful Planet Nine probe.
//!
//! The flux density of a body of radius `R` at heliocentric distance `d`
//! (observed near opposition, Δ ≈ d) radiating as a blackbody of unit
//! emissivity is
//!
//!   F_ν = π B_ν(T) (R / d)²,
//!
//! falling as 1/d². Radius from mass uses the shared Neptune-anchored
//! mass–radius relation in `p9_core::analysis::photometry`, and the effective
//! temperature is the larger of an internal-luminosity floor and the solar
//! equilibrium temperature (the latter is ~10 K beyond a few hundred AU, so
//! the floor dominates) — the same temperature model as the WISE crate, so the
//! two reproductions share a physical body.

use std::f64::consts::PI;

use p9_core::analysis::photometry::mass_radius_neptunian;

/// Planck constant (J·s).
const H_PLANCK: f64 = 6.626_070_15e-34;
/// Boltzmann constant (J/K).
const K_BOLTZ: f64 = 1.380_649e-23;
/// Speed of light (m/s).
const C_LIGHT: f64 = 2.997_924_58e8;
/// Earth radius (m).
const R_EARTH_M: f64 = 6.371e6;
/// 1 AU in meters.
const AU_M: f64 = 1.495_978_707e11;
/// Solar radius (m).
const R_SUN_M: f64 = 6.957e8;
/// Solar effective temperature (K).
const T_SUN: f64 = 5778.0;
/// 1 Jansky in W/m²/Hz.
const JY: f64 = 1.0e-26;

/// Neptune-like geometric/Bond albedo, used only for the solar-equilibrium
/// temperature term.
pub const DEFAULT_ALBEDO: f64 = 0.41;

/// Internal-luminosity effective-temperature floor (K). A several-Earth-mass
/// ice giant retains formation/radiogenic heat; Fortney et al. (2016) and the
/// ACT/SO discussions adopt an effective temperature of ~30–50 K essentially
/// independent of heliocentric distance beyond a few hundred AU. We take 40 K
/// as the nominal floor — identical to the WISE crate so both reproductions
/// model the same body.
pub const INTERNAL_TEMP_FLOOR_K: f64 = 40.0;

/// Physical parameters of a candidate Planet Nine for the mm thermal model.
#[derive(Debug, Clone, Copy)]
pub struct P9Thermal {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Heliocentric distance in AU.
    pub distance_au: f64,
    /// Geometric albedo (solar-equilibrium term only).
    pub albedo: f64,
    /// Internal-luminosity effective-temperature floor in Kelvin.
    pub internal_temp_k: f64,
}

impl P9Thermal {
    /// A Planet Nine at the given mass and heliocentric distance with the
    /// default albedo and internal-temperature floor.
    pub fn new(mass_earth: f64, distance_au: f64) -> Self {
        Self {
            mass_earth,
            distance_au,
            albedo: DEFAULT_ALBEDO,
            internal_temp_k: INTERNAL_TEMP_FLOOR_K,
        }
    }

    /// Physical radius in meters (Neptune-anchored mass–radius relation).
    pub fn radius_m(&self) -> f64 {
        mass_radius_neptunian(self.mass_earth) * R_EARTH_M
    }

    /// Solar equilibrium temperature (K) for a fast rotator radiating from its
    /// full surface: T_eq = T_sun √(R_sun / 2d) (1 − A)^¼.
    pub fn solar_equilibrium_temp(&self) -> f64 {
        let d_m = self.distance_au * AU_M;
        T_SUN * (R_SUN_M / (2.0 * d_m)).sqrt() * (1.0 - self.albedo).powf(0.25)
    }

    /// Effective temperature: the larger of the internal-heat floor and the
    /// solar-equilibrium temperature.
    pub fn effective_temp(&self) -> f64 {
        self.solar_equilibrium_temp().max(self.internal_temp_k)
    }

    /// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
    pub fn planck_bnu(t: f64, nu_hz: f64) -> f64 {
        let x = H_PLANCK * nu_hz / (K_BOLTZ * t);
        if x > 700.0 {
            return 0.0;
        }
        (2.0 * H_PLANCK * nu_hz.powi(3) / (C_LIGHT * C_LIGHT)) / (x.exp() - 1.0)
    }

    /// Rayleigh–Jeans approximation to the Planck function (valid for x ≪ 1):
    ///
    ///   B_ν(T) ≈ 2 ν² k T / c².
    ///
    /// Exposed so tests can pin that the exact Planck flux follows ∝ ν² T at
    /// the SO/ACT bands.
    pub fn rayleigh_jeans_bnu(t: f64, nu_hz: f64) -> f64 {
        2.0 * nu_hz * nu_hz * K_BOLTZ * t / (C_LIGHT * C_LIGHT)
    }

    /// Thermal flux density (Jy) at frequency `nu_hz`, a blackbody emitter of
    /// unit emissivity observed near opposition (Δ ≈ d):
    ///
    ///   F_ν = π B_ν(T) (R / d)².
    pub fn flux_jy(&self, nu_hz: f64) -> f64 {
        let bnu = Self::planck_bnu(self.effective_temp(), nu_hz);
        let r = self.radius_m();
        let d = self.distance_au * AU_M;
        let solid_angle = PI * (r / d).powi(2);
        solid_angle * bnu / JY
    }

    /// Thermal flux density (mJy) at frequency `nu_hz`.
    pub fn flux_mjy(&self, nu_hz: f64) -> f64 {
        self.flux_jy(nu_hz) * 1000.0
    }
}

/// Convenience: thermal flux density (mJy) for a default-floor Planet Nine of
/// `mass_earth` at `distance_au`, observed at frequency `nu_ghz` (GHz).
pub fn flux_mjy(mass_earth: f64, distance_au: f64, nu_ghz: f64) -> f64 {
    P9Thermal::new(mass_earth, distance_au).flux_mjy(nu_ghz * 1e9)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::bands::SO_280;

    #[test]
    fn cold_body_is_in_rayleigh_jeans_regime_at_280ghz() {
        // x = hν/kT ≈ 0.34 at 280 GHz, 40 K — squarely Rayleigh–Jeans (x ≪ 1),
        // the opposite of WISE W1 (x ≈ 107). The exact Planck function and the
        // RJ approximation must agree to better than ~20% here.
        let nu = 280e9;
        let t = 40.0;
        let x = H_PLANCK * nu / (K_BOLTZ * t);
        assert!(x < 0.5, "x = {x:.3} should be ≪ 1 (Rayleigh–Jeans)");
        let exact = P9Thermal::planck_bnu(t, nu);
        let rj = P9Thermal::rayleigh_jeans_bnu(t, nu);
        let rel = (exact - rj).abs() / rj;
        assert!(rel < 0.2, "Planck vs RJ relative error {rel:.3} at 280 GHz");
    }

    #[test]
    fn flux_follows_rayleigh_jeans_nu_squared_t() {
        // In the RJ regime F_ν ∝ ν² T at fixed (R, d). Compare flux at 280 vs
        // 140 GHz: doubling ν should ~quadruple the flux.
        let p = P9Thermal::new(5.0, 500.0);
        let f_hi = p.flux_jy(280e9);
        let f_lo = p.flux_jy(140e9);
        let ratio = f_hi / f_lo;
        // Exact Planck deviates slightly below the pure ν² law at 280 GHz (the
        // x³/(eˣ−1) curvature), but stays close to 4×: a clear ν² signature.
        assert!(
            (3.5..=4.0).contains(&ratio),
            "F(280)/F(140) = {ratio:.3}, expected ≈4 (ν² scaling)"
        );

        // Approximately linear in temperature at fixed frequency (RJ): doubling
        // T roughly doubles the flux. The exact Planck function rises slightly
        // *faster* than linear here (the colder 30 K body is a touch further off
        // the RJ limit, x = 0.45, than the 60 K one), so the ratio sits just
        // above 2 — an honest, small departure from the pure RJ law.
        let cold = P9Thermal {
            internal_temp_k: 30.0,
            ..p
        };
        let warm = P9Thermal {
            internal_temp_k: 60.0,
            ..p
        };
        let tratio = warm.flux_jy(280e9) / cold.flux_jy(280e9);
        assert!(
            (2.0..=2.4).contains(&tratio),
            "F(60K)/F(30K) = {tratio:.3}, expected just above 2 (≈ linear in T)"
        );
    }

    #[test]
    fn flux_falls_as_inverse_distance_squared() {
        let near = P9Thermal::new(5.0, 500.0).flux_jy(SO_280);
        let far = P9Thermal::new(5.0, 1000.0).flux_jy(SO_280);
        let ratio = near / far;
        assert!(
            (3.95..=4.05).contains(&ratio),
            "F(500)/F(1000) = {ratio:.4}, expected 4 (1/d²)"
        );
    }

    #[test]
    fn flux_scale_matches_act_so_regime() {
        // A 5 M⊕ P9 at 500 AU emits ~14 mJy at 280 GHz: above ACT's 4–12 mJy
        // depth (detectable) and far above SO's deeper threshold.
        let f = flux_mjy(5.0, 500.0, 280.0);
        assert!((10.0..20.0).contains(&f), "F(5 M⊕, 500 AU) = {f:.1} mJy");
        // At 900 AU it has fallen to ~4 mJy — at the edge of even SO's reach.
        let f900 = flux_mjy(5.0, 900.0, 280.0);
        assert!(
            (3.0..6.0).contains(&f900),
            "F(5 M⊕, 900 AU) = {f900:.1} mJy"
        );
    }

    #[test]
    fn brighter_with_mass() {
        let light = flux_mjy(5.0, 500.0, 280.0);
        let heavy = flux_mjy(15.0, 500.0, 280.0);
        assert!(heavy > light, "heavier should be brighter");
    }

    #[test]
    fn solar_heating_negligible_floor_dominates() {
        let p = P9Thermal::new(5.0, 500.0);
        assert!(p.solar_equilibrium_temp() < 15.0);
        assert!((p.effective_temp() - INTERNAL_TEMP_FLOOR_K).abs() < 1e-9);
    }
}
