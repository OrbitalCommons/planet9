//! Thermal-evolution effective temperature and far-IR SED of Planet Nine.
//!
//! Fortney et al. (2016), "The Hunt for Planet Nine: Atmosphere, Spectra,
//! Evolution, and Detectability" (arXiv:1604.07424), evolve ice-giant interior
//! models forward to ~4.5 Gyr and find Planet Nine should retain enough
//! internal heat to sit at an **effective temperature of ~35–50 K for masses
//! of 5–20 M⊕** — far above the feeble solar-equilibrium temperature at
//! hundreds of AU. The body therefore radiates almost entirely in the
//! far-infrared / sub-mm, with a blackbody SED peaking at tens of microns, and
//! is essentially invisible in the near-IR WISE bands (W1 3.4 µm, W2 4.6 µm)
//! that probe the Wien tail.
//!
//! # What this crate computes (no hard-coded answers)
//!
//! The headline number is the effective temperature from an energy balance
//! between the internal luminosity `L_int` (a thermal-evolution model input)
//! and the absorbed sunlight `L_absorbed`:
//!
//! ```text
//!   T_eff = [ (L_int + L_absorbed) / (4 π R² σ) ]^{1/4}
//! ```
//!
//! where `R` is the Neptune-anchored mass-radius radius
//! ([`p9_core::analysis::photometry::mass_radius_neptunian`]), `σ` is the
//! Stefan-Boltzmann constant, and
//!
//! ```text
//!   L_absorbed = (1 − A) · L_sun · R² / (4 d²)
//! ```
//!
//! is the starlight the planet intercepts and absorbs at heliocentric distance
//! `d` (Bond albedo `A`). At hundreds of AU `L_absorbed ≪ L_int`, so `T_eff` is
//! nearly distance-independent and pinned at the ~40 K internal floor.
//!
//! The SED peak comes from Wien's displacement law, and we evaluate the Planck
//! flux density in a near-IR band (WISE W2, 4.6 µm) versus a far-IR / mm band
//! (350 µm) to show where the cold body is brightest.
//!
//! # Internal luminosity normalization (honest model input)
//!
//! Fortney et al.'s `L_int` comes from a full coupled atmosphere-interior
//! cooling calculation; we cannot re-derive their cooling tracks here. Instead
//! we **normalize** `L_int` so a nominal 10 M⊕, 4.5 Gyr Planet Nine lands at
//! the middle of the published 35–50 K band, then let the energy balance carry
//! the mass and distance dependence. The normalization constant
//! ([`L_INT_REF_W`] at [`REF_MASS_EARTH`]) is a labelled reference, documented
//! as a model input rather than a derived quantity. The published 35–50 K
//! band and the WISE/Spitzer detectability conclusions are kept as labelled
//! [`published`] constants and checked as residuals.

use std::f64::consts::PI;

use p9_core::analysis::photometry::{mass_radius_neptunian, ALBEDO_NEPTUNE};
use p9_core::analysis::thermal::{planck_bnu, thermal_flux_jy, C_LIGHT, L_SUN_W};
use p9_core::analysis::thermal::{R_EARTH_M, WIEN_B_M_K};
use p9_core::constants::AU_M;
use p9_core::units::{au, earth_masses, meters, Length, Mass};

pub mod published;

/// Stefan-Boltzmann constant (W / m² / K⁴), 2018 CODATA.
pub const SIGMA_SB: f64 = 5.670_374_419e-8;
/// WISE W2 effective wavelength (4.6 µm); the redder of the two short WISE
/// bands and the one the Meisner et al. (2018) shift-and-stack search uses.
pub const W2_WAVELENGTH_M: f64 = 4.6028e-6;
/// WISE W1 effective wavelength (3.4 µm).
pub const W1_WAVELENGTH_M: f64 = 3.3526e-6;
/// A representative far-IR / sub-mm band (350 µm; e.g. Herschel SPIRE /
/// ground-based sub-mm), well past the SED peak of a 40 K blackbody.
pub const FAR_IR_WAVELENGTH_M: f64 = 350.0e-6;

/// Reference mass for the internal-luminosity normalization (M⊕).
pub const REF_MASS_EARTH: f64 = 10.0;

/// Internal luminosity (W) of the reference [`REF_MASS_EARTH`] planet at
/// ~4.5 Gyr. This is the single labelled **model input**: it is chosen so the
/// energy balance places a nominal 10 M⊕ Planet Nine at ~42 K (mid-band of the
/// published 35–50 K), reproducing Fortney et al.'s thermal-evolution result.
/// For reference, a 3.34 R⊕ body radiating as a 42 K blackbody has
/// L = 4πR²σT⁴ ≈ 1.0×10¹⁵ W, the right order for a cooling few-Earth-mass ice
/// giant (Neptune itself has L_int ≈ 3×10¹⁵ W at 17 M⊕).
pub const L_INT_REF_W: f64 = 1.0e15;

/// Power-law index for how internal luminosity scales with mass at fixed age.
/// More massive bodies retain more formation/radiogenic heat and have larger
/// radiating areas; we take `L_int ∝ M` as a documented simple scaling (the
/// detailed exponent depends on the cooling tracks Fortney et al. integrate).
pub const L_INT_MASS_INDEX: f64 = 1.0;

/// Physical parameters of a candidate Planet Nine for the thermal model.
#[derive(Debug, Clone, Copy)]
pub struct P9Thermal {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Heliocentric distance in AU.
    pub distance_au: f64,
    /// Bond albedo for the absorbed-sunlight term.
    pub bond_albedo: f64,
    /// Internal luminosity (W) at the reference mass and ~4.5 Gyr.
    pub l_int_ref_w: f64,
}

impl P9Thermal {
    /// A Planet Nine of `mass_earth` at `distance_au`, with the Neptune-like
    /// Bond albedo and the reference internal luminosity.
    pub fn new(mass_earth: f64, distance_au: f64) -> Self {
        Self {
            mass_earth,
            distance_au,
            // Conflation kept deliberately: core's ALBEDO_NEPTUNE (0.41) is a
            // GEOMETRIC albedo; Neptune's Bond albedo is ~0.29. At 100s of AU
            // the equilibrium-temperature difference is well under a kelvin,
            // so the shared constant is used rather than a second knob.
            bond_albedo: ALBEDO_NEPTUNE,
            l_int_ref_w: L_INT_REF_W,
        }
    }

    /// Physical radius in meters (Neptune-anchored mass-radius relation,
    /// reused from `p9_core`).
    pub fn radius_m(&self) -> f64 {
        mass_radius_neptunian(self.mass_earth) * R_EARTH_M
    }

    /// Mass as a typed [`Mass`] (from the Earth-mass field).
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

    /// Internal luminosity (W) at this mass: the reference value scaled by
    /// `(M / M_ref)^L_INT_MASS_INDEX`.
    pub fn internal_luminosity_w(&self) -> f64 {
        self.l_int_ref_w * (self.mass_earth / REF_MASS_EARTH).powf(L_INT_MASS_INDEX)
    }

    /// Absorbed sunlight (W): the solar flux at distance `d` times the
    /// intercepting disc πR² times the absorbed fraction (1 − A):
    ///
    /// ```text
    ///   L_absorbed = (1 − A) · (L_sun / 4π d²) · π R²
    ///              = (1 − A) · L_sun · R² / (4 d²)
    /// ```
    pub fn absorbed_sunlight_w(&self) -> f64 {
        let d_m = self.distance_au * AU_M;
        let r = self.radius_m();
        (1.0 - self.bond_albedo) * L_SUN_W * r * r / (4.0 * d_m * d_m)
    }

    /// Effective temperature (K) from the energy balance
    /// `4πR²σ T⁴ = L_int + L_absorbed`.
    pub fn effective_temp(&self) -> f64 {
        let r = self.radius_m();
        let l_tot = self.internal_luminosity_w() + self.absorbed_sunlight_w();
        (l_tot / (4.0 * PI * r * r * SIGMA_SB)).powf(0.25)
    }

    /// Pure solar-equilibrium temperature (K) with **no** internal heat — the
    /// temperature a dead, airless body would reach. Used to show how feeble
    /// insolation is at these distances.
    pub fn solar_equilibrium_temp(&self) -> f64 {
        let r = self.radius_m();
        (self.absorbed_sunlight_w() / (4.0 * PI * r * r * SIGMA_SB)).powf(0.25)
    }

    /// Wien-law peak wavelength (m) of a blackbody at the effective
    /// temperature: λ_peak = b / T.
    pub fn sed_peak_wavelength_m(&self) -> f64 {
        WIEN_B_M_K / self.effective_temp()
    }

    /// SED peak wavelength as a typed [`Length`].
    pub fn sed_peak_wavelength(&self) -> Length {
        meters(self.sed_peak_wavelength_m())
    }

    /// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
    pub fn planck_bnu(t: f64, nu_hz: f64) -> f64 {
        planck_bnu(t, nu_hz)
    }

    /// Thermal (blackbody) flux density at `wavelength_m` in Jansky, observed
    /// from heliocentric distance `d` (opposition Δ ≈ d): a unit-emissivity
    /// emitter gives F_ν = π B_ν(T) (R/d)².
    pub fn flux_density_jy(&self, wavelength_m: f64) -> f64 {
        let nu = C_LIGHT / wavelength_m;
        thermal_flux_jy(self.effective_temp(), self.radius_m(), self.distance_au, nu)
    }

    /// Far-IR (350 µm) flux density in Jansky.
    pub fn far_ir_flux_jy(&self) -> f64 {
        self.flux_density_jy(FAR_IR_WAVELENGTH_M)
    }

    /// WISE W2 (4.6 µm) flux density in Jansky.
    pub fn w2_flux_jy(&self) -> f64 {
        self.flux_density_jy(W2_WAVELENGTH_M)
    }

    /// WISE W1 (3.4 µm) flux density in Jansky.
    pub fn w1_flux_jy(&self) -> f64 {
        self.flux_density_jy(W1_WAVELENGTH_M)
    }
}

/// Effective temperature (K) of a default Planet Nine of `mass_earth` at
/// `distance_au`.
pub fn effective_temp(mass_earth: f64, distance_au: f64) -> f64 {
    P9Thermal::new(mass_earth, distance_au).effective_temp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::analysis::thermal::{H_PLANCK, K_BOLTZ};

    #[test]
    fn nominal_p9_effective_temp_in_published_band() {
        // The headline: a 10 M⊕, 4.5 Gyr Planet Nine at the nominal ~500 AU
        // sits in the published 35–50 K band.
        let t = effective_temp(10.0, 500.0);
        assert!(
            t >= published::TEFF_MIN_K && t <= published::TEFF_MAX_K,
            "T_eff(10 M⊕, 500 AU) = {t:.1} K, outside {}–{} K",
            published::TEFF_MIN_K,
            published::TEFF_MAX_K
        );
    }

    #[test]
    fn five_to_twenty_earth_masses_span_the_band() {
        // Fortney et al.: ~35–50 K for 5–20 M⊕. With L_int ∝ M and the
        // mass-radius law, the whole 5–20 M⊕ range lands in (or just at the
        // edge of) the published band.
        for &m in &[5.0, 10.0, 15.0, 20.0] {
            let t = effective_temp(m, 500.0);
            assert!(
                t >= published::TEFF_MIN_K - 3.0 && t <= published::TEFF_MAX_K + 3.0,
                "T_eff({m} M⊕) = {t:.1} K outside the ~35–50 K band"
            );
        }
        // Heavier -> warmer (more retained heat, modest radius growth).
        assert!(effective_temp(20.0, 500.0) > effective_temp(5.0, 500.0));
    }

    #[test]
    fn insolation_is_negligible_beyond_250_au() {
        // The solar-equilibrium temperature (no internal heat) is far below
        // the internal-floor temperature at hundreds of AU.
        let p = P9Thermal::new(10.0, 500.0);
        assert!(
            p.solar_equilibrium_temp() < 12.0,
            "T_solar = {:.1} K should be tiny",
            p.solar_equilibrium_temp()
        );
        // Internal heat dominates the energy budget by orders of magnitude.
        assert!(
            p.internal_luminosity_w() > 100.0 * p.absorbed_sunlight_w(),
            "L_int {:.3e} should swamp L_absorbed {:.3e}",
            p.internal_luminosity_w(),
            p.absorbed_sunlight_w()
        );
    }

    #[test]
    fn teff_nearly_distance_independent_beyond_250_au() {
        // Because internal heat dominates, T_eff barely changes from 250 AU
        // out to 1000 AU.
        let near = effective_temp(10.0, 250.0);
        let far = effective_temp(10.0, 1000.0);
        assert!(
            (near - far).abs() < 0.5,
            "T_eff varies only {:.3} K over 250–1000 AU ({near:.2} → {far:.2})",
            (near - far).abs()
        );
        // And both sit in the published band.
        assert!(near >= published::TEFF_MIN_K && near <= published::TEFF_MAX_K);
        assert!(far >= published::TEFF_MIN_K && far <= published::TEFF_MAX_K);
    }

    #[test]
    fn sed_peak_is_in_the_far_infrared() {
        // A ~42 K blackbody peaks near λ = b/T ≈ 69 µm — squarely in the
        // far-IR, tens of microns, nowhere near the WISE near-IR bands.
        let p = P9Thermal::new(10.0, 500.0);
        let peak_um = p.sed_peak_wavelength_m() * 1e6;
        assert!(
            (40.0..=120.0).contains(&peak_um),
            "SED peak {peak_um:.1} µm should be in the far-IR (tens of µm)"
        );
        // Far redward of WISE W1/W2 (3.4 / 4.6 µm).
        assert!(p.sed_peak_wavelength_m() > 10.0 * W2_WAVELENGTH_M);
    }

    #[test]
    fn far_ir_flux_exceeds_near_ir() {
        // The cold body is vastly brighter in the far-IR/mm than at WISE
        // wavelengths on the Wien tail.
        let p = P9Thermal::new(10.0, 500.0);
        let far = p.far_ir_flux_jy();
        let w2 = p.w2_flux_jy();
        let w1 = p.w1_flux_jy();
        assert!(far > w2, "far-IR {far:.3e} Jy should exceed W2 {w2:.3e} Jy");
        assert!(w2 > w1, "W2 should exceed W1 on the Wien tail");
        // The W1 flux is utterly negligible (deep on the Wien tail).
        assert!(w1 < 1e-20, "W1 flux {w1:.3e} Jy should be vanishing");
        // The contrast spans dozens of orders of magnitude.
        assert!(
            far / w2 > 1e20,
            "far-IR/W2 contrast = {:.3e} should be enormous",
            far / w2
        );
    }

    #[test]
    fn energy_balance_recovers_stefan_boltzmann() {
        // Self-consistency: 4πR²σ T_eff⁴ must equal L_int + L_absorbed.
        let p = P9Thermal::new(12.0, 600.0);
        let r = p.radius_m();
        let l_rad = 4.0 * PI * r * r * SIGMA_SB * p.effective_temp().powi(4);
        let l_tot = p.internal_luminosity_w() + p.absorbed_sunlight_w();
        assert_relative_eq!(l_rad, l_tot, max_relative = 1e-12);
    }

    #[test]
    fn typed_accessors_match_f64() {
        let p = P9Thermal::new(10.0, 500.0);
        assert_relative_eq!(
            (p.mass() / earth_masses(1.0)).value,
            p.mass_earth,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (p.distance() / au(1.0)).value,
            p.distance_au,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (p.radius() / meters(1.0)).value,
            p.radius_m(),
            max_relative = 1e-12
        );
        assert_relative_eq!(
            (p.sed_peak_wavelength() / meters(1.0)).value,
            p.sed_peak_wavelength_m(),
            max_relative = 1e-12
        );
    }

    #[test]
    fn radius_matches_ice_giant_models() {
        // Sanity: the reused mass-radius law gives ~3.3 R⊕ at 10 M⊕.
        let r = P9Thermal::new(10.0, 500.0).radius_m() / R_EARTH_M;
        assert!((3.0..3.6).contains(&r), "R(10 M⊕) = {r:.2} R⊕");
    }

    #[test]
    fn wien_peak_matches_blackbody_planck_maximum() {
        // Cross-check Wien's law against a direct scan of B_λ(T): the peak of
        // the Planck spectral radiance in wavelength must coincide with b/T.
        let t = 42.0;
        let lam_wien = WIEN_B_M_K / t;
        // B_λ ∝ ν³ B_ν · (c/λ²)?  Easier: B_λ(λ) = 2hc²/λ⁵ / (exp(hc/λkT)−1).
        let b_lambda = |lam: f64| {
            let x = H_PLANCK * C_LIGHT / (lam * K_BOLTZ * t);
            (2.0 * H_PLANCK * C_LIGHT * C_LIGHT / lam.powi(5)) / (x.exp() - 1.0)
        };
        let mut best_lam = 0.0;
        let mut best_val = 0.0;
        for k in 1..=20_000 {
            let lam = k as f64 * 1e-7; // 0.1 µm steps up to 2 mm
            let v = b_lambda(lam);
            if v > best_val {
                best_val = v;
                best_lam = lam;
            }
        }
        assert_relative_eq!(best_lam, lam_wien, max_relative = 1e-3);
    }

    #[test]
    fn published_band_residual_is_documented() {
        // Honest residual: our normalized model places the nominal 10 M⊕ body
        // at the mid-band; the spread across 5–20 M⊕ is a model output, not a
        // direct reproduction of Fortney et al.'s cooling tracks.
        let t10 = effective_temp(10.0, 500.0);
        let mid = 0.5 * (published::TEFF_MIN_K + published::TEFF_MAX_K);
        assert!(
            (t10 - mid).abs() < 5.0,
            "nominal T_eff {t10:.1} K within 5 K of band midpoint {mid:.1} K"
        );
    }
}
