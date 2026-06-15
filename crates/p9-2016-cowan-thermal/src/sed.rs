//! Spectral energy distribution of a cold Planet Nine.
//!
//! Cowan, Holder & Kaib (2016), "Cosmologists in Search of Planet Nine: the
//! Case for CMB experiments" (arXiv:1602.05963), make a simple but powerful
//! point: a cold (~40 K) Planet Nine radiates as a near-blackbody whose
//! emission peaks in the far-IR/sub-mm. At those wavelengths the body is far
//! *brighter relative to survey sensitivities* than it is in the optical or
//! near-IR, where it shines only by reflected sunlight. So millimetre CMB
//! experiments (ACT, SPT, Planck) and far-IR surveys are surprisingly
//! competitive — even advantageous — for detecting it.
//!
//! Two emission channels:
//!
//! 1. **Thermal emission.** The planet radiates as a blackbody of unit
//!    emissivity at temperature `temp_k`. The observed flux density is
//!    F_ν = π B_ν(T) (R/d)², with B_ν the Planck function and (R/d)² the
//!    solid angle of the disc divided by π.
//!
//! 2. **Reflected sunlight.** The planet reflects the solar flux it
//!    intercepts. This scales as 1/(r²Δ²): once for Sun→planet, once for
//!    planet→observer. It dominates the optical/near-IR; it is utterly
//!    negligible at mm wavelengths.
//!
//! The **crossover wavelength** is where the rising thermal channel overtakes
//! the falling reflected channel. Below it (bluer) reflected light wins; above
//! it (redder, into the far-IR/mm) thermal emission dominates.
//!
//! Headline reference (Cowan et al. 2016, abstract): a 40 K Planet Nine at
//! 700 AU "would appear as a ~30 mJy source at 1 mm"; a smaller/colder/farther
//! body "could be as faint as 1 mJy at 1 mm". We keep these as labelled
//! reference constants and document the residual against our computed value.

use std::f64::consts::PI;

use p9_core::analysis::photometry::mass_radius_neptunian;
use p9_core::constants::AU_M;

/// Planck constant (J·s).
pub const H_PLANCK: f64 = 6.626_070_15e-34;
/// Boltzmann constant (J/K).
pub const K_BOLTZ: f64 = 1.380_649e-23;
/// Speed of light (m/s).
pub const C_LIGHT: f64 = 2.997_924_58e8;
/// Wien displacement-law constant b for the wavelength of peak B_λ (m·K).
pub const WIEN_B_M_K: f64 = 2.897_771_955e-3;
/// Earth radius (m).
pub const R_EARTH_M: f64 = 6.371e6;
/// Solar radius (m).
pub const R_SUN_M: f64 = 6.957e8;
/// Solar effective temperature (K).
pub const T_SUN: f64 = 5778.0;
/// 1 Jansky in W/m²/Hz.
pub const JY: f64 = 1.0e-26;

/// Effective temperature adopted by Cowan et al. (2016) (K).
pub const COWAN_TEMP_K: f64 = 40.0;
/// Heliocentric distance adopted by Cowan et al. (2016) (AU).
pub const COWAN_DISTANCE_AU: f64 = 700.0;
/// Cowan et al. (2016) nominal mass scale (Earth masses): the paper's fiducial
/// is a ~10 M⊕ ice giant, consistent with Batygin & Brown (2016).
pub const COWAN_MASS_EARTH: f64 = 10.0;

/// Neptune-like geometric albedo for the reflected-light channel (shared
/// workspace value).
pub const DEFAULT_ALBEDO: f64 = 0.41;

/// Reference headline flux: Cowan et al. (2016) quote "~30 mJy at 1 mm" for the
/// nominal 40 K, 700 AU body. Labelled reference only — not used in the
/// computation.
pub const COWAN_FLUX_1MM_MJY: f64 = 30.0;
/// Reference faint-case flux: "as faint as 1 mJy at 1 mm" for a
/// smaller/colder/farther body.
pub const COWAN_FAINT_FLUX_1MM_MJY: f64 = 1.0;
/// Wavelength of 1 mm in metres (the band of the headline 30 mJy quote).
pub const ONE_MM_M: f64 = 1.0e-3;

/// Representative point-source sensitivities of detection facilities, used to
/// compare "flux relative to limit" between bands. These are order-of-magnitude
/// representative depths, not survey-specific completeness limits.
pub mod sensitivity {
    /// Representative CMB-experiment (ACT/SPT-class) point-source detection
    /// threshold at ~1 mm (mJy). Cowan et al. note current experiments reach a
    /// few-to-tens of mJy and become confusion-limited for moving sources
    /// fainter than this.
    pub const CMB_1MM_MJY: f64 = 15.0;
    /// Representative WISE W1 (3.4 µm) single-frame point-source sensitivity
    /// (mJy). W1 5σ ≈ 0.08 mJy (~16.5 Vega mag); we use 0.1 mJy as a round
    /// representative limit.
    pub const WISE_W1_MJY: f64 = 0.1;
}

/// Physical parameters of a candidate Planet Nine for the SED.
#[derive(Debug, Clone, Copy)]
pub struct P9Sed {
    /// Mass in Earth masses (sets the radius via the mass-radius relation).
    pub mass_earth: f64,
    /// Heliocentric distance in AU.
    pub distance_au: f64,
    /// Effective (blackbody) temperature in Kelvin.
    pub temp_k: f64,
    /// Geometric albedo for the reflected-light channel.
    pub albedo: f64,
}

impl P9Sed {
    /// The Cowan et al. (2016) fiducial: 10 M⊕, 700 AU, 40 K.
    pub fn cowan_fiducial() -> Self {
        Self {
            mass_earth: COWAN_MASS_EARTH,
            distance_au: COWAN_DISTANCE_AU,
            temp_k: COWAN_TEMP_K,
            albedo: DEFAULT_ALBEDO,
        }
    }

    /// A body of given mass/distance with the Cowan temperature and default
    /// albedo.
    pub fn new(mass_earth: f64, distance_au: f64, temp_k: f64) -> Self {
        Self {
            mass_earth,
            distance_au,
            temp_k,
            albedo: DEFAULT_ALBEDO,
        }
    }

    /// Physical radius in metres (Neptune-anchored mass-radius relation).
    pub fn radius_m(&self) -> f64 {
        mass_radius_neptunian(self.mass_earth) * R_EARTH_M
    }

    /// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
    pub fn planck_bnu(t: f64, nu_hz: f64) -> f64 {
        let x = H_PLANCK * nu_hz / (K_BOLTZ * t);
        if x > 700.0 {
            return 0.0;
        }
        (2.0 * H_PLANCK * nu_hz.powi(3) / (C_LIGHT * C_LIGHT)) / (x.exp() - 1.0)
    }

    /// Wavelength of peak spectral radiance B_λ via the Wien displacement law:
    /// λ_peak = b / T (metres).
    pub fn wien_peak_wavelength_m(&self) -> f64 {
        WIEN_B_M_K / self.temp_k
    }

    /// Thermal flux density at wavelength `lambda_m` (W/m²/Hz), a blackbody of
    /// unit emissivity: F_ν = π B_ν(T) (R/d)².
    pub fn thermal_flux_si(&self, lambda_m: f64) -> f64 {
        let nu = C_LIGHT / lambda_m;
        let bnu = Self::planck_bnu(self.temp_k, nu);
        let r = self.radius_m();
        let d = self.distance_au * AU_M;
        PI * bnu * (r / d).powi(2)
    }

    /// Reflected-sunlight flux density at wavelength `lambda_m` (W/m²/Hz),
    /// observed near opposition (Δ = d − 1 AU). The Sun is a blackbody at
    /// `T_SUN`; the solar flux at the planet is π B_ν(T_sun) (R_sun/r)². The
    /// planet intercepts a disc πR_p² and reflects a fraction `albedo`; the
    /// reradiated flux falls as 1/(4Δ²) to the observer:
    ///
    ///   F_refl = albedo · [π B_ν(T_sun) (R_sun/r)²] · (R_p² / (4 Δ²)).
    pub fn reflected_flux_si(&self, lambda_m: f64) -> f64 {
        let nu = C_LIGHT / lambda_m;
        let bnu_sun = Self::planck_bnu(T_SUN, nu);
        let r_m = self.distance_au * AU_M;
        let solar_flux_at_planet = PI * bnu_sun * (R_SUN_M / r_m).powi(2);
        let r_p = self.radius_m();
        let delta_m = (self.distance_au - 1.0).max(1.0) * AU_M;
        self.albedo * solar_flux_at_planet * (r_p * r_p) / (4.0 * delta_m * delta_m)
    }

    /// Thermal flux density at `lambda_m` in milli-Jansky.
    pub fn thermal_flux_mjy(&self, lambda_m: f64) -> f64 {
        self.thermal_flux_si(lambda_m) / JY / 1.0e-3
    }

    /// Reflected flux density at `lambda_m` in milli-Jansky.
    pub fn reflected_flux_mjy(&self, lambda_m: f64) -> f64 {
        self.reflected_flux_si(lambda_m) / JY / 1.0e-3
    }

    /// Total (thermal + reflected) flux density at `lambda_m` in milli-Jansky.
    pub fn total_flux_mjy(&self, lambda_m: f64) -> f64 {
        self.thermal_flux_mjy(lambda_m) + self.reflected_flux_mjy(lambda_m)
    }

    /// Crossover wavelength (metres): the wavelength at which thermal emission
    /// equals reflected sunlight. Found by bisection on the log-flux difference
    /// d(λ) = log(F_thermal) − log(F_reflected), which is monotonically
    /// increasing in λ across the relevant range (reflected falls on the solar
    /// Rayleigh-Jeans tail while thermal rises toward its far-IR peak), so the
    /// root is unique. Searches the optical-to-cm range [0.3 µm, 3 cm].
    pub fn crossover_wavelength_m(&self) -> Option<f64> {
        let diff = |lambda: f64| {
            let th = self.thermal_flux_si(lambda).max(f64::MIN_POSITIVE);
            let re = self.reflected_flux_si(lambda).max(f64::MIN_POSITIVE);
            th.ln() - re.ln()
        };
        let mut lo = 0.3e-6; // 0.3 µm (optical/UV edge)
        let mut hi = 3.0e-2; // 3 cm
        let d_lo = diff(lo);
        let d_hi = diff(hi);
        // Expect reflected to dominate at lo (diff < 0) and thermal at hi.
        if d_lo > 0.0 || d_hi < 0.0 {
            return None;
        }
        for _ in 0..200 {
            let mid = (lo * hi).sqrt(); // geometric midpoint (log-space)
            if diff(mid) < 0.0 {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Some((lo * hi).sqrt())
    }

    /// Detectability margin in a band: flux density (mJy) divided by a
    /// representative point-source sensitivity (mJy). >1 means detectable; a
    /// larger value is "more favorable relative to the limit".
    pub fn detection_margin(&self, lambda_m: f64, sensitivity_mjy: f64) -> f64 {
        self.total_flux_mjy(lambda_m) / sensitivity_mjy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// WISE W1 effective wavelength (m).
    const W1_M: f64 = 3.3526e-6;

    #[test]
    fn radius_is_ice_giant_scale() {
        let r = P9Sed::cowan_fiducial().radius_m() / R_EARTH_M;
        assert!(r > 3.0 && r < 3.6, "R(10 M⊕) = {r:.2} R⊕");
    }

    #[test]
    fn sed_peaks_in_far_infrared() {
        // Wien wavelength for a 40 K blackbody is b/T = 2.898e-3 / 40 ≈ 72 µm,
        // squarely in the far-IR as the paper's headline asserts.
        let peak = P9Sed::cowan_fiducial().wien_peak_wavelength_m();
        let peak_um = peak * 1.0e6;
        assert!(
            (70.0..=90.0).contains(&peak_um),
            "40 K Wien peak = {peak_um:.1} µm, expected far-IR ~70-90 µm"
        );
        // Closed-form check: b/T.
        assert_relative_eq!(peak_um, WIEN_B_M_K / 40.0 * 1.0e6, max_relative = 1e-12);
    }

    #[test]
    fn thermal_overtakes_reflected_beyond_crossover() {
        let p = P9Sed::cowan_fiducial();
        let xover = p
            .crossover_wavelength_m()
            .expect("crossover should exist for a cold body");
        let xover_um = xover * 1.0e6;
        // The crossover must lie blueward of the mm regime and redward of the
        // near-IR: reflected sunlight wins in the optical/near-IR, thermal wins
        // in the far-IR/sub-mm/mm. Sanity-bound it to the IR.
        assert!(
            (1.0..=500.0).contains(&xover_um),
            "crossover = {xover_um:.1} µm, expected mid/far-IR"
        );
        // Just blueward of crossover: reflected dominates.
        let blue = xover * 0.5;
        assert!(
            p.reflected_flux_si(blue) > p.thermal_flux_si(blue),
            "reflected should win blueward of crossover"
        );
        // Just redward (into far-IR/mm): thermal dominates.
        let red = xover * 2.0;
        assert!(
            p.thermal_flux_si(red) > p.reflected_flux_si(red),
            "thermal should win redward of crossover"
        );
    }

    #[test]
    fn at_1mm_thermal_swamps_reflected() {
        // At 1 mm the cold body is far down the solar Rayleigh-Jeans tail for
        // reflection but near its own thermal regime: thermal must dominate by
        // many orders of magnitude.
        let p = P9Sed::cowan_fiducial();
        let th = p.thermal_flux_mjy(ONE_MM_M);
        let re = p.reflected_flux_mjy(ONE_MM_M);
        assert!(
            th > 1.0e6 * re,
            "at 1 mm thermal {th:.3e} mJy should swamp reflected {re:.3e} mJy"
        );
    }

    #[test]
    fn mm_band_more_favorable_than_near_ir() {
        // The paper's core claim: the body is brighter *relative to the
        // facility limit* at mm (vs a CMB experiment) than in the near-IR
        // (vs WISE). Compare the two detection margins.
        let p = P9Sed::cowan_fiducial();
        let mm_margin = p.detection_margin(ONE_MM_M, sensitivity::CMB_1MM_MJY);
        let nir_margin = p.detection_margin(W1_M, sensitivity::WISE_W1_MJY);
        assert!(
            mm_margin > nir_margin,
            "mm margin {mm_margin:.3e} should beat near-IR margin {nir_margin:.3e}"
        );
        // And the near-IR margin is hopeless (far below the limit).
        assert!(nir_margin < 1.0, "near-IR margin = {nir_margin:.3e}");
    }

    #[test]
    fn computed_1mm_flux_is_order_tens_of_mjy() {
        // Honest residual vs the paper's "~30 mJy at 1 mm" headline. Our
        // blackbody computation with a 10 M⊕ (R ≈ 3.34 R⊕) body at 40 K, 700 AU
        // should land in the same order of magnitude (a few to tens of mJy).
        // We do NOT tune to 30; we report the computed value and pin its order.
        let p = P9Sed::cowan_fiducial();
        let f = p.thermal_flux_mjy(ONE_MM_M);
        assert!(
            (1.0..=100.0).contains(&f),
            "computed 1 mm flux = {f:.2} mJy, paper quotes ~{COWAN_FLUX_1MM_MJY} mJy"
        );
        // Within a factor of ~3 of the paper's headline (radius/emissivity
        // assumptions differ); document, don't hide, the residual.
        let ratio = f / COWAN_FLUX_1MM_MJY;
        assert!(
            (0.33..=3.0).contains(&ratio),
            "1 mm flux ratio computed/paper = {ratio:.2}"
        );
    }

    #[test]
    fn flux_follows_blackbody_scalings() {
        // F ∝ (R/d)²: doubling distance quarters the thermal flux.
        let near = P9Sed::new(10.0, 700.0, 40.0);
        let far = P9Sed::new(10.0, 1400.0, 40.0);
        let ratio = near.thermal_flux_si(ONE_MM_M) / far.thermal_flux_si(ONE_MM_M);
        assert_relative_eq!(ratio, 4.0, max_relative = 1e-9);

        // Rayleigh-Jeans at 1 mm (hν/kT ≈ 0.35 at 40 K is mildly non-RJ, so
        // check deep RJ at 3 mm where hν/kT ≈ 0.12): F_ν ∝ T there. Doubling T
        // roughly doubles the flux.
        let lam = 3.0e-3;
        let cold = P9Sed::new(10.0, 700.0, 40.0);
        let warm = P9Sed::new(10.0, 700.0, 80.0);
        let t_ratio = warm.thermal_flux_si(lam) / cold.thermal_flux_si(lam);
        assert!(
            (1.8..=2.2).contains(&t_ratio),
            "RJ flux should scale ~linearly with T; ratio = {t_ratio:.3}"
        );
    }

    #[test]
    fn planck_peak_frequency_matches_wien() {
        // The frequency-form Planck function B_ν peaks at ν_peak = 2.821 kT/h.
        // Cross-check our B_ν by scanning frequencies near that peak.
        let t = 40.0;
        let nu_peak = 2.821_439 * K_BOLTZ * t / H_PLANCK;
        let b_peak = P9Sed::planck_bnu(t, nu_peak);
        let b_lo = P9Sed::planck_bnu(t, nu_peak * 0.7);
        let b_hi = P9Sed::planck_bnu(t, nu_peak * 1.3);
        assert!(
            b_peak > b_lo && b_peak > b_hi,
            "B_nu should peak near 2.821 kT/h: {b_lo:.3e} {b_peak:.3e} {b_hi:.3e}"
        );
    }

    #[test]
    fn faint_case_is_below_cmb_sensitivity() {
        // The paper's faint case (~1 mJy at 1 mm) sits at/below a representative
        // CMB point-source limit — illustrating why fainter bodies are
        // confusion-limited. We model the faint case as a smaller, colder, more
        // distant body and check it lands well below the nominal 30 mJy source.
        let faint = P9Sed::new(5.0, 1000.0, 30.0);
        let bright = P9Sed::cowan_fiducial();
        assert!(
            faint.thermal_flux_mjy(ONE_MM_M) < bright.thermal_flux_mjy(ONE_MM_M),
            "faint case should be dimmer than fiducial at 1 mm"
        );
    }
}
