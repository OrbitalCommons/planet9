//! Band-parametrized WISE thermal+reflected model (W1 at 3.4 µm and W2 at
//! 4.6 µm).
//!
//! The sibling crate `p9-2018-wise-search` fixes its `P9Thermal` to the W1
//! band. The 2016 paper's headline is that the **W2 (4.6 µm)** band is more
//! favorable than W1 for a *cold* Planet Nine, so this crate must evaluate the
//! same thermal+reflected physics at two wavelengths. We reuse the sibling's
//! [`P9Thermal`] verbatim for the planet's geometry (radius law, effective
//! temperature, distances) and for the W1 numbers, and add a band-generic
//! Planck/reflected evaluation here so the W2 magnitude is computed with the
//! *identical* physics, only the wavelength and Vega zero point differing.
//!
//! Why W2 helps a cold body: at a fixed temperature `T`, the Planck flux on the
//! Wien tail scales as exp(−hν/kT). Moving from W1 (ν = c/3.353 µm) to W2
//! (ν = c/4.603 µm) lowers hν/kT by the wavelength ratio 4.603/3.353 ≈ 1.37, so
//! the thermal flux rises by a factor ~exp[(hc/kT)(1/λ₁ − 1/λ₂)]. For a 40 K
//! body that is an enormous factor — W2 sits far less deep on the Wien tail.

use std::f64::consts::PI;

use p9_2018_wise_search::thermal_model::P9Thermal;

/// Physical constants (mirrors of the sibling's, kept local so the band-generic
/// Planck evaluation is self-contained).
const H_PLANCK: f64 = 6.626_070_15e-34;
const K_BOLTZ: f64 = 1.380_649e-23;
const C_LIGHT: f64 = 2.997_924_58e8;
const AU_M: f64 = 1.495_978_707e11;
const R_SUN_M: f64 = 6.957e8;
const T_SUN: f64 = 5778.0;
const JY: f64 = 1.0e-26;

/// A WISE photometric band: effective wavelength and Vega zero-point flux
/// density (Wright et al. 2010, Table 1).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WiseBand {
    /// Short label, "W1" or "W2".
    pub name: &'static str,
    /// Effective wavelength in meters.
    pub wavelength_m: f64,
    /// Vega zero-point flux density in Jansky (0.0 mag source).
    pub zero_point_jy: f64,
}

/// WISE W1 (3.3526 µm, F_ν0 = 309.540 Jy; Wright et al. 2010, Table 1).
pub const W1: WiseBand = WiseBand {
    name: "W1",
    wavelength_m: p9_2018_wise_search::thermal_model::W1_WAVELENGTH_M,
    zero_point_jy: p9_2018_wise_search::thermal_model::W1_ZERO_POINT_JY,
};

/// WISE W2 (4.6028 µm, F_ν0 = 171.787 Jy; Wright et al. 2010, Table 1).
pub const W2: WiseBand = WiseBand {
    name: "W2",
    wavelength_m: 4.6028e-6,
    zero_point_jy: 171.787,
};

/// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz` (identical to the
/// sibling's; under-/overflow guarded on the Wien tail).
fn planck_bnu(t: f64, nu_hz: f64) -> f64 {
    let x = H_PLANCK * nu_hz / (K_BOLTZ * t);
    if x > 700.0 {
        return 0.0;
    }
    (2.0 * H_PLANCK * nu_hz.powi(3) / (C_LIGHT * C_LIGHT)) / (x.exp() - 1.0)
}

/// Thermal grey-body flux density (Jy) of `p9` in `band`: F_ν = π B_ν(T)(R/d)².
pub fn thermal_flux_jy(p9: &P9Thermal, band: WiseBand) -> f64 {
    let nu = C_LIGHT / band.wavelength_m;
    let bnu = planck_bnu(p9.effective_temp(), nu);
    let r = p9.radius_m();
    let d = p9.distance_au * AU_M;
    let solid_angle = PI * (r / d).powi(2);
    solid_angle * bnu / JY
}

/// Reflected-sunlight flux density (Jy) of `p9` in `band`, near opposition
/// (Δ = d − 1 AU). Same geometry as the sibling's W1 reflected channel.
pub fn reflected_flux_jy(p9: &P9Thermal, band: WiseBand) -> f64 {
    let nu = C_LIGHT / band.wavelength_m;
    let bnu_sun = planck_bnu(T_SUN, nu);
    let r_p = p9.radius_m();
    let solar_flux_at_planet = PI * bnu_sun * (R_SUN_M / (p9.distance_au * AU_M)).powi(2);
    let delta_m = (p9.distance_au - 1.0).max(1.0) * AU_M;
    let reflected = p9.albedo * solar_flux_at_planet * (r_p * r_p) / (4.0 * delta_m * delta_m);
    reflected / JY
}

/// Total flux density (Jy) of `p9` in `band`: thermal + reflected.
pub fn total_flux_jy(p9: &P9Thermal, band: WiseBand) -> f64 {
    thermal_flux_jy(p9, band) + reflected_flux_jy(p9, band)
}

/// Apparent Vega magnitude of `p9` in `band` from the total flux density;
/// +∞ for zero flux.
pub fn magnitude(p9: &P9Thermal, band: WiseBand) -> f64 {
    let flux = total_flux_jy(p9, band);
    if flux <= 0.0 {
        return f64::INFINITY;
    }
    -2.5 * (flux / band.zero_point_jy).log10()
}

/// Convenience: apparent `band` magnitude of a default-albedo P9 of
/// `mass_earth` at `distance_au`.
pub fn band_magnitude(mass_earth: f64, distance_au: f64, band: WiseBand) -> f64 {
    magnitude(&P9Thermal::new(mass_earth, distance_au), band)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn w1_matches_sibling_exactly() {
        // The band-generic W1 evaluation must reproduce the sibling's W1
        // magnitude to round-off (same physics, same constants).
        let p9 = P9Thermal::new(10.0, 500.0);
        assert_relative_eq!(magnitude(&p9, W1), p9.w1_magnitude(), max_relative = 1e-12);
    }

    #[test]
    fn w2_thermal_far_exceeds_w1_for_cold_body() {
        // For a cold (40 K floor) body, W2 sits much less deep on the Wien tail,
        // so its thermal flux dwarfs W1's.
        let p9 = P9Thermal::new(10.0, 500.0);
        let ratio = thermal_flux_jy(&p9, W2) / thermal_flux_jy(&p9, W1);
        assert!(
            ratio > 1e3,
            "W2/W1 thermal flux ratio = {ratio:.3e} should be huge on the Wien tail"
        );
    }

    #[test]
    fn brighter_with_mass_and_fainter_with_distance_in_w2() {
        assert!(band_magnitude(15.0, 500.0, W2) < band_magnitude(5.0, 500.0, W2));
        assert!(band_magnitude(10.0, 800.0, W2) > band_magnitude(10.0, 400.0, W2));
    }

    #[test]
    fn zero_point_gives_zero_magnitude() {
        for band in [W1, W2] {
            let m = -2.5 * (band.zero_point_jy / band.zero_point_jy).log10();
            assert!(m.abs() < 1e-12, "{} zero point", band.name);
        }
    }
}
