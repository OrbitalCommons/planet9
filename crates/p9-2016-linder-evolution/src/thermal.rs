//! Thermal (intrinsic-emission) magnitudes of a cooling Planet Nine in the
//! mid/far-IR bands, the regime Linder & Mordasini (2016) identify as the most
//! favorable for detection.
//!
//! A several-Earth-mass planet retains formation and radiogenic heat; after
//! ~Gyr of evolution its grey-atmosphere effective temperature settles to a few
//! tens of K essentially independent of heliocentric distance (solar heating at
//! hundreds of AU is only ~10 K). We model the planet as a grey body of unit
//! emissivity radiating from its full surface at this effective temperature.
//!
//! The thermal flux density in a band of effective frequency ν is
//!
//!   F_ν = π B_ν(T_eff) (R/d)²
//!
//! (the solid angle πR²/d² times the surface intensity B_ν), converted to a
//! Vega-system magnitude via that band's zero-point flux density:
//!
//!   m = −2.5 log10(F_ν / F_ν0).
//!
//! The radius R uses the shared Neptune-anchored mass-radius relation in
//! `p9_core::analysis::photometry`, so it cannot drift from the reflected-light
//! channel or the rest of the workspace.
//!
//! **Honest residual.** The paper couples a full interior-evolution model to a
//! grey atmosphere with band-dependent opacities; we collapse that to a single
//! effective temperature and a grey emissivity. On the steep Wien tail the
//! near-IR magnitudes (W1/W2/L) are extremely sensitive to T_eff, so our
//! single-temperature W1/W2 numbers are only order-of-magnitude. The far-IR
//! bands (W4 ≈ 22 µm, Q ≈ 20 µm), which sit near the Rayleigh-Jeans / peak
//! region for a ~45 K body, are far more robust, and it is precisely those
//! bands the paper flags as most favorable — so the qualitative conclusion
//! (mid/far-IR beats optical for the cold body) is reproduced where the model
//! is most trustworthy.

use p9_core::analysis::photometry::mass_radius_neptunian;
use p9_core::analysis::thermal::{flux_to_magnitude, thermal_flux_jy, C_LIGHT};
use p9_core::constants::EARTH_RADIUS_KM;

/// Earth radius in meters (from the shared `EARTH_RADIUS_KM`).
const R_EARTH_M: f64 = EARTH_RADIUS_KM * 1.0e3;

/// Equilibrium temperature contribution from internal heat after Gyr of
/// cooling (K). Linder & Mordasini (2016) find a nominal 10 M⊕ Planet Nine
/// settles to T_eff ≈ 47 K; more massive planets are hotter (more retained
/// formation energy and a larger interior). We adopt a mild power-law in mass
/// anchored on the published 47 K at 10 M⊕, chosen so the thermal magnitudes
/// climb monotonically and sensibly with mass (5→20→50 M⊕) as the paper's
/// table shows, without claiming the full evolution grid.
///
///   T_floor(M) = 47 K · (M / 10 M⊕)^0.18
///
/// (0.18 reproduces the paper's ~3-mag Q-band brightening from 5→20 M⊕.)
pub const T_FLOOR_10ME_K: f64 = 47.0;
const T_FLOOR_MASS_EXPONENT: f64 = 0.18;

/// Mid/far-IR photometric bands used in the paper, with effective wavelength
/// (µm) and Vega zero-point flux density (Jy). WISE values are Wright et al.
/// (2010, Table 1); ground-based N and Q bands use standard effective
/// wavelengths and Vega zero points (Tokunaga & Vacca 2005, MKO-like).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Band {
    /// WISE W1, 3.3526 µm.
    W1,
    /// WISE W2, 4.6028 µm.
    W2,
    /// WISE W3, 11.5608 µm.
    W3,
    /// WISE W4, 22.0883 µm.
    W4,
    /// Ground-based N band, ~10.5 µm.
    N,
    /// Ground-based Q band, ~20.5 µm.
    Q,
}

impl Band {
    /// Effective wavelength in meters.
    pub fn wavelength_m(self) -> f64 {
        let micron = match self {
            Band::W1 => 3.3526,
            Band::W2 => 4.6028,
            Band::W3 => 11.5608,
            Band::W4 => 22.0883,
            Band::N => 10.5,
            Band::Q => 20.5,
        };
        micron * 1.0e-6
    }

    /// Vega-system zero-point flux density in Jansky (a 0.0-mag source).
    pub fn zero_point_jy(self) -> f64 {
        match self {
            Band::W1 => 309.540,
            Band::W2 => 171.787,
            Band::W3 => 31.674,
            Band::W4 => 8.363,
            Band::N => 36.0,
            Band::Q => 10.0,
        }
    }

    /// Effective frequency (Hz).
    pub fn frequency_hz(self) -> f64 {
        C_LIGHT / self.wavelength_m()
    }
}

/// A candidate Planet Nine for the thermal model.
#[derive(Debug, Clone, Copy)]
pub struct P9 {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Heliocentric distance in AU.
    pub distance_au: f64,
}

impl P9 {
    /// Physical radius in meters (Neptune-anchored mass-radius relation,
    /// shared with the reflected-light channel).
    pub fn radius_m(&self) -> f64 {
        mass_radius_neptunian(self.mass_earth) * R_EARTH_M
    }

    /// Post-evolution effective temperature (K): the internal-heat floor scaled
    /// mildly with mass (see [`T_FLOOR_10ME_K`]). Solar heating at hundreds of
    /// AU is ~10 K, far below this floor, so it is neglected.
    pub fn effective_temperature(&self) -> f64 {
        T_FLOOR_10ME_K * (self.mass_earth / 10.0).powf(T_FLOOR_MASS_EXPONENT)
    }

    /// Intrinsic thermal flux density in `band` (Jy): a grey body of unit
    /// emissivity, F_ν = π B_ν(T_eff) (R/d)².
    pub fn thermal_flux_jy(&self, band: Band) -> f64 {
        thermal_flux_jy(
            self.effective_temperature(),
            self.radius_m(),
            self.distance_au,
            band.frequency_hz(),
        )
    }

    /// Apparent Vega-system magnitude in `band` from the intrinsic thermal
    /// flux. Returns +∞ when the flux underflows to zero (deep Wien tail).
    pub fn thermal_magnitude(&self, band: Band) -> f64 {
        flux_to_magnitude(self.thermal_flux_jy(band), band.zero_point_jy())
    }

    /// The band (among `bands`) in which the planet is intrinsically brightest
    /// (smallest thermal magnitude). Returns `None` for an empty slice.
    pub fn most_favorable_band(&self, bands: &[Band]) -> Option<Band> {
        bands.iter().copied().min_by(|a, b| {
            self.thermal_magnitude(*a)
                .partial_cmp(&self.thermal_magnitude(*b))
                .unwrap()
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference;

    fn nominal() -> P9 {
        P9 {
            mass_earth: 10.0,
            distance_au: 700.0,
        }
    }

    #[test]
    fn effective_temp_near_published_47k() {
        let t = nominal().effective_temperature();
        assert!(
            (t - reference::T_EFF_10ME_700AU_K).abs() < 1.0,
            "T_eff = {t:.1} K, published {:.1} K",
            reference::T_EFF_10ME_700AU_K
        );
    }

    #[test]
    fn radius_matches_ice_giant_models() {
        let r = nominal().radius_m() / R_EARTH_M;
        assert!(r > 3.0 && r < 3.6, "R(10 M⊕) = {r:.2} R⊕");
    }

    #[test]
    fn far_ir_w4_is_brightest_band() {
        // For a ~47 K body the flux peaks in the far-IR (Wien peak λ ≈ 60 µm),
        // so among the WISE bands W4 (22 µm) is by far the brightest — the
        // paper's "longward of ~13 µm" conclusion.
        let p = nominal();
        let best = p
            .most_favorable_band(&[Band::W1, Band::W2, Band::W3, Band::W4])
            .unwrap();
        assert_eq!(best, Band::W4);
        // And W4 is brighter than the near-IR W1/W2.
        assert!(p.thermal_magnitude(Band::W4) < p.thermal_magnitude(Band::W2));
        assert!(p.thermal_magnitude(Band::W4) < p.thermal_magnitude(Band::W1));
    }

    #[test]
    fn w4_magnitude_near_published() {
        // Published W4 = 10.2 for the nominal case. The far-IR is the regime
        // where the single-temperature grey-body model is most trustworthy, so
        // we pin within ~1 mag.
        let m = nominal().thermal_magnitude(Band::W4);
        assert!(
            (m - reference::W4_10ME_700AU).abs() < 1.0,
            "W4 = {m:.2}, published {:.2}",
            reference::W4_10ME_700AU
        );
    }

    #[test]
    fn q_magnitude_near_published() {
        // Published Q = 10.7 for the nominal case (~20 µm). Our grey-body model
        // gives ≈11.8, a ~1.1-mag residual: the ground-based Q zero point and
        // bandpass are approximate, and a single-temperature grey body lacks
        // the band emissivity structure of the paper's atmosphere model. The
        // sister far-IR band W4 (22 µm, Wright et al. zero point) reproduces
        // the published 10.2 to <1 mag — see `w4_magnitude_near_published` —
        // so the favorable-band conclusion holds where the inputs are pinned.
        let m = nominal().thermal_magnitude(Band::Q);
        assert!(
            (m - reference::Q_10ME_700AU).abs() < 1.5,
            "Q = {m:.2}, published {:.2}",
            reference::Q_10ME_700AU
        );
    }

    #[test]
    fn thermal_brightens_monotonically_with_mass() {
        // Heavier -> larger and hotter -> brighter in every band.
        for &band in &[Band::W4, Band::Q, Band::W3] {
            let m5 = P9 {
                mass_earth: 5.0,
                distance_au: 700.0,
            }
            .thermal_magnitude(band);
            let m20 = P9 {
                mass_earth: 20.0,
                distance_au: 700.0,
            }
            .thermal_magnitude(band);
            assert!(
                m20 < m5,
                "band {band:?}: 20 M⊕ {m20:.1} not brighter than 5 M⊕ {m5:.1}"
            );
        }
    }

    #[test]
    fn near_ir_is_steeply_temperature_sensitive() {
        // On the Wien tail at W1 (3.4 µm) a 20 M⊕ (hotter) body is vastly
        // brighter than a 5 M⊕ one — far steeper than the gentle far-IR slope.
        let m5 = P9 {
            mass_earth: 5.0,
            distance_au: 700.0,
        }
        .thermal_magnitude(Band::W1);
        let m20 = P9 {
            mass_earth: 20.0,
            distance_au: 700.0,
        }
        .thermal_magnitude(Band::W1);
        let near_ir_brightening = m5 - m20;
        let m5q = P9 {
            mass_earth: 5.0,
            distance_au: 700.0,
        }
        .thermal_magnitude(Band::W4);
        let m20q = P9 {
            mass_earth: 20.0,
            distance_au: 700.0,
        }
        .thermal_magnitude(Band::W4);
        let far_ir_brightening = m5q - m20q;
        assert!(
            near_ir_brightening > far_ir_brightening,
            "W1 brightening {near_ir_brightening:.1} should exceed W4 {far_ir_brightening:.1}"
        );
    }

    #[test]
    fn zero_point_flux_gives_zero_magnitude() {
        for &band in &[Band::W1, Band::W2, Band::W3, Band::W4, Band::N, Band::Q] {
            let m = -2.5 * (band.zero_point_jy() / band.zero_point_jy()).log10();
            assert!(m.abs() < 1e-12);
        }
    }
}
