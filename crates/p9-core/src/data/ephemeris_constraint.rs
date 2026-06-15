//! Cassini / Iorio favored true-anomaly zone for Planet Nine.
//!
//! Cassini Earth–Saturn range residuals (Fienga et al. 2016) and Saturn-
//! precession bounds favor Planet Nine being near true anomaly `ν ≈ 108–129°`
//! on the Brown & Batygin (2016) reference orbit *right now* (and disfavor
//! large complementary arcs). This module is the shared, documented home for
//! that reference orbit and the favored-ν constants; the paper-specific range-
//! residual physics lives in the reproduction crates.

use crate::constants::DEG2RAD;
use crate::types::P9Params;

/// Planet Nine on the Brown & Batygin (2016) orbit, exactly as tabulated by
/// Fienga et al. (2016), Table 1: `a = 700 AU`, `e = 0.6`, `i = 30°`,
/// `ω = 150°`, `Ω = 113°` (IERS ecliptic). The mean anomaly is irrelevant for
/// the favored-ν analysis — callers sweep the true anomaly explicitly. Note
/// `Ω = 113°` differs from `P9Params::nominal_2016()` (`Ω = 100°`), so this
/// paper-reference orbit is built directly.
pub fn brown_batygin_orbit() -> P9Params {
    P9Params {
        mass_earth: P9_MASS_EARTH,
        a: 700.0,
        e: 0.6,
        i: 30.0 * DEG2RAD,
        omega: 150.0 * DEG2RAD,
        omega_big: 113.0 * DEG2RAD,
        mean_anomaly: 0.0,
    }
}

/// Planet Nine mass assumed by the Brown & Batygin orbit (Earth masses).
pub const P9_MASS_EARTH: f64 = 10.0;

/// Most-probable true anomaly of P9 that REDUCES the Cassini Earth–Saturn
/// range residuals (degrees): `v = 117.8°₋₁₀⁺¹¹` (Fienga et al. 2016, §4).
pub const PREFERRED_TRUE_ANOMALY_DEG: f64 = 117.8;

/// Favored interval for the true anomaly (degrees): `v ∈ [108°, 129°]`
/// (Fienga et al. 2016, Conclusions / Fig 6 green zone).
pub const FAVORED_INTERVAL_DEG: (f64, f64) = (108.0, 129.0);

/// Tolerance (deg) within which an analytic-proxy favored anomaly should land
/// of the published preferred true anomaly.
pub const TRUE_ANOMALY_TOLERANCE_DEG: f64 = 30.0;

/// Published reference values from Siraj, Chyba & Tremaine (2024),
/// "Orbit of a Possible Planet X" (arXiv:2410.18170). This independent
/// inference yields a markedly different perturber than Brown & Batygin (2021):
/// a lower mass, a smaller semi-major axis, and a lower inclination.
/// Inferred perturber mass (Earth masses).
pub const SIRAJ_2024_MASS_EARTH: f64 = 4.4;

/// Inferred perturber semi-major axis (AU).
pub const SIRAJ_2024_A_AU: f64 = 290.0;

/// Inferred perturber inclination (degrees).
pub const SIRAJ_2024_I_DEG: f64 = 6.8;

/// Approximate 1σ uncertainty on the inferred mass (Earth masses). Used as the
/// Siraj-side scatter when computing the tension against Brown & Batygin 2021.
pub const SIRAJ_2024_MASS_SIGMA_EARTH: f64 = 1.0;

/// Approximate 1σ uncertainty on the inferred semi-major axis (AU).
pub const SIRAJ_2024_A_SIGMA_AU: f64 = 30.0;

/// Approximate 1σ uncertainty on the inferred inclination (degrees).
pub const SIRAJ_2024_I_SIGMA_DEG: f64 = 1.7;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_siraj_2024_reference_constants() {
        assert_eq!(SIRAJ_2024_MASS_EARTH, 4.4);
        assert_eq!(SIRAJ_2024_A_AU, 290.0);
        assert!((SIRAJ_2024_I_DEG - 6.8).abs() < 1e-9);
        assert_eq!(SIRAJ_2024_MASS_SIGMA_EARTH, 1.0);
        assert_eq!(SIRAJ_2024_A_SIGMA_AU, 30.0);
        assert!((SIRAJ_2024_I_SIGMA_DEG - 1.7).abs() < 1e-9);
    }
}
