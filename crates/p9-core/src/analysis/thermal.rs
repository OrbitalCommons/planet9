//! Shared blackbody / grey-body thermal model primitives for the IR–mm Planet
//! Nine detectability reproductions.
//!
//! Several crates model the thermal (and reflected-sunlight) emission of a cold
//! Planet Nine across the infrared and millimeter bands. They share the same
//! underlying physics — the Planck function, an effective-temperature energy
//! balance (internal-luminosity floor vs. solar equilibrium), the
//! solid-angle × brightness flux law, the conversion of a flux density to a
//! Vega-system magnitude, and a bisection for the maximum detectable distance.
//!
//! This module hosts the **generic, parameterized** versions of those
//! primitives. Every crate-specific number — band wavelength, Vega zero point,
//! geometric albedo, internal-temperature floor, mass–radius anchoring — stays
//! in the calling crate and is passed in. Only the band-agnostic math lives
//! here.

use std::f64::consts::PI;

/// Planck constant (J·s).
pub const H_PLANCK: f64 = 6.626_070_15e-34;
/// Boltzmann constant (J/K).
pub const K_BOLTZ: f64 = 1.380_649e-23;
/// Speed of light (m/s).
pub const C_LIGHT: f64 = 2.997_924_58e8;
/// Solar radius (m).
pub const R_SUN_M: f64 = 6.957e8;
/// Solar effective temperature (K).
pub const T_SUN: f64 = 5778.0;
/// Solar bolometric luminosity (W).
pub const L_SUN_W: f64 = 3.828e26;
/// 1 AU in meters.
pub const AU_M: f64 = 1.495_978_707e11;
/// 1 Jansky in W/m²/Hz.
pub const JY: f64 = 1.0e-26;

/// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
///
/// The dimensionless argument x = hν/kT is guarded: for x > 700 (deep on the
/// Wien tail) the spectral radiance underflows to 0, which also avoids
/// overflowing `exp`.
pub fn planck_bnu(t: f64, nu_hz: f64) -> f64 {
    let x = H_PLANCK * nu_hz / (K_BOLTZ * t);
    if x > 700.0 {
        return 0.0;
    }
    (2.0 * H_PLANCK * nu_hz.powi(3) / (C_LIGHT * C_LIGHT)) / (x.exp() - 1.0)
}

/// Rayleigh–Jeans brightness B_ν ≈ 2 ν² k T / c² (W/m²/Hz/sr). The low-energy
/// (hν ≪ kT) limit of [`planck_bnu`]; explicitly ∝ ν²·T.
pub fn rayleigh_jeans_bnu(t: f64, nu_hz: f64) -> f64 {
    2.0 * nu_hz * nu_hz * K_BOLTZ * t / (C_LIGHT * C_LIGHT)
}

/// Dimensionless photon-energy parameter x = hν/kT. The Rayleigh–Jeans regime
/// is x ≪ 1, the Wien regime x ≫ 1.
pub fn photon_energy_ratio(t: f64, nu_hz: f64) -> f64 {
    H_PLANCK * nu_hz / (K_BOLTZ * t)
}

/// Solar equilibrium temperature (K) for a fast rotator radiating from its full
/// surface at heliocentric distance `distance_au` with geometric/Bond albedo
/// `albedo`: T_eq = T_sun √(R_sun / 2d) (1 − A)^¼.
pub fn solar_equilibrium_temp(distance_au: f64, albedo: f64) -> f64 {
    let d_m = distance_au * AU_M;
    T_SUN * (R_SUN_M / (2.0 * d_m)).sqrt() * (1.0 - albedo).powf(0.25)
}

/// Effective temperature (K): the larger of the internal-heat floor
/// `internal_temp_k` and the solar-equilibrium temperature at `distance_au`
/// with `albedo`. Beyond a few hundred AU the solar term is ~10 K, so the
/// internal floor dominates.
pub fn effective_temp(distance_au: f64, albedo: f64, internal_temp_k: f64) -> f64 {
    solar_equilibrium_temp(distance_au, albedo).max(internal_temp_k)
}

/// Thermal grey-body flux density (Jy) of a unit-emissivity emitter of physical
/// radius `radius_m` at heliocentric distance `distance_au`, evaluated at
/// frequency `nu_hz`: F_ν = π B_ν(T) (R/d)².
pub fn thermal_flux_jy(t: f64, radius_m: f64, distance_au: f64, nu_hz: f64) -> f64 {
    let bnu = planck_bnu(t, nu_hz);
    let d = distance_au * AU_M;
    let solid_angle = PI * (radius_m / d).powi(2);
    solid_angle * bnu / JY
}

/// Reflected-sunlight flux density (Jy) of a body of radius `radius_m` and
/// reflectance `albedo` at heliocentric distance `distance_au`, observed near
/// opposition (Δ = d − 1 AU, floored at 1 AU), evaluated at frequency `nu_hz`.
///
/// The Sun is a blackbody at [`T_SUN`]; the solar flux at the planet is
/// π B_ν(T_sun) (R_sun / r)². The planet intercepts a disc πR_p², reflects a
/// fraction `albedo` into ~2π sr (Lambert factor absorbed into the geometric
/// albedo), and that reradiated flux falls as 1/Δ² to the observer:
///
///   F_refl = albedo · [π B_ν(T_sun) (R_sun/r)²] · (R_p² / (4 Δ²)).
pub fn reflected_flux_jy(albedo: f64, radius_m: f64, distance_au: f64, nu_hz: f64) -> f64 {
    let bnu_sun = planck_bnu(T_SUN, nu_hz);
    let solar_flux_at_planet = PI * bnu_sun * (R_SUN_M / (distance_au * AU_M)).powi(2);
    let delta_m = (distance_au - 1.0).max(1.0) * AU_M;
    let reflected =
        albedo * solar_flux_at_planet * (radius_m * radius_m) / (4.0 * delta_m * delta_m);
    reflected / JY
}

/// Apparent Vega-system magnitude from a flux density `flux` and a band Vega
/// zero-point flux `zero_point` (same units). Fainter (larger) numbers for
/// smaller flux; +∞ when the flux is non-positive.
pub fn flux_to_magnitude(flux: f64, zero_point: f64) -> f64 {
    if flux <= 0.0 {
        return f64::INFINITY;
    }
    -2.5 * (flux / zero_point).log10()
}

/// Maximum heliocentric distance (AU) at which `magnitude_at(distance)` reaches
/// `depth`, found by bisection on the monotone (fainting-with-distance)
/// magnitude–distance relation over `[d_min, d_max]`.
///
/// Returns `None` if the source is already fainter than `depth` at `d_min`
/// (undetectable anywhere in range); returns `Some(d_max)` if it is still
/// brighter than `depth` at `d_max` (detectable to the far edge of the range).
pub fn max_detectable_distance(
    d_min: f64,
    d_max: f64,
    depth: f64,
    magnitude_at: impl Fn(f64) -> f64,
) -> Option<f64> {
    if magnitude_at(d_min) > depth {
        return None;
    }
    if magnitude_at(d_max) <= depth {
        return Some(d_max);
    }

    let mut lo = d_min;
    let mut hi = d_max;
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if magnitude_at(mid) <= depth {
            lo = mid;
        } else {
            hi = mid;
        }
        if hi - lo < 1e-4 {
            break;
        }
    }
    Some(0.5 * (lo + hi))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn planck_underflows_on_deep_wien_tail() {
        // A 40 K body at WISE W1 (3.3526 µm) has x ≈ 107 ≫ 700? No — but the
        // guard returns exactly 0 once x exceeds 700, and the flux is otherwise
        // a tiny positive number on the Wien tail.
        let nu_w1 = C_LIGHT / 3.3526e-6;
        let b = planck_bnu(40.0, nu_w1);
        assert!(b >= 0.0 && b < 1e-40, "B_nu(40 K, W1) = {b:.3e}");
        // Hard guard: an extreme x returns exactly 0.
        assert_eq!(planck_bnu(1.0, nu_w1), 0.0);
    }

    #[test]
    fn rayleigh_jeans_matches_planck_at_low_x() {
        // At 280 GHz and 40 K, x ≈ 0.34: the RJ approximation tracks the exact
        // Planck function to better than 20%.
        let nu = 280e9;
        let t = 40.0;
        assert!(photon_energy_ratio(t, nu) < 0.5);
        let exact = planck_bnu(t, nu);
        let rj = rayleigh_jeans_bnu(t, nu);
        assert!((exact - rj).abs() / rj < 0.2);
    }

    #[test]
    fn solar_equilibrium_is_cold_beyond_a_few_hundred_au() {
        // At 500 AU the solar-equilibrium temperature is ~10 K, below a 40 K
        // internal floor, so the effective temperature is the floor.
        assert!(solar_equilibrium_temp(500.0, 0.41) < 15.0);
        assert!((effective_temp(500.0, 0.41, 40.0) - 40.0).abs() < 1e-9);
    }

    #[test]
    fn flux_falls_as_inverse_distance_squared() {
        let near = thermal_flux_jy(40.0, 2.0e7, 400.0, 150e9);
        let far = thermal_flux_jy(40.0, 2.0e7, 800.0, 150e9);
        assert!((near / far - 4.0).abs() < 1e-9);
    }

    #[test]
    fn flux_scales_with_radius_squared() {
        let a = thermal_flux_jy(40.0, 1.0e7, 500.0, 150e9);
        let b = thermal_flux_jy(40.0, 2.0e7, 500.0, 150e9);
        assert!((b / a - 4.0).abs() < 1e-9);
    }

    #[test]
    fn zero_point_gives_zero_magnitude() {
        let zp = 309.540;
        assert!(flux_to_magnitude(zp, zp).abs() < 1e-12);
    }

    #[test]
    fn magnitude_is_infinite_for_zero_flux() {
        assert!(flux_to_magnitude(0.0, 309.540).is_infinite());
    }

    #[test]
    fn max_detectable_distance_brackets_the_crossover() {
        // A toy magnitude that faints by 5 log10(d): pick a depth and confirm
        // the returned distance reproduces it.
        let depth = 16.5;
        let mag = |d: f64| 5.0 * d.log10();
        let d = max_detectable_distance(50.0, 50_000.0, depth, mag).expect("detectable");
        assert!((mag(d) - depth).abs() < 1e-3);
    }

    #[test]
    fn max_detectable_distance_none_when_undetectable_at_near_edge() {
        let mag = |d: f64| 5.0 * d.log10();
        // Depth far brighter than even the near-edge magnitude → never detectable.
        assert!(max_detectable_distance(50.0, 50_000.0, 1.0, mag).is_none());
    }
}
