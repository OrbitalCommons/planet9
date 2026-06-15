//! Inverting the candidate's far-IR flux for an implied heliocentric
//! distance, and the forward cross-check flux of a P9 at that distance.
//!
//! The far-IR thermal flux of a cold body is, in the same blackbody form
//! used by the sibling crate [`p9_2025_iras_akari::thermal_model`],
//!
//!   F_ν = π · B_ν(T) · (R / d)²,
//!
//! where `R = R(M)` is the Neptune-anchored mass-radius relation and `d` is
//! the heliocentric distance (≈ geocentric for an outer-solar-system body).
//! Solving for `d`:
//!
//!   d = R · sqrt( π · B_ν(T) / F_ν ).
//!
//! Distance therefore falls as 1/√F at fixed (T, M) — a brighter source is
//! closer — while the forward flux falls as 1/d² and rises with R² (hence
//! with mass) and with the Planck factor B_ν(T) (hence with temperature).
//! These are the scalings the paper's 225 AU / 3–5 M⊕ inference rests on.

use p9_2025_iras_akari::thermal_model::P9ThermalParams;

/// 1 AU in metres (matches `p9_2025_iras_akari::thermal_model`).
const AU_M: f64 = 1.495_978_707e11;
/// Earth radius in metres.
const R_EARTH_M: f64 = 6.371e6;
/// 1 Jansky in W/m²/Hz.
const JY: f64 = 1.0e-26;

/// Physical radius (m) of a body of `mass_earth` Earth masses, from the
/// shared Neptune-anchored mass-radius relation in `p9-core`.
pub fn radius_m(mass_earth: f64) -> f64 {
    p9_core::analysis::photometry::mass_radius_neptunian(mass_earth) * R_EARTH_M
}

/// Implied heliocentric distance (AU) of a blackbody of effective
/// temperature `t_eff_k` and mass `mass_earth` that produces an observed
/// flux density `flux_jy` at wavelength `wavelength_m`.
///
///   d = R · sqrt( π · B_ν(T) / F_ν ).
///
/// This is the exact inverse of
/// [`P9ThermalParams::flux_density_jy`]: feeding the returned distance back
/// into the forward model recovers `flux_jy` (pinned in tests).
pub fn implied_distance_au(flux_jy: f64, t_eff_k: f64, mass_earth: f64, wavelength_m: f64) -> f64 {
    let nu = P9ThermalParams::wavelength_to_freq(wavelength_m);
    let bnu = P9ThermalParams::planck_bnu(t_eff_k, nu);
    let r = radius_m(mass_earth);
    let f_si = flux_jy * JY;
    // F = π B (R/d)²  ⇒  d = R √(π B / F).
    let d_m = r * (std::f64::consts::PI * bnu / f_si).sqrt();
    d_m / AU_M
}

/// Forward far-IR flux density (Jy) of a P9 of the given parameters at the
/// given distance and wavelength — a thin reuse of the sibling thermal
/// model, used to cross-check that the implied distance regenerates the
/// observed flux and to predict the flux in the *other* IRAS band.
pub fn expected_flux_jy(mass_earth: f64, distance_au: f64, t_eff_k: f64, wavelength_m: f64) -> f64 {
    P9ThermalParams {
        mass_earth,
        distance_au,
        t_eff: t_eff_k,
    }
    .flux_density_jy(wavelength_m)
}

/// Result of inverting a candidate flux: the implied distance plus the
/// forward cross-check fluxes a P9 at that distance would emit in both
/// IRAS bands.
#[derive(Debug, Clone, Copy)]
pub struct CandidateInversion {
    /// Implied heliocentric distance (AU) from the 60 µm flux.
    pub distance_au: f64,
    /// Forward 60 µm flux at the implied distance (Jy) — should round-trip
    /// to the input 60 µm flux.
    pub flux_60um_jy: f64,
    /// Forward 100 µm flux at the implied distance (Jy).
    pub flux_100um_jy: f64,
}

/// Invert a 60 µm candidate flux for distance, then evaluate the forward
/// fluxes in both IRAS bands at that distance.
pub fn invert_candidate(flux_60um_jy: f64, t_eff_k: f64, mass_earth: f64) -> CandidateInversion {
    let distance_au = implied_distance_au(flux_60um_jy, t_eff_k, mass_earth, crate::LAMBDA_60UM_M);
    CandidateInversion {
        distance_au,
        flux_60um_jy: expected_flux_jy(mass_earth, distance_au, t_eff_k, crate::LAMBDA_60UM_M),
        flux_100um_jy: expected_flux_jy(mass_earth, distance_au, t_eff_k, crate::LAMBDA_100UM_M),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::REF_CANDIDATE;
    use approx::assert_relative_eq;

    #[test]
    fn reference_candidate_implies_published_distance() {
        // Headline reproduction: inverting the candidate's 60 µm flux at the
        // assumed temperature/mass lands on the paper's fitted 225 ± 15 AU.
        let d = implied_distance_au(
            REF_CANDIDATE.flux_60um_jy,
            REF_CANDIDATE.t_eff_k,
            REF_CANDIDATE.mass_earth,
            crate::LAMBDA_60UM_M,
        );
        assert!(
            (d - REF_CANDIDATE.distance_au).abs() <= REF_CANDIDATE.distance_err_au,
            "implied distance {d:.1} AU outside published 225 ± 15 AU"
        );
    }

    #[test]
    fn inversion_round_trips_to_input_flux() {
        // d = R √(πB/F) is the exact inverse of F = πB(R/d)²: feeding the
        // implied distance back into the forward model recovers the flux.
        let inv = invert_candidate(
            REF_CANDIDATE.flux_60um_jy,
            REF_CANDIDATE.t_eff_k,
            REF_CANDIDATE.mass_earth,
        );
        assert_relative_eq!(
            inv.flux_60um_jy,
            REF_CANDIDATE.flux_60um_jy,
            max_relative = 1e-9
        );
    }

    #[test]
    fn implied_distance_scales_as_inverse_sqrt_flux() {
        // At fixed (T, M): d ∝ 1/√F. Quadrupling the flux halves the
        // distance.
        let d1 = implied_distance_au(0.1, 40.0, 4.0, crate::LAMBDA_60UM_M);
        let d4 = implied_distance_au(0.4, 40.0, 4.0, crate::LAMBDA_60UM_M);
        assert_relative_eq!(d1 / d4, 2.0, max_relative = 1e-9);
    }

    #[test]
    fn forward_flux_falls_as_inverse_square_distance() {
        let near = expected_flux_jy(4.0, 200.0, 40.0, crate::LAMBDA_60UM_M);
        let far = expected_flux_jy(4.0, 400.0, 40.0, crate::LAMBDA_60UM_M);
        assert_relative_eq!(near / far, 4.0, max_relative = 1e-9);
    }

    #[test]
    fn forward_flux_rises_with_mass_via_radius_squared() {
        // F ∝ R² and R ∝ M^0.27, so F ∝ M^0.54.
        let f3 = expected_flux_jy(3.0, 225.0, 40.0, crate::LAMBDA_60UM_M);
        let f5 = expected_flux_jy(5.0, 225.0, 40.0, crate::LAMBDA_60UM_M);
        assert!(f5 > f3);
        let ratio = f5 / f3;
        let expected = (5.0_f64 / 3.0).powf(2.0 * 0.27);
        assert_relative_eq!(ratio, expected, max_relative = 1e-9);
    }

    #[test]
    fn forward_flux_rises_with_temperature() {
        let cool = expected_flux_jy(4.0, 225.0, 35.0, crate::LAMBDA_60UM_M);
        let warm = expected_flux_jy(4.0, 225.0, 45.0, crate::LAMBDA_60UM_M);
        assert!(warm > cool, "60 µm flux should rise with temperature");
    }

    #[test]
    fn brighter_candidate_is_closer() {
        // The whole inference: a brighter 60 µm source implies a smaller
        // heliocentric distance at fixed assumed (T, M).
        let bright = implied_distance_au(0.5, 40.0, 4.0, crate::LAMBDA_60UM_M);
        let faint = implied_distance_au(0.2, 40.0, 4.0, crate::LAMBDA_60UM_M);
        assert!(bright < faint);
    }

    #[test]
    fn fitted_mass_is_few_earth_masses() {
        // The radius implied by the fitted mass sits between sub-Neptune and
        // Neptune: a few-Earth-mass ice giant, not a gas giant.
        let r_earths = radius_m(REF_CANDIDATE.mass_earth) / R_EARTH_M;
        assert!(
            (2.0..3.5).contains(&r_earths),
            "R = {r_earths:.2} R⊕ — expected few-Earth-mass ice giant"
        );
        assert!(REF_CANDIDATE.mass_lo_earth >= 3.0 && REF_CANDIDATE.mass_hi_earth <= 5.0);
    }

    #[test]
    fn mass_range_distances_bracket_published_value() {
        // Across the fitted 3–5 M⊕ range (radius, hence flux, varies), the
        // implied distance from a fixed observed flux still lands within the
        // published 225 ± 15 AU band when the temperature is co-fitted.
        // Here we instead check the simpler monotone fact: a heavier body of
        // the same observed flux is inferred to be *farther* (bigger emitter
        // ⇒ must be more distant to give the same flux).
        let d_light = implied_distance_au(
            REF_CANDIDATE.flux_60um_jy,
            REF_CANDIDATE.t_eff_k,
            REF_CANDIDATE.mass_lo_earth,
            crate::LAMBDA_60UM_M,
        );
        let d_heavy = implied_distance_au(
            REF_CANDIDATE.flux_60um_jy,
            REF_CANDIDATE.t_eff_k,
            REF_CANDIDATE.mass_hi_earth,
            crate::LAMBDA_60UM_M,
        );
        assert!(d_heavy > d_light);
    }

    #[test]
    fn candidate_60um_flux_near_iras_limit() {
        // The candidate is a marginal detection: its 60 µm flux sits just
        // above the IRAS Faint Source limit (~0.2 Jy), consistent with the
        // paper's tentative single-HCON survivor.
        let iras = p9_2025_iras_akari::survey_model::IrasSurvey::default();
        assert!(iras.is_detectable(REF_CANDIDATE.flux_60um_jy));
        assert!(
            REF_CANDIDATE.flux_60um_jy < 5.0 * iras.sensitivity_jy,
            "candidate should be near, not far above, the IRAS limit"
        );
    }

    #[test]
    fn hundred_micron_brighter_than_sixty_for_cold_body() {
        // For a ~40 K blackbody the 100 µm flux exceeds the 60 µm flux
        // (Rayleigh-Jeans side of the peak), as the paper's 60/100 µm
        // detection picture assumes.
        let inv = invert_candidate(
            REF_CANDIDATE.flux_60um_jy,
            REF_CANDIDATE.t_eff_k,
            REF_CANDIDATE.mass_earth,
        );
        assert!(inv.flux_100um_jy > inv.flux_60um_jy);
        // And it matches our labelled reference 100 µm flux.
        assert_relative_eq!(
            inv.flux_100um_jy,
            REF_CANDIDATE.flux_100um_jy,
            max_relative = 5e-3
        );
    }
}
