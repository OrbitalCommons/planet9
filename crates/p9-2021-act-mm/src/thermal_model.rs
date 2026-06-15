//! Millimeter-wave thermal emission of a cold Planet Nine.
//!
//! Naess et al. (2021) point out that a cold (~30–50 K) Planet Nine, although
//! hopeless to detect thermally in the near-IR WISE bands (a 40 K body sits far
//! down the Wien tail at 3.4 µm; see `p9-2018-wise-search`), is comparatively
//! *bright* at millimeter wavelengths. At ACT's 150 GHz the photon energy is
//!
//!   hν/kT = h(150 GHz)/(k·40 K) ≈ 0.18,
//!
//! so the body is firmly on the **Rayleigh–Jeans** tail (hν ≪ kT), where the
//! Planck function reduces to B_ν ≈ 2 ν² k T / c² and the flux scales as ν²·T
//! rather than the exp(−hν/kT) cliff of the IR. The flux density of an
//! unresolved grey body of physical radius R at distance d is
//!
//!   F_ν = π B_ν(T) (R/d)²        [W/m²/Hz],
//!
//! the same solid-angle × Planck-brightness law used by the WISE crate's
//! [`P9Thermal::thermal_flux_w1_jy`]. This module reuses that crate's
//! [`P9Thermal`] for the cold-body effective temperature (internal-heat floor
//! dominating solar heating beyond a few hundred AU) and the Neptune-anchored
//! mass–radius relation, and evaluates the same Planck law at mm frequencies.

use std::f64::consts::PI;

use p9_2018_wise_search::thermal_model::P9Thermal;
use p9_core::constants::AU_M;

/// Planck constant (J·s).
pub const H_PLANCK: f64 = 6.626_070_15e-34;
/// Boltzmann constant (J/K).
pub const K_BOLTZ: f64 = 1.380_649e-23;
/// Speed of light (m/s).
pub const C_LIGHT: f64 = 2.997_924_58e8;

/// 1 Jansky in W/m²/Hz.
pub const JY: f64 = 1.0e-26;
/// 1 milliJansky in W/m²/Hz.
pub const MJY: f64 = 1.0e-29;

/// Planck function B_ν(T) in W/m²/Hz/sr at frequency `nu_hz`.
///
/// This is the exact Planck law; at ACT mm frequencies and ~40 K it is well
/// approximated by [`rayleigh_jeans_bnu`] (the two agree to a few percent).
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

/// A cold Planet Nine viewed at millimeter wavelengths. Wraps the WISE crate's
/// [`P9Thermal`] (mass, distance, internal-temperature floor, Neptune-anchored
/// radius) so the temperature and radius physics is shared, not duplicated.
#[derive(Debug, Clone, Copy)]
pub struct P9Millimeter {
    /// Underlying thermal body (sets effective temperature and radius).
    pub body: P9Thermal,
}

impl P9Millimeter {
    /// A Planet Nine of `mass_earth` at heliocentric `distance_au` with the
    /// default cold internal-temperature floor (≈40 K) from the WISE crate.
    pub fn new(mass_earth: f64, distance_au: f64) -> Self {
        Self {
            body: P9Thermal::new(mass_earth, distance_au),
        }
    }

    /// Override the internal-luminosity effective-temperature floor (K).
    pub fn with_temperature(mut self, temp_k: f64) -> Self {
        self.body.internal_temp_k = temp_k;
        self
    }

    /// Effective temperature (K): the larger of solar heating and the internal
    /// floor. At ACT-relevant distances this is the cold internal floor.
    pub fn effective_temp(&self) -> f64 {
        self.body.effective_temp()
    }

    /// Physical radius (m), Neptune-anchored mass–radius relation.
    pub fn radius_m(&self) -> f64 {
        self.body.radius_m()
    }

    /// Thermal mm flux density (W/m²/Hz) at frequency `nu_hz`, unit-emissivity
    /// grey body: F_ν = π B_ν(T) (R/d)².
    pub fn flux_density_si(&self, nu_hz: f64) -> f64 {
        let bnu = planck_bnu(self.effective_temp(), nu_hz);
        let r = self.radius_m();
        let d = self.body.distance_au * AU_M;
        PI * (r / d).powi(2) * bnu
    }

    /// Thermal mm flux density in milliJansky at `nu_hz` — the unit ACT
    /// point-source sensitivities are quoted in.
    pub fn flux_density_mjy(&self, nu_hz: f64) -> f64 {
        self.flux_density_si(nu_hz) / MJY
    }

    /// Flux density in Jansky at `nu_hz`.
    pub fn flux_density_jy(&self, nu_hz: f64) -> f64 {
        self.flux_density_si(nu_hz) / JY
    }

    /// Rayleigh–Jeans approximation to the mm flux density (mJy), using
    /// [`rayleigh_jeans_bnu`]. Differs from [`flux_density_mjy`] by the small
    /// (few-percent) RJ correction at ACT frequencies.
    pub fn flux_density_rj_mjy(&self, nu_hz: f64) -> f64 {
        let bnu = rayleigh_jeans_bnu(self.effective_temp(), nu_hz);
        let r = self.radius_m();
        let d = self.body.distance_au * AU_M;
        PI * (r / d).powi(2) * bnu / MJY
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_2018_wise_search::thermal_model::INTERNAL_TEMP_FLOOR_K;

    /// ACT 150 GHz in Hz.
    const NU_150: f64 = 150.0e9;

    #[test]
    fn cold_body_is_on_the_rayleigh_jeans_tail_at_mm() {
        // At 150 GHz and the 40 K floor, hν/kT ≈ 0.18 ≪ 1: deep in the RJ
        // regime, the opposite of WISE W1 where hν/kT ≈ 107 (Wien cliff).
        let x = photon_energy_ratio(INTERNAL_TEMP_FLOOR_K, NU_150);
        assert!(x < 0.25, "hν/kT = {x:.3} should be ≪ 1 (RJ regime)");

        // The WISE-W1 ratio for the same 40 K body is enormous.
        let nu_w1 = C_LIGHT / 3.3526e-6;
        let x_w1 = photon_energy_ratio(INTERNAL_TEMP_FLOOR_K, nu_w1);
        assert!(x_w1 > 100.0, "W1 hν/kT = {x_w1:.0} should be ≫ 1 (Wien)");
    }

    #[test]
    fn planck_matches_rayleigh_jeans_at_mm() {
        // On the RJ tail the exact Planck and RJ brightness agree to a few
        // percent (the leading correction is (1 − x/2), x ≈ 0.18).
        let p = P9Millimeter::new(10.0, 500.0);
        let exact = p.flux_density_mjy(NU_150);
        let rj = p.flux_density_rj_mjy(NU_150);
        let rel = (rj - exact).abs() / exact;
        assert!(rel < 0.12, "RJ vs Planck relative diff = {rel:.3}");
        // RJ slightly over-predicts (it omits the −x/2 term).
        assert!(rj > exact);
    }

    #[test]
    fn flux_follows_rayleigh_jeans_nu_squared_times_t() {
        // RJ brightness ∝ ν²·T at fixed geometry. Doubling ν quadruples flux;
        // doubling T doubles it.
        let p = P9Millimeter::new(10.0, 500.0);
        let f1 = p.flux_density_rj_mjy(NU_150);
        let f2 = p.flux_density_rj_mjy(2.0 * NU_150);
        assert!((f2 / f1 - 4.0).abs() < 1e-9, "ν² scaling: {:.4}", f2 / f1);

        let cold = P9Millimeter::new(10.0, 500.0).with_temperature(30.0);
        let warm = P9Millimeter::new(10.0, 500.0).with_temperature(60.0);
        let ratio = warm.flux_density_rj_mjy(NU_150) / cold.flux_density_rj_mjy(NU_150);
        assert!((ratio - 2.0).abs() < 1e-9, "T scaling: {ratio:.4}");
    }

    #[test]
    fn flux_falls_as_inverse_distance_squared() {
        let near = P9Millimeter::new(10.0, 400.0).flux_density_mjy(NU_150);
        let far = P9Millimeter::new(10.0, 800.0).flux_density_mjy(NU_150);
        // Double the distance → quarter the flux (radius/temperature fixed).
        assert!(
            (near / far - 4.0).abs() < 1e-6,
            "1/d² scaling: {:.4}",
            near / far
        );
    }

    #[test]
    fn flux_scales_with_radius_squared() {
        // F ∝ R². The Neptune relation gives R ∝ M^0.27, so flux ∝ M^0.54;
        // here we verify the explicit R² law by comparing equal-T bodies of
        // different radii at fixed distance.
        let a = P9Millimeter::new(10.0, 500.0).with_temperature(40.0);
        let b = P9Millimeter::new(20.0, 500.0).with_temperature(40.0);
        let ratio_flux = b.flux_density_mjy(NU_150) / a.flux_density_mjy(NU_150);
        let ratio_r2 = (b.radius_m() / a.radius_m()).powi(2);
        assert!(
            (ratio_flux - ratio_r2).abs() < 1e-6,
            "flux ratio {ratio_flux:.4} should equal (R_b/R_a)² = {ratio_r2:.4}"
        );
    }

    #[test]
    fn mm_flux_of_nominal_p9_is_a_few_mjy() {
        // A 10 M⊕, 40 K Planet Nine at 500 AU has a 150 GHz flux of order a
        // few mJy — comparable to ACT's 4–12 mJy point-source limit, which is
        // exactly why the mm search is feasible while WISE is not.
        let f = P9Millimeter::new(10.0, 500.0).flux_density_mjy(NU_150);
        assert!(
            (0.5..50.0).contains(&f),
            "F_150(10 M⊕, 500 AU) = {f:.2} mJy"
        );
    }
}
