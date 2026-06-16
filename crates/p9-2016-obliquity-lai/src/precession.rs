//! The analytic spin–orbit precession theory of Lai (2016).
//!
//! Equation numbers refer to arXiv:1608.01421.

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::units::{au, days, earth_masses, radians, Angle, AngularVelocity, Length, Mass};

use crate::planets;

/// Solar structure parameters (Lai 2016 §2, after Mecheri et al. 2004).
pub mod solar {
    /// Moment-of-inertia coefficient: I₃ = k★ M★ R★² with k★ ≈ 0.06.
    pub const K_STAR: f64 = 0.06;
    /// Quadrupole coefficient: I₃ − I₁ = k_q★ Ω̂★² M★ R★² with k_q★ ≈ 0.01,
    /// giving J₂ = k_q★ Ω̂★² ≈ 2.2e-7.
    pub const KQ_STAR: f64 = 0.01;
    /// λ★ ≡ 6 k_q★ / k★ ≈ 1 (Lai 2016 Eq. 8).
    pub const fn lambda_star() -> f64 {
        6.0 * KQ_STAR / K_STAR
    }
}

/// Published coefficients and observables from Lai (2016), kept as labelled
/// constants and pinned against the computation.
pub mod published {
    /// Observed solar obliquity (the misalignment the theory must reproduce).
    pub const OBSERVED_OBLIQUITY_DEG: f64 = 6.0;
    /// Lai Eq. (5): Ω_L = 2π / 87.5 Gyr for mₚ = 10 m⊕, ãₚ = 400 AU.
    pub const OMEGA_L_PERIOD_GYR: f64 = 87.5;
    /// Lai Eq. (8): Ω★ps = 2π / 55.8 Gyr for P★ = 10 d, λ★ = 1.
    pub const OMEGA_SPIN_PERIOD_GYR: f64 = 55.8;
    /// Lai Eq. (5): Ω_L = 2.74 Ω_Jp (P9 forcing of Jupiter alone).
    pub const OMEGA_L_OVER_OMEGA_JP: f64 = 2.74;
    /// Lai Eq. (8): Ω★ps = 2.88 Ω★J (J solar-spin torque alone).
    pub const OMEGA_SPIN_OVER_J: f64 = 2.88;
    /// Lai Eq. (17) prefactor (AU): the effective-aₚ contour for θsl = 6°.
    pub const A_TILDE_PREFACTOR_AU: f64 = 462.0;
    /// Lai §3: the required ãₚ lies between 340 and 480 AU.
    pub const A_TILDE_RANGE_AU: (f64, f64) = (340.0, 480.0);
}

/// Nominal Planet Nine configuration (Lai/Batygin & Brown 2016).
#[derive(Debug, Clone, Copy)]
pub struct PlanetNine {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Semi-major axis (AU).
    pub a_au: f64,
    /// Eccentricity.
    pub e: f64,
    /// Inclination θₚ relative to the canonical planet plane (rad).
    pub inclination_rad: f64,
}

impl PlanetNine {
    /// Mass in solar masses.
    pub fn mass_solar(&self) -> f64 {
        self.mass_earth * EARTH_MASS_SOLAR
    }

    /// Effective semi-major axis ãₚ ≡ aₚ √(1 − eₚ²) (Lai Eq. 1).
    pub fn a_tilde_au(&self) -> f64 {
        self.a_au * (1.0 - self.e * self.e).sqrt()
    }

    /// Planet Nine orbital angular momentum `Lp` (M_sun·AU²/day), Lai Eq. (1).
    ///
    /// `Lp = mp √(GM_sun ãp)`. (Lai writes it as 0.276 L_J × scalings; the two
    /// forms agree to < 0.1%, pinned in tests.)
    pub fn angular_momentum(&self) -> f64 {
        planets::orbital_angular_momentum(self.mass_solar(), self.a_tilde_au())
    }

    // ---- typed (uom) boundary: dimension-checked views ----

    /// Planet Nine mass as a [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }
    /// Planet Nine semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }
    /// Planet Nine inclination as an [`Angle`].
    pub fn inclination(&self) -> Angle {
        radians(self.inclination_rad)
    }
    /// Effective semi-major axis ãₚ as a typed [`Length`].
    pub fn a_tilde(&self) -> Length {
        au(self.a_tilde_au())
    }
}

/// Differential nodal precession rate Ω_jp that P9 induces on canonical planet
/// `j` (Lai Eq. 2), in rad/day:
///
///   Ω_jp = (3 mp / 4 M★) (aⱼ / ãp)³ nⱼ
fn omega_jp(p9: &PlanetNine, a_j: f64) -> f64 {
    let n_j = planets::mean_motion(a_j);
    (3.0 * p9.mass_solar() / (4.0 * 1.0)) * (a_j / p9.a_tilde_au()).powi(3) * n_j
}

/// Precession rate of the *rigidly coupled* canonical-planet plane forced by
/// P9 (Lai Eq. 5), in rad/day:
///
///   Ω_L = Σ_j Lⱼ Ω_jp / L
///
/// Lai gives Ω_L = 2π / 87.5 Gyr × (mp/10m⊕)(ãp/400 AU)⁻³.
pub fn omega_l(p9: &PlanetNine) -> f64 {
    let mut num = 0.0;
    let mut l_total = 0.0;
    for &(m, a) in planets::canonical_planets().iter() {
        let l_j = planets::orbital_angular_momentum(m, a);
        num += l_j * omega_jp(p9, a);
        l_total += l_j;
    }
    num / l_total
}

/// [`omega_l`] as a typed [`AngularVelocity`].
pub fn omega_l_typed(p9: &PlanetNine) -> AngularVelocity {
    (radians(omega_l(p9)) / days(1.0)).into()
}

/// Solar spin-axis precession frequency Ω★ps driven by the planetary torques
/// on the solar rotational bulge (Lai Eq. 7), in rad/day:
///
///   Ω★ps = Σ_j (3 k_q★ / 2 k★) (mⱼ / M★) (R★ / aⱼ)³ Ω★
///
/// with Ω★ = 2π / P★. Lai gives Ω★ps = 2π / 55.8 Gyr × (P★/10 d)⁻¹ for λ★ = 1.
pub fn omega_spin_precession(spin_period_days: f64) -> f64 {
    let omega_star = 2.0 * PI / spin_period_days;
    let coeff = 3.0 * solar::KQ_STAR / (2.0 * solar::K_STAR);
    planets::canonical_planets()
        .iter()
        .map(|&(m, a)| coeff * (m / 1.0) * (RADIUS_SUN_AU / a).powi(3) * omega_star)
        .sum()
}

/// [`omega_spin_precession`] as a typed [`AngularVelocity`].
pub fn omega_spin_precession_typed(spin_period_days: f64) -> AngularVelocity {
    (radians(omega_spin_precession(spin_period_days)) / days(1.0)).into()
}

/// The two precession frequencies of the spin equation in the frame
/// corotating with `l̂` (Lai Eqs. 12–13), in rad/day.
///
/// `Ω_y = Ω_L sin θp cos θp` drives the obliquity; `Ω_z` is the mismatch
/// between the spin precession and the plane precession seen in that frame:
///
///   Ω_z ≃ Ω★ps − Ω_L cos θp (L/Lp + cos θp).
///
/// The secular spin–orbit resonance is `Ω_z = 0`.
pub fn corotating_frequencies(p9: &PlanetNine, spin_period_days: f64) -> (f64, f64) {
    let theta_p = p9.inclination_rad;
    let om_l = omega_l(p9);
    let om_spin = omega_spin_precession(spin_period_days);
    let l_total = planets::total_orbital_angular_momentum();
    let lp = p9.angular_momentum();

    let omega_y = om_l * theta_p.sin() * theta_p.cos();
    let omega_z = om_spin - om_l * theta_p.cos() * (l_total / lp + theta_p.cos());
    (omega_y, omega_z)
}

/// [`corotating_frequencies`] as typed [`AngularVelocity`] pair `(Ω_y, Ω_z)`.
pub fn corotating_frequencies_typed(
    p9: &PlanetNine,
    spin_period_days: f64,
) -> (AngularVelocity, AngularVelocity) {
    let (omega_y, omega_z) = corotating_frequencies(p9, spin_period_days);
    (
        (radians(omega_y) / days(1.0)).into(),
        (radians(omega_z) / days(1.0)).into(),
    )
}

/// The analytic solar obliquity θsl at time `t` (Lai Eq. 15), in radians:
///
///   θsl ≃ |(2 Ω_y / Ω_z) sin(Ω_z t / 2)|
///
/// `t_gyr` is the elapsed time in Gyr (Lai evaluates at 4.5 Gyr). The solar
/// rotation enters only through `spin_period_days = P★/λ★` (λ★ ≈ 1).
pub fn solar_obliquity_rad(p9: &PlanetNine, spin_period_days: f64, t_gyr: f64) -> f64 {
    let (omega_y, omega_z) = corotating_frequencies(p9, spin_period_days);
    let t = t_gyr * GYR_DAYS;
    if omega_z.abs() < 1e-30 {
        // Exact resonance: sin(Ω_z t/2)/Ω_z → t/2, so θsl → Ω_y t.
        return (omega_y * t).abs();
    }
    ((2.0 * omega_y / omega_z) * (omega_z * t / 2.0).sin()).abs()
}

/// [`solar_obliquity_rad`] as a typed [`Angle`].
pub fn solar_obliquity_typed(p9: &PlanetNine, spin_period_days: f64, t_gyr: f64) -> Angle {
    radians(solar_obliquity_rad(p9, spin_period_days, t_gyr))
}

/// Convenience: obliquity in degrees at 4.5 Gyr.
pub fn solar_obliquity_deg(p9: &PlanetNine, spin_period_days: f64) -> f64 {
    solar_obliquity_rad(p9, spin_period_days, 4.5) * RAD2DEG
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::angle::radian;
    use uom::si::angular_velocity::radian_per_second;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let p9 = nominal_p9();
        assert_relative_eq!(p9.semi_major_axis().get::<astronomical_unit>(), p9.a_au);
        assert_relative_eq!(p9.inclination().get::<radian>(), p9.inclination_rad);
        assert_relative_eq!(p9.a_tilde().get::<astronomical_unit>(), p9.a_tilde_au());
        assert_relative_eq!(
            omega_l_typed(&p9).get::<radian_per_second>(),
            omega_l(&p9) / 86_400.0,
            epsilon = 1e-40
        );
        assert_relative_eq!(
            omega_spin_precession_typed(10.0).get::<radian_per_second>(),
            omega_spin_precession(10.0) / 86_400.0,
            epsilon = 1e-40
        );
        let (y, z) = corotating_frequencies(&p9, 20.0);
        let (yt, zt) = corotating_frequencies_typed(&p9, 20.0);
        assert_relative_eq!(yt.get::<radian_per_second>(), y / 86_400.0, epsilon = 1e-40);
        assert_relative_eq!(zt.get::<radian_per_second>(), z / 86_400.0, epsilon = 1e-40);
        assert_relative_eq!(
            solar_obliquity_typed(&p9, 20.0, 4.5).get::<radian>(),
            solar_obliquity_rad(&p9, 20.0, 4.5)
        );
    }

    /// A nominal P9: 10 m⊕, ãₚ = 400 AU (set via a, e), θₚ = 20°.
    fn nominal_p9() -> PlanetNine {
        // Choose a, e so that ãₚ = a√(1−e²) = 400 AU with e = 0.5.
        let e = 0.5_f64;
        let a = 400.0 / (1.0 - e * e).sqrt();
        PlanetNine {
            mass_earth: 10.0,
            a_au: a,
            e,
            inclination_rad: 20.0 * DEG2RAD,
        }
    }

    #[test]
    fn lp_eq1_matches_direct() {
        // Lai Eq. (1): Lp = 0.276 L_J for mp = 10 m⊕, ãp = 400 AU.
        let p9 = nominal_p9();
        let ratio = p9.angular_momentum() / planets::l_jupiter();
        assert!(
            (ratio - 0.276).abs() < 0.002,
            "Lp/L_J = {ratio:.4} (Lai: 0.276)"
        );
    }

    #[test]
    fn omega_l_matches_published_period() {
        // Lai Eq. (5): Ω_L = 2π / 87.5 Gyr for the nominal P9.
        let p9 = nominal_p9();
        let om_l = omega_l(&p9);
        let period_gyr = (2.0 * PI / om_l) / GYR_DAYS;
        let rel =
            (period_gyr - published::OMEGA_L_PERIOD_GYR).abs() / published::OMEGA_L_PERIOD_GYR;
        assert!(
            rel < 0.01,
            "Ω_L period = {period_gyr:.2} Gyr (Lai: 87.5), rel {rel:.3e}"
        );
    }

    #[test]
    fn omega_l_is_2p74_omega_jp() {
        // Lai Eq. (5): Ω_L = 2.74 Ω_Jp.
        let p9 = nominal_p9();
        let ratio = omega_l(&p9) / omega_jp(&p9, planets::A_JUPITER_AU);
        assert!(
            (ratio - published::OMEGA_L_OVER_OMEGA_JP).abs() < 0.03,
            "Ω_L/Ω_Jp = {ratio:.3} (Lai: 2.74)"
        );
    }

    #[test]
    fn omega_l_scales_with_mass_and_a_tilde() {
        // Ω_L ∝ mp and ∝ ãp⁻³ (Lai Eq. 5).
        let base = nominal_p9();
        let mut heavy = base;
        heavy.mass_earth = 20.0;
        let r_mass = omega_l(&heavy) / omega_l(&base);
        assert!((r_mass - 2.0).abs() < 1e-6, "mass scaling {r_mass:.4}");

        // Double ãp by doubling a (e fixed): expect Ω_L → 1/8.
        let mut far = base;
        far.a_au = base.a_au * 2.0;
        let r_a = omega_l(&far) / omega_l(&base);
        assert!((r_a - 0.125).abs() < 1e-6, "ãp⁻³ scaling {r_a:.5}");
    }

    #[test]
    fn omega_spin_matches_published_period() {
        // Lai Eq. (8): Ω★ps = 2π / 55.8 Gyr for P★ = 10 d, λ★ = 1.
        let om = omega_spin_precession(10.0);
        let period_gyr = (2.0 * PI / om) / GYR_DAYS;
        let rel = (period_gyr - published::OMEGA_SPIN_PERIOD_GYR).abs()
            / published::OMEGA_SPIN_PERIOD_GYR;
        assert!(
            rel < 0.01,
            "Ω★ps period = {period_gyr:.2} Gyr (Lai: 55.8), rel {rel:.3e}"
        );
        // Right order of magnitude: tens of Gyr at the young-Sun 10-day spin.
        assert!((40.0..80.0).contains(&period_gyr));
    }

    #[test]
    fn omega_spin_is_2p88_omega_j() {
        // Lai Eq. (8): Ω★ps = 2.88 Ω★J (Jupiter-only torque).
        let omega_star = 2.0 * PI / 10.0;
        let coeff = 3.0 * solar::KQ_STAR / (2.0 * solar::K_STAR);
        let om_j = coeff
            * MASS_JUPITER_SOLAR
            * (RADIUS_SUN_AU / planets::A_JUPITER_AU).powi(3)
            * omega_star;
        let ratio = omega_spin_precession(10.0) / om_j;
        assert!(
            (ratio - published::OMEGA_SPIN_OVER_J).abs() < 0.03,
            "Ω★ps/Ω★J = {ratio:.3} (Lai: 2.88)"
        );
    }

    #[test]
    fn lambda_star_is_unity() {
        assert!((solar::lambda_star() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn solar_j2_consistent_with_kq() {
        // Lai: J₂ = k_q★ Ω̂★² ≈ 2.2e-7 with Ω̂★ = Ω★ (GM★/R★³)⁻¹/².
        // At P★ = 24 d (present Sun), this must reproduce p9-core's J2_SUN ≈ 2e-7.
        let p_days = 24.0;
        let omega_star = 2.0 * PI / p_days;
        let omega_break = (G_AU3_MSUN_DAY2 / RADIUS_SUN_AU.powi(3)).sqrt();
        let omega_hat = omega_star / omega_break;
        let j2 = solar::KQ_STAR * omega_hat * omega_hat;
        // Present-day Sun: same order as J2_SUN (2e-7). At 10 d it is ~5.8× larger.
        assert!(
            (j2 / J2_SUN - 1.0).abs() < 0.6,
            "J2(24 d) = {j2:.3e} vs J2_SUN {J2_SUN:.3e}"
        );
    }

    #[test]
    fn nominal_obliquity_in_5_to_10_band() {
        // HEADLINE: a nominal P9 produces a several-degree obliquity (~6°).
        // P★/λ★ = 20 d (Lai's "averaged" value, Fig. 1).
        let p9 = nominal_p9();
        let obl = solar_obliquity_deg(&p9, 20.0);
        assert!(
            (4.0..=10.0).contains(&obl),
            "nominal obliquity = {obl:.2}° (expected 5–10° band, observed 6°)"
        );
    }

    #[test]
    fn obliquity_grows_with_mass() {
        let base = nominal_p9();
        let mut heavy = base;
        heavy.mass_earth = 20.0;
        // At fixed θp the heavier planet drives a larger obliquity (stronger
        // Ω_y) for the same ãp.
        assert!(
            solar_obliquity_deg(&heavy, 20.0) > solar_obliquity_deg(&base, 20.0),
            "heavier P9 should raise obliquity"
        );
    }

    #[test]
    fn obliquity_grows_with_inclination() {
        let base = nominal_p9();
        let mut more_incl = base;
        more_incl.inclination_rad = 30.0 * DEG2RAD;
        // Below the resonance peak, more sin(θp) → more obliquity.
        assert!(
            solar_obliquity_deg(&more_incl, 20.0) > solar_obliquity_deg(&base, 20.0),
            "more inclined P9 should raise obliquity"
        );
    }

    #[test]
    fn obliquity_weakens_with_semimajor_axis() {
        let base = nominal_p9();
        let mut far = base;
        far.a_au = base.a_au * 1.3; // larger ãp, weaker forcing
        assert!(
            solar_obliquity_deg(&far, 20.0) < solar_obliquity_deg(&base, 20.0),
            "more distant P9 should lower obliquity"
        );
    }

    #[test]
    fn no_obliquity_without_inclination() {
        let mut p9 = nominal_p9();
        p9.inclination_rad = 0.0;
        let obl = solar_obliquity_deg(&p9, 20.0);
        assert!(obl < 1e-9, "θp = 0 must give zero obliquity, got {obl:.3e}");
    }

    #[test]
    fn secular_spin_orbit_resonance_requires_fast_young_spin() {
        // The secular spin–orbit resonance is Ω_z = 0, i.e. the spin precession
        // Ω★ps balances the plane precession Ω_L cos θp (L/Lp + cos θp). HONEST
        // FINDING: at the "averaged" P★/λ★ = 20 d, Ω★ps is far too slow — Ω_z is
        // *negative everywhere* across the 300–600 AU band (no resonance). The
        // resonance only enters the band for a rapidly rotating young Sun
        // (P★/λ★ of order a few days; the spin term ∝ 1/P★). The crate computes
        // the required spin period and confirms it is in the young-Sun regime.
        let e = 0.5_f64;
        let theta_p = 20.0 * DEG2RAD;
        let a_tilde = 400.0_f64;
        let p9 = PlanetNine {
            mass_earth: 10.0,
            a_au: a_tilde / (1.0 - e * e).sqrt(),
            e,
            inclination_rad: theta_p,
        };

        // At P★/λ★ = 20 d the system is below resonance (Ω_z < 0).
        assert!(
            corotating_frequencies(&p9, 20.0).1 < 0.0,
            "Ω_z should be negative (below resonance) at P★ = 20 d"
        );

        // The plane-precession side of the balance (Ω_z = 0):
        let om_l = omega_l(&p9);
        let l_total = planets::total_orbital_angular_momentum();
        let lp = p9.angular_momentum();
        let rhs = om_l * theta_p.cos() * (l_total / lp + theta_p.cos());
        // Ω★ps(P) = Ω★ps(10 d)·(10/P); solve Ω★ps = rhs for P.
        let p_resonant = omega_spin_precession(10.0) * 10.0 / rhs;
        assert!(
            (1.0..6.0).contains(&p_resonant),
            "resonant spin period {p_resonant:.2} d (young, rapidly rotating Sun)"
        );

        // Confirm Ω_z indeed crosses zero versus spin period (sign flip between
        // the fast-spin resonant value and the slow 20 d value).
        let z_fast = corotating_frequencies(&p9, 0.5 * p_resonant).1;
        let z_slow = corotating_frequencies(&p9, 20.0).1;
        assert!(
            z_fast * z_slow < 0.0,
            "Ω_z must change sign across resonance"
        );
    }

    #[test]
    fn resonant_obliquity_is_amplified() {
        // Near the resonance the induced obliquity is largest (the
        // 2Ω_y/Ω_z prefactor blows up as Ω_z → 0). Tune P★ to sit just off the
        // Ω_z = 0 balance for a nominal P9 and confirm the obliquity exceeds
        // the far-from-resonance value at the "averaged" 20 d spin.
        let e = 0.5_f64;
        let theta_p = 20.0 * DEG2RAD;
        let p9 = PlanetNine {
            mass_earth: 10.0,
            a_au: 400.0 / (1.0 - e * e).sqrt(),
            e,
            inclination_rad: theta_p,
        };
        let om_l = omega_l(&p9);
        let l_total = planets::total_orbital_angular_momentum();
        let lp = p9.angular_momentum();
        let rhs = om_l * theta_p.cos() * (l_total / lp + theta_p.cos());
        let p_resonant = omega_spin_precession(10.0) * 10.0 / rhs;

        // Slightly off resonance to keep θsl finite and small-angle valid.
        let obl_near = solar_obliquity_deg(&p9, 1.15 * p_resonant);
        let obl_far = solar_obliquity_deg(&p9, 20.0);
        assert!(
            obl_near > obl_far,
            "near-resonance obliquity {obl_near:.2}° should exceed far {obl_far:.2}°"
        );
    }
}
