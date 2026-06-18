//! Spin-axis precession constant α for Uranus (Lu & Laughlin 2022, Eq. 2–3).
//!
//! A planet's spin axis precesses about its orbit normal at the rate
//!
//!   α = (3 n² / 2 ω) · (J₂ + q) / (λ + l)
//!
//! where
//!   n  = heliocentric mean motion,
//!   ω  = planetary spin frequency,
//!   J₂ = gravitational quadrupole of the planet,
//!   λ  = C/(Mₚ Rₚ²), the normalized polar moment of inertia,
//!   q  = ½ Σᵢ (Mᵢ/Mₚ)(aᵢ/Rₚ)²        (satellite oblateness contribution),
//!   l  = (1/Mₚ Rₚ² ω) Σᵢ Mᵢ aᵢ² nᵢ   (satellite spin-angular-momentum term).
//!
//! For Uranus the regular satellites dominate q and l: the moons contribute
//! roughly five times the bare J₂ to the effective quadrupole, which is why
//! Uranus' present precession is so slow (Tα ≈ 1.7×10² Myr). The paper quotes
//! α ≈ 0.045 arcsec/yr; this module computes it from first principles.

use std::f64::consts::PI;

use p9_core::constants::{J2_URANUS, RAD2DEG, RADIUS_URANUS_AU};
use p9_core::units::{arcseconds, julian_year, km, AngularVelocity, Length};

/// A regular satellite: (mass in kg, semi-major axis in km).
#[derive(Debug, Clone, Copy)]
pub struct Satellite {
    pub mass_kg: f64,
    pub a_km: f64,
}

impl Satellite {
    /// Orbital semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        km(self.a_km)
    }
}

/// Uranus' polar mass in kg (JPL).
pub const URANUS_MASS_KG: f64 = 8.6810e25;

/// Uranus equatorial radius in km (matches `RADIUS_URANUS_AU`).
pub const URANUS_RADIUS_KM: f64 = 25_559.0;

/// Uranus spin period in hours (sidereal, 17.24 h).
pub const URANUS_SPIN_PERIOD_HR: f64 = 17.24;

/// Normalized polar moment of inertia C/(Mₚ Rₚ²) for Uranus.
/// Lu & Laughlin adopt λ = 0.225.
pub const URANUS_LAMBDA: f64 = 0.225;

/// Uranus orbital semi-major axis (AU).
pub const URANUS_A_AU: f64 = 19.2184;

/// The five major regular satellites of Uranus (mass kg, a km). These set the
/// q and l terms; the minor inner moons are negligible for both sums.
pub const URANUS_MOONS: [Satellite; 5] = [
    Satellite {
        mass_kg: 6.59e19,
        a_km: 129_900.0,
    }, // Miranda
    Satellite {
        mass_kg: 1.353e21,
        a_km: 190_900.0,
    }, // Ariel
    Satellite {
        mass_kg: 1.172e21,
        a_km: 266_000.0,
    }, // Umbriel
    Satellite {
        mass_kg: 3.527e21,
        a_km: 436_300.0,
    }, // Titania
    Satellite {
        mass_kg: 3.014e21,
        a_km: 583_500.0,
    }, // Oberon
];

/// Heliocentric mean motion n (rad/yr) from Kepler's third law: P = a^{3/2} yr.
pub fn mean_motion_rad_per_yr(a_au: f64) -> f64 {
    2.0 * PI / a_au.powf(1.5)
}

/// Planetary spin frequency ω (rad/yr) from the spin period in hours.
pub fn spin_frequency_rad_per_yr(period_hr: f64) -> f64 {
    let period_yr = period_hr / 24.0 / 365.25;
    2.0 * PI / period_yr
}

/// Satellite oblateness contribution q = ½ Σ (Mᵢ/Mₚ)(aᵢ/Rₚ)² (Eq. 3).
pub fn satellite_q(moons: &[Satellite], m_planet_kg: f64, r_planet_km: f64) -> f64 {
    0.5 * moons
        .iter()
        .map(|s| (s.mass_kg / m_planet_kg) * (s.a_km / r_planet_km).powi(2))
        .sum::<f64>()
}

/// Satellite spin-angular-momentum contribution
/// l = (1 / Mₚ Rₚ² ω) Σ Mᵢ aᵢ² nᵢ  (Eq. 3), made dimensionless by using the
/// same length unit (km) and frequency unit (rad/yr) throughout. nᵢ is each
/// moon's orbital mean motion about the planet, GMₚ = n² aᵢ³.
pub fn satellite_l(
    moons: &[Satellite],
    m_planet_kg: f64,
    r_planet_km: f64,
    omega_rad_per_yr: f64,
) -> f64 {
    // GM of the planet in km^3 / yr^2, from the outermost moon as a calibrator
    // is unnecessary — we get nᵢ directly from each moon's own GM via Kepler.
    // GMₚ (km^3/yr^2): use G·Mₚ with G in SI then convert.
    let g_si = 6.674_30e-11; // m^3 kg^-1 s^-2
    let sec_per_yr = 365.25 * 86_400.0;
    // GMₚ in km^3/yr^2: (m^3/s^2) -> (km^3/yr^2): /1e9 * sec_per_yr^2
    let gm_planet = g_si * m_planet_kg / 1.0e9 * sec_per_yr * sec_per_yr;
    let denom = m_planet_kg * r_planet_km * r_planet_km * omega_rad_per_yr;
    moons
        .iter()
        .map(|s| {
            let n_i = (gm_planet / s.a_km.powi(3)).sqrt(); // rad/yr
            s.mass_kg * s.a_km * s.a_km * n_i
        })
        .sum::<f64>()
        / denom
}

/// The full spin-axis precession constant α (rad/yr) for Uranus today
/// (Eq. 2). `cos θ` is folded in: α as defined is the precession-constant
/// prefactor; the precession rate is α cos θ. Here we return the bare α.
pub fn alpha_rad_per_yr() -> f64 {
    let n = mean_motion_rad_per_yr(URANUS_A_AU);
    let omega = spin_frequency_rad_per_yr(URANUS_SPIN_PERIOD_HR);
    let q = satellite_q(&URANUS_MOONS, URANUS_MASS_KG, URANUS_RADIUS_KM);
    let l = satellite_l(&URANUS_MOONS, URANUS_MASS_KG, URANUS_RADIUS_KM, omega);
    (3.0 * n * n / (2.0 * omega)) * (J2_URANUS + q) / (URANUS_LAMBDA + l)
}

/// α in arcsec/yr (the paper's units).
pub fn alpha_arcsec_per_yr() -> f64 {
    alpha_rad_per_yr() * RAD2DEG * 3600.0
}

/// α as a typed [`AngularVelocity`] (arcsec per Julian year).
pub fn alpha_typed() -> AngularVelocity {
    (arcseconds(alpha_arcsec_per_yr()) / julian_year()).into()
}

/// Spin-axis precession period Tα = 2π / (α cos θ) in Myr at obliquity θ.
pub fn precession_period_myr(theta_rad: f64) -> f64 {
    let rate = alpha_rad_per_yr() * theta_rad.cos().abs();
    (2.0 * PI / rate) / 1.0e6
}

/// Convenience: a sanity check that `RADIUS_URANUS_AU` and the km radius agree.
pub fn radius_consistency_ratio() -> f64 {
    let au_km = 1.495_978_707e8;
    (RADIUS_URANUS_AU * au_km) / URANUS_RADIUS_KM
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn typed_accessors_match_f64() {
        use uom::si::angular_velocity::radian_per_second;
        use uom::si::length::kilometer;
        let moon = URANUS_MOONS[0];
        assert_relative_eq!(
            moon.semi_major_axis().get::<kilometer>(),
            moon.a_km,
            epsilon = 1e-9
        );
        // α typed (arcsec/yr) read back in rad/s matches the f64 source.
        let alpha_rad_s = alpha_typed().get::<radian_per_second>();
        let expect = alpha_rad_per_yr() / (365.25 * 86_400.0);
        assert_relative_eq!(alpha_rad_s, expect, max_relative = 1e-9);
    }

    #[test]
    fn radius_km_matches_core_radius_au() {
        // Sanity: our km radius reproduces p9-core's RADIUS_URANUS_AU.
        assert_relative_eq!(radius_consistency_ratio(), 1.0, epsilon = 1e-6);
    }

    #[test]
    fn satellite_q_dominates_bare_j2() {
        // Lu & Laughlin note the moons dominate the effective quadrupole.
        let q = satellite_q(&URANUS_MOONS, URANUS_MASS_KG, URANUS_RADIUS_KM);
        assert!(q > 0.0);
        assert!(
            q > 3.0 * J2_URANUS,
            "satellite q = {q:.4e} should exceed bare J2 = {J2_URANUS:.4e} severalfold"
        );
        // Roughly q ~ 0.016 from the Uranian system.
        assert!(
            (0.008..0.030).contains(&q),
            "q = {q:.4e} out of expected band"
        );
    }

    #[test]
    fn mean_motion_and_spin_orders_of_magnitude() {
        let n = mean_motion_rad_per_yr(URANUS_A_AU);
        let omega = spin_frequency_rad_per_yr(URANUS_SPIN_PERIOD_HR);
        // Uranus orbital period ~84 yr -> n ~ 0.075 rad/yr.
        assert_relative_eq!(n, 2.0 * PI / 84.0, epsilon = 0.01);
        // Spin period 17.24 h -> ω huge compared to n.
        assert!(omega / n > 1e4, "ω/n = {:.3e}", omega / n);
    }

    #[test]
    fn alpha_matches_paper_present_value() {
        // HEADLINE PIN: Lu & Laughlin (2022) quote α ≈ 0.045 arcsec/yr today.
        let alpha = alpha_arcsec_per_yr();
        assert!(
            (0.030..0.070).contains(&alpha),
            "α = {alpha:.4} arcsec/yr; paper quotes ≈0.045"
        );
        // Tight pin: within 20% of 0.045.
        assert_relative_eq!(alpha, 0.045, max_relative = 0.20);
    }

    #[test]
    fn precession_period_right_order_of_magnitude() {
        // Paper: Tα ≈ 169 Myr at the present 98° obliquity. cos(98°) is small,
        // so the precession about the orbit normal is slow; we land within a
        // factor of a few — see crate-level note on residuals.
        let t = precession_period_myr(crate::uranus_obliquity_rad());
        assert!(t > 50.0 && t < 2000.0, "Tα = {t:.1} Myr (paper ~169)");
    }
}
