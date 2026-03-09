//! Solar interior model and spin-down from Bailey+ (2016).
//!
//! The Sun is modeled as an n=3 polytrope with:
//!   - Dimensionless moment of inertia I_hat = 0.08
//!   - Apsidal motion constant k₂ = 0.01 (Love number)
//!   - Initial rotation period ~10 days
//!   - Skumanich spin-down: ω ∝ t^{-1/2}

use std::f64::consts::PI;

use p9_core::constants::*;

/// Dimensionless moment of inertia for n=3 polytrope (Table 1 of Bailey+2016).
pub const I_HAT: f64 = 0.08;

/// Apsidal motion constant (Love number k₂) for the Sun.
pub const K2_SUN: f64 = 0.01;

/// Solar radius in AU.
pub const R_SUN_AU: f64 = RADIUS_SUN_AU;

/// Solar mass in solar masses (trivially 1.0).
pub const M_SUN_SOLAR: f64 = 1.0;

/// Present-day solar rotation period in days (~25.4 day sidereal).
pub const P_SUN_PRESENT_DAYS: f64 = 25.4;

/// Initial solar rotation period in days (~10 days for young Sun).
pub const P_SUN_INITIAL_DAYS: f64 = 10.0;

/// Angular velocity of the Sun from rotation period (rad/day).
pub fn omega_sun(period_days: f64) -> f64 {
    2.0 * PI / period_days
}

/// Skumanich spin-down: ω(t) = ω₀ * (t₀/t)^{1/2}
///
/// We normalize so that ω(4.5 Gyr) = present-day ω.
/// Then ω(t) = ω_present * sqrt(t_age / t) for t > 0.
///
/// `t` is time since formation in days.
/// `t_age` is the total age (4.5 Gyr in days).
pub fn solar_omega_at_time(t: f64, t_age: f64) -> f64 {
    let omega_present = omega_sun(P_SUN_PRESENT_DAYS);
    if t <= 0.0 {
        return omega_sun(P_SUN_INITIAL_DAYS);
    }
    omega_present * (t_age / t).sqrt()
}

/// Spin angular momentum of the Sun: L = I_hat * M * R² * ω
///
/// Returns value in units of M_sun * AU² / day.
pub fn solar_spin_angular_momentum(omega: f64) -> f64 {
    I_HAT * M_SUN_SOLAR * R_SUN_AU * R_SUN_AU * omega
}

/// Effective semi-major axis of the "solar spin ring" (Eq. 9 of Bailey+2016):
///
///   a_tilde = [ 16 * ω² * k₂² * R⁶ / (9 * I_hat² * G * M_sun) ]^{1/3}
///
/// This maps the rotational bulge's gravitational effect to an equivalent wire.
/// Returns value in AU.
pub fn solar_ring_semimajor_axis(omega: f64) -> f64 {
    let r6 = R_SUN_AU.powi(6);
    let numerator = 16.0 * omega * omega * K2_SUN * K2_SUN * r6;
    let denominator = 9.0 * I_HAT * I_HAT * G_AU3_MSUN_DAY2 * M_SUN_SOLAR;
    (numerator / denominator).powf(1.0 / 3.0)
}

/// Combined angular momentum of the four giant planets.
///
/// L_planets = Σ m_i * sqrt(G * M_sun * a_i) for each giant planet.
/// Returns in M_sun * AU² / day.
pub fn giant_planet_angular_momentum() -> f64 {
    let planets = [
        (MASS_JUPITER_SOLAR, 5.2026),  // Jupiter a in AU
        (MASS_SATURN_SOLAR, 9.5549),   // Saturn
        (MASS_URANUS_SOLAR, 19.2184),  // Uranus
        (MASS_NEPTUNE_SOLAR, 30.1104), // Neptune
    ];

    planets
        .iter()
        .map(|&(m, a)| m * (G_AU3_MSUN_DAY2 * M_SUN_SOLAR * a).sqrt())
        .sum()
}

/// Effective semi-major axis for the combined giant planet ring.
///
/// The combined ring has mass m_eff and semi-major axis a_eff such that
/// L = m_eff * sqrt(G * M_sun * a_eff).
///
/// We compute a_eff as the angular-momentum-weighted mean.
pub fn giant_planet_effective_orbit() -> (f64, f64) {
    let planets = [
        (MASS_JUPITER_SOLAR, 5.2026),
        (MASS_SATURN_SOLAR, 9.5549),
        (MASS_URANUS_SOLAR, 19.2184),
        (MASS_NEPTUNE_SOLAR, 30.1104),
    ];

    let m_total: f64 = planets.iter().map(|&(m, _)| m).sum();
    let l_total = giant_planet_angular_momentum();

    // a_eff from L = m_total * sqrt(GM * a_eff)
    let a_eff = (l_total / m_total).powi(2) / (G_AU3_MSUN_DAY2 * M_SUN_SOLAR);

    (m_total, a_eff)
}

/// Planet Nine angular momentum magnitude.
///
/// L₉ = m₉ * sqrt(G * M_sun * a₉) * sqrt(1 - e₉²)
pub fn p9_angular_momentum(mass_solar: f64, a: f64, e: f64) -> f64 {
    mass_solar * (G_AU3_MSUN_DAY2 * M_SUN_SOLAR * a).sqrt() * (1.0 - e * e).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_skumanich_spindown() {
        let t_age = 4.5 * GYR_DAYS;
        let omega_now = solar_omega_at_time(t_age, t_age);
        let omega_present = omega_sun(P_SUN_PRESENT_DAYS);

        // At present age, should match present-day omega
        assert!(
            (omega_now - omega_present).abs() / omega_present < 1e-10,
            "ω at present: {:.6e} vs expected {:.6e}",
            omega_now,
            omega_present
        );

        // At half the age, should be faster by sqrt(2)
        let omega_half = solar_omega_at_time(t_age / 2.0, t_age);
        let ratio = omega_half / omega_now;
        assert!(
            (ratio - std::f64::consts::SQRT_2).abs() < 1e-10,
            "ratio = {:.6}, expected sqrt(2)",
            ratio
        );
    }

    #[test]
    fn test_solar_ring_semimajor_axis() {
        let omega = omega_sun(P_SUN_PRESENT_DAYS);
        let a_tilde = solar_ring_semimajor_axis(omega);

        // Should be very small (inside the Sun)
        assert!(a_tilde > 0.0 && a_tilde < R_SUN_AU);
    }

    #[test]
    fn test_giant_planet_angular_momentum() {
        let l = giant_planet_angular_momentum();
        // Should be positive and dominated by Jupiter
        assert!(l > 0.0);

        let l_jup = MASS_JUPITER_SOLAR * (G_AU3_MSUN_DAY2 * 5.2026_f64).sqrt();
        // Jupiter should be > 50% of total
        assert!(l_jup / l > 0.5);
    }

    #[test]
    fn test_p9_angular_momentum() {
        let m9 = 10.0 * EARTH_MASS_SOLAR;
        let l9 = p9_angular_momentum(m9, 700.0, 0.6);
        let l_planets = giant_planet_angular_momentum();

        // P9's angular momentum is significant but smaller than giant planets
        // At a=700 AU, L_9/L_gp ~ 0.2 (large a compensates small mass)
        assert!(l9 < l_planets);
        assert!(l9 > 0.0);
    }
}
