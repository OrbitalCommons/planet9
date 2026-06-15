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

/// Skumanich spin-down: ω(t) = ω₀ * (t₀/t)^{1/2}, capped at the initial
/// rotation rate.
///
/// We normalize so that ω(t_age) = present-day ω, and cap the divergent
/// t → 0 behavior at the P = 10 d initial rotation:
///
///   ω(t) = min( ω(P = 10 d), ω_present √(t_age/t) ).
///
/// Without the cap the first integration step sees ω ≈ 300× present (beyond
/// breakup), and since the spin-orbit coupling is largest early, the whole
/// obliquity excitation is inflated.
///
/// `t` is time since formation in days.
/// `t_age` is the total age (4.5 Gyr in days).
pub fn solar_omega_at_time(t: f64, t_age: f64) -> f64 {
    let omega_initial = omega_sun(P_SUN_INITIAL_DAYS);
    if t <= 0.0 {
        return omega_initial;
    }
    let omega_present = omega_sun(P_SUN_PRESENT_DAYS);
    (omega_present * (t_age / t).sqrt()).min(omega_initial)
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
    fn test_spindown_capped_at_initial_rotation() {
        let t_age = 4.5 * GYR_DAYS;
        let omega_initial = omega_sun(P_SUN_INITIAL_DAYS);

        // Early times: the divergent √(t_age/t) is capped at ω(P = 10 d).
        for &t in &[0.0, 1.0, 1e3, 1e6, 0.01 * t_age] {
            let w = solar_omega_at_time(t, t_age);
            assert!(
                w <= omega_initial * (1.0 + 1e-12),
                "ω(t={t:.1e}) = {w:.3e} exceeds the initial rotation"
            );
        }

        // The cap releases exactly at t_c = t_age (ω_present/ω_initial)².
        let t_c = t_age * (P_SUN_INITIAL_DAYS / P_SUN_PRESENT_DAYS).powi(2);
        let w_after = solar_omega_at_time(1.01 * t_c, t_age);
        assert!(w_after < omega_initial);

        // Monotonically non-increasing
        let mut prev = f64::INFINITY;
        for k in 0..100 {
            let w = solar_omega_at_time(k as f64 * t_age / 100.0, t_age);
            assert!(w <= prev + 1e-15);
            prev = w;
        }
    }

    #[test]
    fn test_solar_ring_semimajor_axis() {
        let omega = omega_sun(P_SUN_PRESENT_DAYS);
        let a_tilde = solar_ring_semimajor_axis(omega);

        // Should be very small (inside the Sun)
        assert!(a_tilde > 0.0 && a_tilde < R_SUN_AU);
    }
}
