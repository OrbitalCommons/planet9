//! Linear-theory growth rate of the inclination instability.
//!
//! Madigan & McCourt (2016) show that a self-gravitating disc of eccentric
//! planetesimals develops an exponentially growing inclination mode on a
//! *secular* timescale. They measure the instability against the disc's
//! secular time, the orbital period scaled by the Sun-to-disc mass ratio
//! (their Section 2):
//!
//!   t_sec = (M_sun / M_d) * P,    P = 2 pi sqrt(a^3 / GM_sun).
//!
//! This is the natural clock of self-gravity-driven secular precession in a
//! near-Keplerian disc. The disc secular precession *frequency* is the angular
//! analogue,
//!
//!   omega_sec = (M_d / M_sun) * n,    n = sqrt(GM_sun / a^3),
//!
//! (standard secular disc theory, e.g. Murray & Dermott 1999, Ch. 7); this is
//! the precession rate that competes with the giant planets in the Das &
//! Batygin (2023) suppression criterion (see [`crate::suppression`]).
//!
//! The inclination e-folds in a small order-unity number of disc secular
//! times, so the e-folding (growth) rate is
//!
//!   gamma = (M_d/M_sun) * n / C_growth,   tau = 1/gamma = (C_growth/2 pi) t_sec,
//!
//! where `C_growth` is an eigenvalue of the linearised secular operator. All
//! masses are in solar masses, `a` in AU; `n` and `P` come from `p9_core`, so
//! `gamma` is in 1/day. E-folding times are converted to Myr/Gyr for reporting.

use p9_core::constants::{GM_SUN, GYR_DAYS, TWO_PI, YEAR_DAYS};
use p9_core::types::OrbitalElements;

/// Dimensionless secular-eigenvalue prefactor in gamma = (M_d/M_sun) n /
/// C_GROWTH, equivalently tau = (C_GROWTH/2pi) * t_sec with t_sec = (M_sun/M_d)P.
///
/// Madigan & McCourt (2016) find the inclination e-folds on the order of one
/// disc secular time t_sec (their Fig. 3 / Section 4). Choosing C_GROWTH = 2 pi
/// sets one e-fold per disc secular *period* (tau = 1 * t_sec), the published
/// order. This is a *labelled* reference value tied to the paper's reported
/// secular-time scaling, not a free fit; the resulting fiducial timescale is
/// pinned in the tests and documented in the crate's residual notes.
pub const C_GROWTH: f64 = TWO_PI;

/// Mean motion n = sqrt(GM_sun / a^3) in rad/day for an orbit of semi-major
/// axis `a` (AU), reusing `p9_core`'s element machinery so the Kepler
/// convention is shared with the rest of the workspace.
pub fn mean_motion(a: f64) -> f64 {
    let elem = OrbitalElements {
        a,
        e: 0.0,
        i: 0.0,
        omega_big: 0.0,
        omega: 0.0,
        mean_anomaly: 0.0,
    };
    elem.mean_motion(GM_SUN)
}

/// Orbital period P = 2 pi sqrt(a^3 / GM_sun) in days.
pub fn orbital_period(a: f64) -> f64 {
    TWO_PI / mean_motion(a)
}

/// Disc self-gravity secular precession frequency omega_sec = (M_d/M_sun) * n
/// [rad/day]. This is the rate compared against the giant planets in the
/// suppression criterion. `m_disk_solar`: disc mass in solar masses;
/// `a`: characteristic semi-major axis (AU).
pub fn secular_frequency(m_disk_solar: f64, a: f64) -> f64 {
    m_disk_solar * mean_motion(a)
}

/// Disc secular time t_sec = (M_sun/M_d) * P [days].
pub fn secular_time_days(m_disk_solar: f64, a: f64) -> f64 {
    orbital_period(a) / m_disk_solar
}

/// Linear-theory growth rate gamma = (M_d/M_sun) * n / C_GROWTH  [1/day].
/// Monotonically increasing in `m_disk_solar`.
pub fn growth_rate(m_disk_solar: f64, a: f64) -> f64 {
    secular_frequency(m_disk_solar, a) / C_GROWTH
}

/// E-folding time tau = 1/gamma  [days]. Decreases with disc mass.
pub fn efolding_time_days(m_disk_solar: f64, a: f64) -> f64 {
    1.0 / growth_rate(m_disk_solar, a)
}

/// E-folding time in millions of years.
pub fn efolding_time_myr(m_disk_solar: f64, a: f64) -> f64 {
    efolding_time_days(m_disk_solar, a) / (YEAR_DAYS * 1.0e6)
}

/// E-folding time in billions of years (Gyr).
pub fn efolding_time_gyr(m_disk_solar: f64, a: f64) -> f64 {
    efolding_time_days(m_disk_solar, a) / GYR_DAYS
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::EARTH_MASS_SOLAR;

    /// Madigan & McCourt's fiducial disc: ~tens of Earth masses spread over
    /// a few hundred AU. We use 20 M_earth at a = 250 AU as a fiducial massive
    /// disc (within their explored 1-100 M_earth, 100s of AU range).
    const FIDUCIAL_MASS_EARTH: f64 = 20.0;
    const FIDUCIAL_A_AU: f64 = 250.0;

    fn fiducial_mass_solar() -> f64 {
        FIDUCIAL_MASS_EARTH * EARTH_MASS_SOLAR
    }

    #[test]
    fn test_growth_rate_monotonic_in_mass() {
        let a = FIDUCIAL_A_AU;
        let masses = [1.0, 5.0, 10.0, 20.0, 50.0];
        let mut prev = 0.0;
        for &m_e in &masses {
            let g = growth_rate(m_e * EARTH_MASS_SOLAR, a);
            assert!(g > prev, "growth rate not increasing at {m_e} M_earth");
            prev = g;
        }
    }

    #[test]
    fn test_efolding_time_decreases_with_mass() {
        let a = FIDUCIAL_A_AU;
        let tau_light = efolding_time_myr(5.0 * EARTH_MASS_SOLAR, a);
        let tau_heavy = efolding_time_myr(50.0 * EARTH_MASS_SOLAR, a);
        assert!(
            tau_heavy < tau_light,
            "heavier disc should e-fold faster: {tau_heavy} vs {tau_light}"
        );
        // Inverse-proportional in mass: 10x mass -> 1/10 the time.
        let ratio = tau_light / tau_heavy;
        assert!(
            (ratio - 10.0).abs() < 1e-9,
            "expected exact 10x scaling, got {ratio}"
        );
    }

    #[test]
    fn test_fiducial_timescale_published_order() {
        // The instability e-folds in ~one disc secular time t_sec =
        // (M_sun/M_d) P, so the absolute timescale runs from tens of Myr (the
        // heaviest ~tens-M_earth discs) up to ~Gyr (the lightest few-M_earth
        // discs Das & Batygin flag as marginal). For the fiducial 20 M_earth
        // disc at 250 AU we compute tau ~ 66 Myr (one secular time). A 5
        // M_earth disc gives ~260 Myr; the band spans the published "secular
        // (hundreds of Myr to Gyr)" descriptor across the disc-mass range.
        // RESIDUAL: a single fiducial mass lands at the low end (66 Myr); the
        // hundreds-of-Myr figure corresponds to lighter discs.
        let tau_fid = efolding_time_myr(fiducial_mass_solar(), FIDUCIAL_A_AU);
        assert!(
            (20.0..=2000.0).contains(&tau_fid),
            "fiducial e-folding time {tau_fid} Myr outside published order"
        );
        // A few-Earth-mass disc reaches the hundreds-of-Myr regime.
        let tau_light = efolding_time_myr(5.0 * EARTH_MASS_SOLAR, FIDUCIAL_A_AU);
        assert!(
            (100.0..=5000.0).contains(&tau_light),
            "5 M_earth e-folding time {tau_light} Myr outside published order"
        );
    }

    #[test]
    fn test_efolding_time_is_one_secular_time() {
        // With C_GROWTH = 2 pi the e-folding time equals one disc secular
        // time t_sec = (M_sun/M_d) P exactly.
        let m = fiducial_mass_solar();
        let a = FIDUCIAL_A_AU;
        let tau = efolding_time_days(m, a);
        let t_sec = secular_time_days(m, a);
        assert!((tau / t_sec - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_secular_frequency_scales_with_orbital_time() {
        // omega_sec = (M_d/M_sun) n: doubling the period (a -> a*2^(2/3))
        // halves the frequency at fixed mass.
        let m = fiducial_mass_solar();
        let a1 = 250.0;
        let a2 = a1 * 2.0_f64.powf(2.0 / 3.0); // period doubles
        let w1 = secular_frequency(m, a1);
        let w2 = secular_frequency(m, a2);
        assert!((w1 / w2 - 2.0).abs() < 1e-9);
    }
}
