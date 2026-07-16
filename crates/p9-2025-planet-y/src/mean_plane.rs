//! Observable mean-plane warp and a simple debiasing of an inclination-limited
//! sample.
//!
//! Siraj et al. (2025) infer the local Laplace plane from the *measured* mean
//! plane of Kuiper-belt objects binned in semi-major axis. The forced
//! inclination computed in [`crate::laplace_plane`] is the underlying signal;
//! what a survey records is the mean pole of objects that pass an
//! inclination-dependent detectability cut. A flat, isotropic-in-node belt
//! whose particles librate about the forced plane has a mean pole *at* the
//! forced plane — but detection efficiency that favours low ecliptic
//! inclination biases the estimate toward the ecliptic. This module computes
//! the forced mean-plane pole and applies a first-order debiasing correction.
//!
//! The "mean plane" here is expressed as a pole vector via p9-core's pole
//! geometry conventions, reusing the same (sin i sin Ω, −sin i cos Ω, cos i)
//! normal used in `p9-2025-clustering::orbital_poles`.

use nalgebra::Vector3;
use p9_core::constants::DEG2RAD;
use p9_core::coords::sky::angular_distance;

use crate::laplace_plane::{forced_inclination, PlanetY};

/// Pole unit vector of the forced (local Laplace) plane at semi-major axis `a`.
///
/// The forced plane is tilted by `i_forced` toward Planet Y's plane. Taking
/// Planet Y's node as the reference longitude (Ω = 0), the forced pole is
/// n̂ = (sin i_f sin Ω, −sin i_f cos Ω, cos i_f) with Ω = 0, i.e.
/// n̂ = (0, −sin i_f, cos i_f).
pub fn forced_pole(a_au: f64, planet_y: &PlanetY) -> Vector3<f64> {
    let i_f = forced_inclination(a_au, planet_y);
    Vector3::new(0.0, -i_f.sin(), i_f.cos())
}

/// Reference (inner invariable) pole = ecliptic z-axis in this reduced frame.
pub fn reference_pole() -> Vector3<f64> {
    Vector3::new(0.0, 0.0, 1.0)
}

/// Tilt (deg) of the local mean plane from the reference plane at `a`,
/// measured as the true spherical angle between poles (reusing p9-core's
/// `angular_distance`, the same metric the clustering crate uses for poles).
pub fn mean_plane_tilt_deg(a_au: f64, planet_y: &PlanetY) -> f64 {
    angular_distance(&forced_pole(a_au, planet_y), &reference_pole()) / DEG2RAD
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::laplace_plane::published;
    use approx::assert_relative_eq;

    fn nominal() -> PlanetY {
        PlanetY::new(0.3, 150.0, published::I_Y_DEG)
    }

    #[test]
    fn test_forced_pole_is_unit() {
        let p = forced_pole(70.0, &nominal());
        assert_relative_eq!(p.norm(), 1.0, epsilon = 1e-12);
    }

    #[test]
    fn test_tilt_matches_forced_inclination() {
        // With Ω = 0 the pole tilt equals i_forced exactly.
        let py = nominal();
        let tilt = mean_plane_tilt_deg(72.0, &py);
        let i_f = forced_inclination(72.0, &py) / DEG2RAD;
        assert_relative_eq!(tilt, i_f, max_relative = 1e-9);
    }

    #[test]
    fn test_tilt_rises_across_band() {
        let py = nominal();
        let t50 = mean_plane_tilt_deg(50.0, &py);
        let t80 = mean_plane_tilt_deg(80.0, &py);
        assert!(
            t80 > t50,
            "mean-plane tilt should rise: {t50:.3} → {t80:.3}"
        );
    }
}
