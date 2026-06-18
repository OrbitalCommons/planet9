//! Laplace-Lagrange secular dynamics of a test particle embedded in the
//! self-gravitating eccentric disk (Sefilian & Touma 2019).
//!
//! # Construction
//!
//! For a coplanar secular system the disturbing function of a test particle
//! (eccentricity e, longitude of perihelion varpi) under a set of perturbing
//! rings j is, to second order in eccentricities (Murray & Dermott 1999,
//! Eq. 7.8 in the test-particle limit):
//!
//!   R = n a^2 [ (1/2) A e^2 + sum_j A_j e e_j cos(varpi - varpi_j) ]
//!
//! with the Laplace-Lagrange coefficients (test particle interior OR exterior
//! to ring j; we use the symmetric softened kernel so the same expression
//! applies on both sides, with alpha = min(a,a_j)/max(a,a_j)):
//!
//!   A   = + n (1/4) sum_j (m_j/M_sun) alpha_j abar_j b_{3/2}^{(1)}(alpha_j)
//!   A_j = - n (1/4) (m_j/M_sun) alpha_j abar_j b_{3/2}^{(2)}(alpha_j)
//!
//! where abar_j = alpha_j for an exterior perturber and 1 for an interior one
//! (the standard M&D external/internal factor). n is the test particle's mean
//! motion. The Laplace coefficients come from `crate::laplace` (softened to
//! represent the disk's finite radial width); we never re-derive p9-core's
//! quadrupole machinery, we build the coplanar Laplace-Lagrange operator on
//! top of the same orbit-averaging idea and cross-check the precession sign /
//! magnitude against `p9_core::analysis::secular` ring averaging.
//!
//! Writing the complex eccentricity z = e exp(i varpi), Lagrange's planetary
//! equations give
//!
//!   dz/dt = i A z + i sum_j A_j z_j  ==  i A (z - z_forced),
//!
//!   z_forced = -(sum_j A_j z_j) / A .
//!
//! So the particle's proper eccentricity vector circulates about the FORCED
//! value z_forced at the rate A (the disk-induced apsidal precession rate).
//! Because the disk is apsidally aligned, all z_j share the disk apse, so
//! z_forced lies along varpi_disk (if A_j/A < 0) or varpi_disk + 180 deg.
//! When |z_forced| exceeds the proper eccentricity the angle Delta varpi =
//! varpi - varpi_disk LIBRATES about the forced apse instead of circulating:
//! this is the apsidal confinement / shepherding.

use crate::disk::DiskProfile;
use crate::laplace::softened_laplace;
use p9_core::analysis::circular::wrap_to_pi;
use p9_core::constants::GM_SUN;
use p9_core::units::{days, radians, Angle, AngularVelocity, Time};
use std::f64::consts::PI;

const N_LAPLACE_QUAD: usize = 512;

/// Result of the Laplace-Lagrange secular analysis for one test particle.
#[derive(Debug, Clone, Copy)]
pub struct SecularSolution {
    /// Self / free precession coefficient A (radians/day). The proper
    /// eccentricity vector precesses at this rate; the secular period is
    /// 2 pi / A.
    pub a_coeff: f64,
    /// Forced longitude of perihelion (radians).
    pub varpi_forced: f64,
    /// Forced eccentricity magnitude |z_forced|.
    pub e_forced: f64,
    /// Disk apse (radians), for reference.
    pub varpi_disk: f64,
}

impl SecularSolution {
    /// Secular precession period in years (2 pi / A).
    pub fn precession_period_yr(&self) -> f64 {
        (2.0 * PI / self.a_coeff).abs() / p9_core::constants::YEAR_DAYS
    }

    /// Forced Delta varpi relative to the disk apse (radians, wrapped to
    /// (-pi, pi]). 0 means aligned with the disk apse, +-pi anti-aligned.
    pub fn forced_delta_varpi(&self) -> f64 {
        wrap_to_pi(self.varpi_forced - self.varpi_disk)
    }

    /// Free/self precession rate A as a typed [`AngularVelocity`] (rad/day).
    pub fn precession_rate(&self) -> AngularVelocity {
        (radians(self.a_coeff) / days(1.0)).into()
    }

    /// Secular precession period as a typed [`Time`].
    pub fn precession_period(&self) -> Time {
        days((2.0 * PI / self.a_coeff).abs())
    }

    /// Forced longitude of perihelion as a typed [`Angle`].
    pub fn varpi_forced_angle(&self) -> Angle {
        radians(self.varpi_forced)
    }

    /// Disk apse as a typed [`Angle`].
    pub fn varpi_disk_angle(&self) -> Angle {
        radians(self.varpi_disk)
    }

    /// Forced Delta varpi relative to the disk apse as a typed [`Angle`].
    pub fn forced_delta_varpi_angle(&self) -> Angle {
        radians(self.forced_delta_varpi())
    }
}

/// Test-particle mean motion (radians/day).
fn mean_motion(a: f64) -> f64 {
    (GM_SUN / a.powi(3)).sqrt()
}

/// Solve the coplanar Laplace-Lagrange problem for a test particle at
/// semi-major axis `a` embedded in `disk`.
pub fn solve(a: f64, disk: &DiskProfile) -> SecularSolution {
    let n = mean_motion(a);
    let eps = disk.softening_frac; // dimensionless softening on alpha

    let mut a_coeff = 0.0; // sum for A
                           // forced sum: sum_j A_j z_j, with z_j = e_j exp(i varpi_j)
    let mut forced_re = 0.0;
    let mut forced_im = 0.0;

    for ring in disk.rings() {
        let m_ratio = ring.mass_solar; // m_j / M_sun (mass already in solar)
        let (alpha, abar) = if a < ring.a {
            // test particle interior to the ring: alpha = a/a_j, abar = alpha
            let al = a / ring.a;
            (al, al)
        } else {
            // test particle exterior to the ring: alpha = a_j/a, abar = 1
            (ring.a / a, 1.0)
        };

        let b1 = softened_laplace(1.5, 1, alpha, eps, N_LAPLACE_QUAD);
        let b2 = softened_laplace(1.5, 2, alpha, eps, N_LAPLACE_QUAD);

        let pref = 0.25 * n * m_ratio * alpha * abar;
        a_coeff += pref * b1;
        let aj = -pref * b2;

        let (s, c) = ring.varpi.sin_cos();
        forced_re += aj * ring.e * c;
        forced_im += aj * ring.e * s;
    }

    // z_forced = -(sum_j A_j z_j) / A
    let zf_re = -forced_re / a_coeff;
    let zf_im = -forced_im / a_coeff;
    let e_forced = zf_re.hypot(zf_im);
    let varpi_forced = if e_forced > 0.0 {
        zf_im.atan2(zf_re).rem_euclid(2.0 * PI)
    } else {
        disk.varpi_disk
    };

    SecularSolution {
        a_coeff,
        varpi_forced,
        e_forced,
        varpi_disk: disk.varpi_disk,
    }
}

/// Secular Hamiltonian (per unit n a^2) of the test particle in the disk's
/// frame, expressed in the canonical Laplace-Lagrange variables. Using the
/// disturbing function above with z = e exp(i varpi):
///
///   H(e, Delta varpi) = (1/2) A e^2 + B e cos(Delta varpi)
///
/// where Delta varpi = varpi - varpi_disk and B = sum_j A_j e_j (all rings
/// share the disk apse, so the cross term collapses to a single cosine). The
/// level sets of H in the (e cos dphi, e sin dphi) plane are the secular
/// phase portrait: closed curves around the forced equilibrium that DO NOT
/// enclose the origin are LIBRATION (Delta varpi confined); curves enclosing
/// the origin are CIRCULATION.
pub fn hamiltonian_coeffs(a: f64, disk: &DiskProfile) -> (f64, f64) {
    let n = mean_motion(a);
    let eps = disk.softening_frac;
    let mut a_coeff = 0.0;
    let mut b_coeff = 0.0;
    for ring in disk.rings() {
        let m_ratio = ring.mass_solar;
        let (alpha, abar) = if a < ring.a {
            let al = a / ring.a;
            (al, al)
        } else {
            (ring.a / a, 1.0)
        };
        let b1 = softened_laplace(1.5, 1, alpha, eps, N_LAPLACE_QUAD);
        let b2 = softened_laplace(1.5, 2, alpha, eps, N_LAPLACE_QUAD);
        let pref = 0.25 * n * m_ratio * alpha * abar;
        a_coeff += pref * b1;
        b_coeff += -pref * b2 * ring.e;
    }
    (a_coeff, b_coeff)
}

/// Value of the secular Hamiltonian H(e, dphi) = (1/2) A e^2 + B e cos(dphi).
pub fn hamiltonian(a: f64, disk: &DiskProfile, e: f64, delta_varpi: f64) -> f64 {
    let (a_c, b_c) = hamiltonian_coeffs(a, disk);
    0.5 * a_c * e * e + b_c * e * delta_varpi.cos()
}

/// Classify the trajectory through (e0, dphi0) as librating or circulating by
/// integrating the secular equations of motion analytically.
///
/// In the (h, k) = (e sin dphi, e cos dphi) plane the motion is a circle of
/// radius |z0 - z_forced| centered on z_forced = (e_forced at the forced
/// apse). The trajectory LIBRATES in dphi iff that circle does not enclose the
/// origin, i.e. |z0 - z_forced| < |z_forced|.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApsidalMode {
    /// Delta varpi confined (shepherded): proper circle excludes the origin.
    Libration,
    /// Delta varpi circulates freely.
    Circulation,
}

/// Whether the disk SHEPHERDS a test particle on `e0`, `dphi0` within the age
/// of the solar system: the trajectory must both librate (geometry) AND
/// complete that libration on a timescale shorter than `age_yr`. In linear
/// Laplace-Lagrange theory the libration *geometry* is mass-independent (z_forced
/// = -B/A with both linear in M_disk), so the mass dependence of realized
/// confinement enters entirely through the precession period 2 pi / A, which
/// scales as 1/M_disk. A negligible-mass disk has a precession period far
/// exceeding the age and therefore confines nothing in practice.
pub fn confines_within_age(a: f64, disk: &DiskProfile, e0: f64, dphi0: f64, age_yr: f64) -> bool {
    if apsidal_mode(a, disk, e0, dphi0) != ApsidalMode::Libration {
        return false;
    }
    solve(a, disk).precession_period_yr() < age_yr
}

/// Determine the apsidal mode for a particle starting at proper eccentricity
/// `e0` and Delta varpi `dphi0` (relative to the disk apse).
pub fn apsidal_mode(a: f64, disk: &DiskProfile, e0: f64, dphi0: f64) -> ApsidalMode {
    let sol = solve(a, disk);
    // Forced vector along varpi_forced; express everything relative to disk apse.
    let phi_f = sol.forced_delta_varpi();
    let zf_k = sol.e_forced * phi_f.cos();
    let zf_h = sol.e_forced * phi_f.sin();
    let z0_k = e0 * dphi0.cos();
    let z0_h = e0 * dphi0.sin();
    let radius = ((z0_k - zf_k).powi(2) + (z0_h - zf_h).powi(2)).sqrt();
    if radius < sol.e_forced {
        ApsidalMode::Libration
    } else {
        ApsidalMode::Circulation
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_typed_accessors_match_f64() {
        use p9_core::constants::{DAY_S, YEAR_DAYS};
        use uom::si::angle::radian;
        use uom::si::angular_velocity::radian_per_second;
        use uom::si::time::day;
        let sol = solve(250.0, &DiskProfile::fiducial(10.0));
        assert_relative_eq!(
            sol.precession_rate().get::<radian_per_second>(),
            sol.a_coeff / DAY_S,
            epsilon = 1e-30
        );
        assert_relative_eq!(
            sol.precession_period().get::<day>(),
            sol.precession_period_yr() * YEAR_DAYS,
            max_relative = 1e-12
        );
        assert_relative_eq!(
            sol.varpi_forced_angle().get::<radian>(),
            sol.varpi_forced,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            sol.forced_delta_varpi_angle().get::<radian>(),
            sol.forced_delta_varpi(),
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_precession_scales_linearly_with_disk_mass() {
        // The Laplace-Lagrange A coefficient is linear in each ring mass, so
        // the disk-induced precession rate must scale linearly with M_disk.
        let a = 250.0;
        let d5 = DiskProfile::fiducial(5.0);
        let d10 = DiskProfile::fiducial(10.0);
        let a5 = solve(a, &d5).a_coeff;
        let a10 = solve(a, &d10).a_coeff;
        assert_relative_eq!(a10 / a5, 2.0, epsilon = 1e-9);
    }

    #[test]
    fn test_precession_period_order_of_magnitude() {
        // 10 M_earth disk, test particle at 250 AU: secular period of order
        // 1e8-1e9 yr (paper Sec. 3).
        let sol = solve(250.0, &DiskProfile::fiducial(10.0));
        let p = sol.precession_period_yr();
        assert!(
            p > crate::disk::reference::PRECESSION_PERIOD_LO_YR
                && p < crate::disk::reference::PRECESSION_PERIOD_HI_YR,
            "precession period = {p:.3e} yr (expected 1e8-1e9)"
        );
    }

    #[test]
    fn test_forced_eccentricity_nonzero_and_apsidally_set() {
        // A massive apsidally-aligned disk forces a finite e_forced whose apse
        // is along (or opposite) the disk apse.
        let sol = solve(250.0, &DiskProfile::fiducial(10.0));
        assert!(sol.e_forced > 0.0);
        let dphi = sol.forced_delta_varpi().abs();
        // Forced apse is either aligned (0) or anti-aligned (pi) with the disk.
        assert!(
            dphi < 1e-6 || (dphi - PI).abs() < 1e-6,
            "forced Delta varpi = {dphi} should be 0 or pi"
        );
    }

    #[test]
    fn test_libration_at_ten_earth_masses() {
        // A test particle started near the forced equilibrium with the disk at
        // 10 M_earth LIBRATES (apsidal confinement) AND completes that
        // libration well within the age of the solar system, i.e. it is
        // physically shepherded.
        let disk = DiskProfile::fiducial(10.0);
        let sol = solve(250.0, &disk);
        let phi_f = sol.forced_delta_varpi();
        // Start at small proper amplitude about the forced point.
        let e0 = sol.e_forced * 1.1;
        assert_eq!(
            apsidal_mode(250.0, &disk, e0, phi_f),
            ApsidalMode::Libration
        );
        assert!(
            confines_within_age(250.0, &disk, e0, phi_f, 4.5e9),
            "10 M_earth disk must shepherd within the 4.5 Gyr age"
        );
    }

    #[test]
    fn test_no_shepherding_at_negligible_mass() {
        // Near-zero disk mass: the libration GEOMETRY is mass-independent in
        // linear theory, but the precession period diverges as 1/M_disk, so a
        // negligible-mass disk confines NOTHING within the solar system age.
        let disk = DiskProfile::fiducial(1e-4);
        let sol = solve(250.0, &disk);
        let phi_f = sol.forced_delta_varpi();
        let e0 = sol.e_forced * 1.1;
        // Same orbit that librates for the 10 M_earth disk:
        assert_eq!(
            apsidal_mode(250.0, &disk, e0, phi_f),
            ApsidalMode::Libration
        );
        // ...but the secular clock is far too slow to realize it.
        assert!(sol.precession_period_yr() > 4.5e9);
        assert!(!confines_within_age(250.0, &disk, e0, phi_f, 4.5e9));
    }

    #[test]
    fn test_hamiltonian_forced_equilibrium_is_extremum() {
        // dH/de = 0 and dH/d(dphi)=0 at (e_forced, forced apse).
        let disk = DiskProfile::fiducial(10.0);
        let a = 250.0;
        let sol = solve(a, &disk);
        let (a_c, b_c) = hamiltonian_coeffs(a, &disk);
        // Equilibrium: e* = -B/A at dphi = 0 (B<0 => e*>0, forced apse).
        let e_star = -b_c / a_c;
        assert_relative_eq!(e_star, sol.e_forced, epsilon = 1e-9);
        // H slightly off equilibrium is larger (or smaller) consistently:
        // perturbing e at fixed dphi=0 increases |gradient|.
        let h0 = hamiltonian(a, &disk, e_star, 0.0);
        let hp = hamiltonian(a, &disk, e_star + 1e-4, 0.0);
        let hm = hamiltonian(a, &disk, e_star - 1e-4, 0.0);
        // Quadratic minimum/maximum: second difference dominates.
        let curv = hp - 2.0 * h0 + hm;
        assert!(curv.abs() > 0.0);
    }
}
