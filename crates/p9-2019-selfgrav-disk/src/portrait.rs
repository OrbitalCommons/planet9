//! Secular phase portrait in the (Delta varpi, e) plane, demonstrating the
//! libration island that confines test ETNOs to the disk apse.
//!
//! The level sets of H(e, Delta varpi) = (1/2) A e^2 + B e cos(Delta varpi)
//! (see `crate::secular`) separate librating from circulating trajectories.
//! The separatrix passes through the origin (e = 0). Any contour with
//! H < H_separatrix = 0 (for B<0, A>0) encircles only the forced equilibrium
//! and is a libration island.

use crate::disk::DiskProfile;
use crate::secular::{hamiltonian, hamiltonian_coeffs};
use std::f64::consts::PI;

/// Grid of H over (Delta varpi, e). Returns (dphi_vals, e_vals, grid) with
/// grid[i][j] = H(e_i, dphi_j).
pub fn phase_portrait(
    a: f64,
    disk: &DiskProfile,
    n_e: usize,
    n_phi: usize,
    e_max: f64,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let e_vals: Vec<f64> = (0..n_e)
        .map(|i| (i as f64 + 0.5) / n_e as f64 * e_max)
        .collect();
    let phi_vals: Vec<f64> = (0..n_phi)
        .map(|j| j as f64 / n_phi as f64 * 2.0 * PI - PI)
        .collect();
    let mut grid = vec![vec![0.0; n_phi]; n_e];
    for (i, &e) in e_vals.iter().enumerate() {
        for (j, &phi) in phi_vals.iter().enumerate() {
            grid[i][j] = hamiltonian(a, disk, e, phi);
        }
    }
    (phi_vals, e_vals, grid)
}

/// Fraction of a uniform grid of initial conditions (in the e <= e_max,
/// Delta varpi in (-pi,pi] box) whose apsidal mode is LIBRATION (pure secular
/// geometry). In linear Laplace-Lagrange theory this fraction depends only on
/// the forced-eccentricity geometry z_forced = -B/A, which is INDEPENDENT of
/// the disk mass (B and A both scale linearly with M_disk). It measures the
/// size of the libration island, not whether it is dynamically realized.
pub fn libration_fraction(a: f64, disk: &DiskProfile, n_e: usize, n_phi: usize, e_max: f64) -> f64 {
    use crate::secular::{apsidal_mode, ApsidalMode};
    let mut lib = 0;
    let mut tot = 0;
    for i in 0..n_e {
        let e = (i as f64 + 0.5) / n_e as f64 * e_max;
        for j in 0..n_phi {
            let phi = j as f64 / n_phi as f64 * 2.0 * PI - PI;
            if apsidal_mode(a, disk, e, phi) == ApsidalMode::Libration {
                lib += 1;
            }
            tot += 1;
        }
    }
    lib as f64 / tot as f64
}

/// Fraction of the same grid that is actually SHEPHERDED within `age_yr`:
/// librating AND precessing fast enough (period < age). This is the
/// physically realized confinement and GROWS with disk mass (the period
/// scales as 1/M_disk), vanishing as M_disk -> 0.
pub fn shepherded_fraction(
    a: f64,
    disk: &DiskProfile,
    n_e: usize,
    n_phi: usize,
    e_max: f64,
    age_yr: f64,
) -> f64 {
    use crate::secular::confines_within_age;
    let mut shep = 0;
    let mut tot = 0;
    for i in 0..n_e {
        let e = (i as f64 + 0.5) / n_e as f64 * e_max;
        for j in 0..n_phi {
            let phi = j as f64 / n_phi as f64 * 2.0 * PI - PI;
            if confines_within_age(a, disk, e, phi, age_yr) {
                shep += 1;
            }
            tot += 1;
        }
    }
    shep as f64 / tot as f64
}

/// Location of the separatrix energy: H at the origin, H(0, *) = 0. Provided
/// for documentation of the libration/circulation split.
pub fn separatrix_energy(a: f64, disk: &DiskProfile) -> f64 {
    let _ = hamiltonian_coeffs(a, disk);
    0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_portrait_shape_and_finite() {
        let disk = DiskProfile::fiducial(10.0);
        let (phi, e, grid) = phase_portrait(250.0, &disk, 10, 16, 0.6);
        assert_eq!(phi.len(), 16);
        assert_eq!(e.len(), 10);
        assert_eq!(grid.len(), 10);
        for row in &grid {
            for &h in row {
                assert!(h.is_finite());
            }
        }
    }

    #[test]
    fn test_libration_island_exists_and_is_mass_independent() {
        // The libration island (pure secular geometry) is present for a finite
        // eccentric disk and its SIZE does not depend on disk mass in linear
        // Laplace-Lagrange theory (z_forced = -B/A, both linear in M_disk).
        let a = 250.0;
        let f_heavy = libration_fraction(a, &DiskProfile::fiducial(10.0), 12, 24, 0.6);
        let f_light = libration_fraction(a, &DiskProfile::fiducial(0.1), 12, 24, 0.6);
        assert!(
            f_heavy > 0.0,
            "10 M_earth disk must produce a libration island"
        );
        assert!(
            (f_heavy - f_light).abs() < 1e-9,
            "island size is mass-independent"
        );
    }

    #[test]
    fn test_shepherded_fraction_grows_with_mass() {
        // The physically REALIZED confinement (libration completed within the
        // age) grows with disk mass, because the precession period ~ 1/M_disk.
        let a = 250.0;
        let f_heavy = shepherded_fraction(a, &DiskProfile::fiducial(10.0), 12, 24, 0.6, 4.5e9);
        let f_light = shepherded_fraction(a, &DiskProfile::fiducial(0.1), 12, 24, 0.6, 4.5e9);
        assert!(
            f_heavy > f_light,
            "shepherded fraction should grow with disk mass: heavy {f_heavy}, light {f_light}"
        );
        assert!(
            f_heavy > 0.0,
            "10 M_earth disk must shepherd a finite fraction"
        );
    }

    #[test]
    fn test_nothing_shepherded_at_negligible_mass() {
        // Vanishing disk: the precession clock is far slower than the age, so
        // nothing is shepherded even though the island geometry persists.
        let f = shepherded_fraction(250.0, &DiskProfile::fiducial(1e-6), 12, 24, 0.6, 4.5e9);
        assert_eq!(f, 0.0);
    }
}
