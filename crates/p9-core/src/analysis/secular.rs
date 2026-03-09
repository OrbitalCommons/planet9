//! Secular perturbation theory.
//!
//! Implements the orbit-averaged Hamiltonian to quadrupole and octupole order
//! for generating phase-space portraits of the (e, Δϖ) dynamics.
//!
//! Reference: Batygin & Brown (2016), Equations 1-5

use std::f64::consts::PI;

/// Quadrupole-order secular Hamiltonian for a test particle perturbed by Planet Nine.
///
/// H_quad = -C * [2 + 3*e^2 - 3*(1 + 4/3*e^2)*sin^2(i)*sin^2(Δω)
///               - 5*e^2*cos^2(Δϖ)*(1 - sin^2(i)*sin^2(Δω))]
///
/// Simplified for the coplanar case (i = 0):
///   H_quad = -C * [2 + 3*e^2 - 5*e^2*cos^2(Δϖ)]
///
/// where C depends on the perturber mass, semi-major axis, eccentricity, and the
/// test particle's semi-major axis.
///
/// `a`: test particle semi-major axis (AU)
/// `e`: test particle eccentricity
/// `delta_varpi`: longitude of perihelion difference (radians)
/// `i`: test particle inclination (radians)
/// `delta_omega`: argument of perihelion difference (radians)
/// `a_p`: perturber semi-major axis (AU)
/// `e_p`: perturber eccentricity
/// `gm_p`: perturber GM (AU^3/day^2)
pub fn quadrupole_hamiltonian(
    a: f64,
    e: f64,
    delta_varpi: f64,
    i: f64,
    delta_omega: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
) -> f64 {
    let alpha = a / a_p;
    let e2 = e * e;
    let sin_i2 = i.sin().powi(2);
    let sin_dw2 = delta_omega.sin().powi(2);
    let cos_dv2 = delta_varpi.cos().powi(2);

    // Coupling constant
    let c = gm_p * alpha * alpha / (4.0 * a_p * (1.0 - e_p * e_p).powf(1.5));

    -c * (2.0 + 3.0 * e2
        - 3.0 * (1.0 + 4.0 / 3.0 * e2) * sin_i2 * sin_dw2
        - 5.0 * e2 * cos_dv2 * (1.0 - sin_i2 * sin_dw2))
}

/// Coplanar quadrupole Hamiltonian (i = 0).
/// H = -C * [2 + 3e^2 - 5e^2 cos^2(Δϖ)]
/// This is what generates the anti-aligned libration islands.
pub fn coplanar_quadrupole(a: f64, e: f64, delta_varpi: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    quadrupole_hamiltonian(a, e, delta_varpi, 0.0, 0.0, a_p, e_p, gm_p)
}

/// Generate a phase-space portrait grid: evaluate the Hamiltonian at a grid
/// of (e, Δϖ) values for a fixed semi-major axis.
///
/// Returns a 2D array where portrait[i][j] = H(e_i, dvarpi_j).
pub fn phase_portrait(
    a: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    n_e: usize,
    n_dvarpi: usize,
) -> (Vec<f64>, Vec<f64>, Vec<Vec<f64>>) {
    let e_vals: Vec<f64> = (0..n_e)
        .map(|i| (i as f64 + 0.5) / n_e as f64 * 0.95)
        .collect();

    let dvarpi_vals: Vec<f64> = (0..n_dvarpi)
        .map(|j| (j as f64) / n_dvarpi as f64 * 2.0 * PI - PI)
        .collect();

    let mut portrait = vec![vec![0.0; n_dvarpi]; n_e];

    for (i, &e) in e_vals.iter().enumerate() {
        for (j, &dv) in dvarpi_vals.iter().enumerate() {
            portrait[i][j] = coplanar_quadrupole(a, e, dv, a_p, e_p, gm_p);
        }
    }

    (e_vals, dvarpi_vals, portrait)
}
