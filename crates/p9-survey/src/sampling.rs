//! Monte Carlo draws of a catalogued orbit solution.
//!
//! One draw perturbs each element by its 1σ spread (clamped to the physical
//! range the sky maps use) and places the planet at a uniformly random mean
//! anomaly — equal probability per unit time, so the draws form a
//! dwell-weighted "where is it now" prior.

use p9_core::constants::DEG2RAD;
use p9_core::types::OrbitalElements;
use rand::Rng;
use rand_distr::{Distribution, Normal};

use crate::schema::OrbitSolution;

/// Element clamps shared by every sky map in the workspace.
pub const A_RANGE_AU: (f64, f64) = (150.0, 1500.0);
pub const E_RANGE: (f64, f64) = (0.01, 0.9);
pub const I_RANGE_DEG: (f64, f64) = (0.0, 60.0);

/// Gaussian draw of the orientation angles (ω, Ω) of `sol`, in radians.
pub fn draw_orientation<R: Rng>(sol: &OrbitSolution, rng: &mut R) -> (f64, f64) {
    let omega = Normal::new(sol.omega_deg, sol.omega_sigma_deg)
        .expect("positive ω spread")
        .sample(rng);
    let omega_big = Normal::new(sol.omega_big_deg, sol.omega_big_sigma_deg)
        .expect("positive Ω spread")
        .sample(rng);
    (omega * DEG2RAD, omega_big * DEG2RAD)
}

/// One draw of the full element set of `sol` at a random orbital phase.
pub fn draw_elements<R: Rng>(sol: &OrbitSolution, rng: &mut R) -> OrbitalElements {
    let normal = |mean: f64, sigma: f64| Normal::new(mean, sigma).expect("positive spread");
    let a = normal(sol.a_au, sol.a_sigma_au)
        .sample(rng)
        .clamp(A_RANGE_AU.0, A_RANGE_AU.1);
    let e = normal(sol.e, sol.e_sigma)
        .sample(rng)
        .clamp(E_RANGE.0, E_RANGE.1);
    let i = normal(sol.i_deg, sol.i_sigma_deg)
        .sample(rng)
        .clamp(I_RANGE_DEG.0, I_RANGE_DEG.1);
    let (omega, omega_big) = draw_orientation(sol, rng);
    OrbitalElements {
        a,
        e,
        i: i * DEG2RAD,
        omega,
        omega_big,
        mean_anomaly: rng.gen_range(0.0..std::f64::consts::TAU),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::studies::catalog;
    use rand::SeedableRng;

    #[test]
    fn draws_respect_the_clamps_and_centre_on_the_solution() {
        let sol = catalog()
            .into_iter()
            .find(|s| s.name.contains("2021"))
            .unwrap();
        let mut rng = rand::rngs::StdRng::seed_from_u64(9);
        let n = 20_000;
        let mut a_sum = 0.0;
        for _ in 0..n {
            let el = draw_elements(&sol, &mut rng);
            assert!((A_RANGE_AU.0..=A_RANGE_AU.1).contains(&el.a));
            assert!((E_RANGE.0..=E_RANGE.1).contains(&el.e));
            assert!(el.i >= 0.0 && el.i <= I_RANGE_DEG.1 * DEG2RAD);
            a_sum += el.a;
        }
        assert!((a_sum / n as f64 - sol.a_au).abs() < 5.0);
    }
}
