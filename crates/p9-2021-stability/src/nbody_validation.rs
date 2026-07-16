//! N-body validation of the analytic diffusion model (M4).
//!
//! Integrates test particles with Neptune on a circular orbit using the
//! p9-core Wisdom-Holman integrator and feeds the resulting a(t) series into
//! the MSD diffusion estimator, so the analytic
//! `p9_core::analysis::resonance::neptune_diffusion_coefficient` is checked against real dynamics
//! instead of against itself (the previous "validation" set
//! `d_measured = d_analytical` and asserted ratio ≈ 1).
//!
//! Reduced scale by default (few particles, ~1 Myr); the full-scale
//! comparison lives behind `#[ignore]`.

use nalgebra::Vector3;
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use p9_core::constants::{A_NEPTUNE_AU, DEG2RAD, GM_NEPTUNE, GM_SUN, MASS_NEPTUNE_SOLAR, TWO_PI};
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::{cartesian_to_elements, MassiveBody, OrbitalElements, SimConfig, StateVector};

/// Neptune as a massive body on a circular, planar orbit at a_N.
pub fn neptune_circular() -> MassiveBody {
    let v = ((GM_SUN + GM_NEPTUNE) / A_NEPTUNE_AU).sqrt();
    MassiveBody {
        name: "Neptune".to_string(),
        gm: GM_NEPTUNE,
        mass: MASS_NEPTUNE_SOLAR,
        state: StateVector::new(
            Vector3::new(A_NEPTUNE_AU, 0.0, 0.0),
            Vector3::new(0.0, v, 0.0),
        ),
        radius_au: 1.655e-4,
        j2: None,
        j4: None,
    }
}

/// Integrate `n_particles` test particles started at (a0, q0) with random
/// phases under Sun + Neptune, recording a(t) every `snapshot_every` steps.
///
/// Returns one a(t) series per particle (truncated if the particle is
/// removed). Seeded and deterministic.
pub fn neptune_scattering_series(
    a0: f64,
    q0: f64,
    n_particles: usize,
    t_total_days: f64,
    dt_days: f64,
    snapshot_every: usize,
    seed: u64,
) -> Vec<Vec<f64>> {
    let mut rng = StdRng::seed_from_u64(seed);
    let e0 = 1.0 - q0 / a0;

    let mut particles: Vec<StateVector> = (0..n_particles)
        .map(|_| {
            let elem = OrbitalElements {
                a: a0,
                e: e0,
                i: rng.gen_range(0.0..5.0) * DEG2RAD,
                omega: rng.gen_range(0.0..TWO_PI),
                omega_big: rng.gen_range(0.0..TWO_PI),
                mean_anomaly: rng.gen_range(0.0..TWO_PI),
            };
            elem.to_state_vector(GM_SUN)
        })
        .collect();
    let mut active = vec![true; n_particles];
    let mut bodies = vec![neptune_circular()];

    let config = SimConfig {
        dt: dt_days,
        t_start: 0.0,
        t_end: t_total_days,
        removal_inner_au: 1.0,
        removal_outer_au: 1.0e5,
        snapshot_interval_days: f64::INFINITY,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };
    let integrator = WhmIntegrator::new();

    let n_steps = (t_total_days / dt_days).round() as usize;
    let mut series: Vec<Vec<f64>> = vec![vec![a0]; n_particles];

    for step in 0..n_steps {
        integrator.step(&mut bodies, &mut particles, &mut active, dt_days, &config);
        if (step + 1) % snapshot_every == 0 {
            for (k, p) in particles.iter().enumerate() {
                if active[k] {
                    let elem = cartesian_to_elements(p, GM_SUN);
                    if elem.e < 1.0 && elem.a > 0.0 {
                        series[k].push(elem.a);
                    } else {
                        active[k] = false;
                    }
                }
            }
        }
    }

    series
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diffusion::{measure_diffusion, measurement_from_ensemble};
    use p9_core::constants::YEAR_DAYS;

    /// Reduced-scale check: real Neptune scattering at (a, q) = (250, 33) AU
    /// produces measurable semi-major-axis spreading, and the MSD estimator
    /// agrees with the analytic D_a(q) at the order-of-magnitude level the
    /// paper claims for its diffusion model.
    #[test]
    fn test_nbody_diffusion_vs_analytic_reduced() {
        let dt = 2500.0; // days; Neptune period / dt ≈ 24 steps/orbit
        let t_total = 1.0e6 * YEAR_DAYS; // 1 Myr
        let snapshot_every = 40; // every ~274 yr
        let ensemble = neptune_scattering_series(250.0, 33.0, 8, t_total, dt, snapshot_every, 11);

        let dt_yr = dt * snapshot_every as f64 / YEAR_DAYS;
        let m = measurement_from_ensemble(250.0, 33.0, &ensemble, dt_yr)
            .expect("ensemble long enough to fit");
        assert!(m.d_measured > 0.0, "no diffusion measured");
        let ratio = m.d_measured / m.d_analytical;
        assert!(
            (0.02..50.0).contains(&ratio),
            "D_measured/D_analytic = {ratio:.3} (measured {:.3e}, analytic {:.3e})",
            m.d_measured,
            m.d_analytical
        );
    }

    /// High-q control: at q = 42 AU (above the chaotic boundary at this a)
    /// the measured spreading must be far smaller than at q = 33 AU.
    #[test]
    fn test_nbody_diffusion_suppressed_at_high_q() {
        let dt = 2500.0;
        let t_total = 5.0e5 * YEAR_DAYS;
        let low_q = neptune_scattering_series(250.0, 33.0, 6, t_total, dt, 40, 13);
        let high_q = neptune_scattering_series(250.0, 42.0, 6, t_total, dt, 40, 13);

        let dt_yr = dt * 40.0 / YEAR_DAYS;
        let d_low = measure_diffusion(&low_q, dt_yr).unwrap().d;
        let d_high = measure_diffusion(&high_q, dt_yr).unwrap().d;
        assert!(
            d_low > 10.0 * d_high,
            "expected strong q suppression: D(33) = {d_low:.3e}, D(42) = {d_high:.3e}"
        );
    }

    /// Full-scale comparison (16 particles, 4 Myr): D within one dex of the
    /// analytic coefficient. Slow — run with `cargo test -- --ignored`.
    #[test]
    #[ignore]
    fn test_nbody_diffusion_vs_analytic_full() {
        let dt = 2500.0;
        let t_total = 4.0e6 * YEAR_DAYS;
        let ensemble = neptune_scattering_series(250.0, 33.0, 16, t_total, dt, 40, 17);
        let dt_yr = dt * 40.0 / YEAR_DAYS;
        let m = measurement_from_ensemble(250.0, 33.0, &ensemble, dt_yr).unwrap();
        let ratio = m.d_measured / m.d_analytical;
        assert!(
            (0.1..10.0).contains(&ratio),
            "D_measured/D_analytic = {ratio:.3}"
        );
    }
}
