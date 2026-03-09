//! Phase-space portrait generation — analytical and N-body.
//!
//! Reproduces Figures 3 and 4 of Batygin & Brown (2016):
//! - Figure 3: analytical secular phase portraits at 6 semi-major axes
//! - Figure 4: N-body simulated phase portraits showing anti-aligned orbits

use std::f64::consts::PI;

use p9_core::constants::*;
use p9_core::integrator::whm::WhmIntegrator;
use p9_core::types::*;

/// Result of tracing a single trajectory in (e, Δϖ) space.
#[derive(Debug, Clone)]
pub struct TrajectoryPoint {
    pub t: f64,
    pub e: f64,
    pub delta_varpi: f64,
    pub a: f64,
}

/// Run an N-body phase-space portrait for a single semi-major axis.
///
/// Creates `n_trajectories` test particles at the given semi-major axis
/// with eccentricities linearly spaced from 0.025 to 0.95, half starting
/// at Δϖ=0 and half at Δϖ=π. Integrates for `t_total` days.
///
/// Returns the trajectory points for all particles.
pub fn nbody_phase_portrait(
    a_test: f64,
    p9: &P9Params,
    bodies: &mut Vec<MassiveBody>,
    n_e_steps: usize,
    t_total: f64,
    dt: f64,
    snapshot_interval: f64,
) -> Vec<TrajectoryPoint> {
    let p9_body = p9.to_body();
    let varpi_p9 = p9.omega + p9.omega_big;

    // Create test particles at two starting Δϖ values: 0 and π
    let mut particles = Vec::new();
    let mut active = Vec::new();

    for i in 0..n_e_steps {
        let e = (i as f64 + 0.5) / n_e_steps as f64 * 0.95;

        for &dv_start in &[0.0, PI] {
            let varpi = varpi_p9 + dv_start;
            let omega = varpi; // For coplanar, Ω=0, so ω=ϖ
            let elem = OrbitalElements {
                a: a_test,
                e,
                i: 0.0,
                omega_big: 0.0,
                omega,
                mean_anomaly: 0.0,
            };
            particles.push(elements_to_cartesian(&elem, GM_SUN));
            active.push(true);
        }
    }

    // Add P9 to bodies
    bodies.push(p9_body);

    let config = SimConfig {
        dt,
        t_start: 0.0,
        t_end: t_total,
        removal_inner_au: 5.0,
        removal_outer_au: 10_000.0,
        snapshot_interval_days: snapshot_interval,
        hybrid_changeover_hill: 3.0,
        bs_epsilon: 1e-11,
    };

    let integrator = WhmIntegrator::new();
    let mut results = Vec::new();
    let mut t = 0.0;
    let mut next_snapshot = 0.0;

    let n_steps = (t_total / dt).ceil() as usize;
    for _ in 0..n_steps {
        integrator.step(bodies, &mut particles, &mut active, dt, &config);
        t += dt;

        if t >= next_snapshot {
            // Record orbital elements for all active particles
            for (i, particle) in particles.iter().enumerate() {
                if !active[i] {
                    continue;
                }
                let elem = cartesian_to_elements(particle, GM_SUN);

                // Compute Δϖ relative to Planet Nine's current position
                let p9_elem = cartesian_to_elements(&bodies.last().unwrap().state, GM_SUN);
                let varpi_p9_current = p9_elem.omega + p9_elem.omega_big;
                let varpi_test = elem.omega + elem.omega_big;
                let mut dv = varpi_test - varpi_p9_current;
                while dv > PI {
                    dv -= 2.0 * PI;
                }
                while dv < -PI {
                    dv += 2.0 * PI;
                }

                results.push(TrajectoryPoint {
                    t,
                    e: elem.e,
                    delta_varpi: dv,
                    a: elem.a,
                });
            }
            next_snapshot += snapshot_interval;
        }
    }

    // Remove P9 from bodies (restore original)
    bodies.pop();

    results
}

/// Classify trajectory points into aligned vs anti-aligned populations.
///
/// Anti-aligned: |Δϖ| > π/2 (perihelion on opposite side from P9)
/// Aligned: |Δϖ| < π/2 (perihelion on same side as P9)
pub fn classify_alignment(
    points: &[TrajectoryPoint],
) -> (Vec<&TrajectoryPoint>, Vec<&TrajectoryPoint>) {
    let mut aligned = Vec::new();
    let mut anti_aligned = Vec::new();

    for p in points {
        if p.delta_varpi.abs() > PI / 2.0 {
            anti_aligned.push(p);
        } else {
            aligned.push(p);
        }
    }

    (aligned, anti_aligned)
}

/// Run phase portraits at the 6 semi-major axes from the paper.
/// Returns a Vec of (a, trajectory_points) for each.
pub fn paper_phase_portraits(
    p9: &P9Params,
    t_total: f64,
    dt: f64,
) -> Vec<(f64, Vec<TrajectoryPoint>)> {
    let test_axes = [50.0, 150.0, 250.0, 350.0, 450.0, 550.0];

    let mut results = Vec::new();

    for &a in &test_axes {
        // Start with no other bodies (pure P9 + Sun for the phase portrait)
        let mut bodies = Vec::new();

        let points = nbody_phase_portrait(a, p9, &mut bodies, 20, t_total, dt, t_total / 100.0);

        results.push((a, points));
    }

    results
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nbody_portrait_runs() {
        let p9 = P9Params::nominal_2016();
        let mut bodies = Vec::new();

        // Short integration for testing
        let points = nbody_phase_portrait(
            250.0,
            &p9,
            &mut bodies,
            5,     // 5 eccentricity steps (10 particles total)
            1e6,   // 1 Myr
            300.0, // 300 day timestep
            1e5,   // snapshot every 100 kyr
        );

        assert!(!points.is_empty(), "Should produce trajectory points");

        for p in &points {
            assert!(
                p.e >= 0.0 && p.e < 1.0,
                "Eccentricity out of range: {}",
                p.e
            );
            assert!(
                p.delta_varpi >= -PI && p.delta_varpi <= PI,
                "Δϖ out of range: {}",
                p.delta_varpi
            );
        }
    }

    #[test]
    fn test_alignment_classification() {
        let points = vec![
            TrajectoryPoint {
                t: 0.0,
                e: 0.5,
                delta_varpi: 0.1,
                a: 300.0,
            },
            TrajectoryPoint {
                t: 0.0,
                e: 0.5,
                delta_varpi: 2.5,
                a: 300.0,
            },
            TrajectoryPoint {
                t: 0.0,
                e: 0.5,
                delta_varpi: -2.5,
                a: 300.0,
            },
        ];

        let (aligned, anti) = classify_alignment(&points);
        assert_eq!(aligned.len(), 1);
        assert_eq!(anti.len(), 2);
    }
}
