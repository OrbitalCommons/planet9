//! Cross-validation of p9-core's Kepler propagation and element conversions
//! against starfield (a SPICE/Skyfield-parity port), used as a permanent
//! independent test oracle.
//!
//! - `kepler_drift` vs `starfield::keplerlib` (via `KeplerOrbit::at`, the
//!   public wrapper around the prop2b universal-variable propagator).
//! - `cartesian_to_elements` / `elements_to_cartesian` vs
//!   `starfield::elementslib::OsculatingElements`.

use nalgebra::Vector3;
use starfield::keplerlib::KeplerOrbit;
use starfield::time::Timescale;

use p9_core::constants::{AU_KM, DAY_S, GM_SUN, TWO_PI};
use p9_core::integrator::kepler_step::kepler_drift;
use p9_core::types::{cartesian_to_elements, elements_to_cartesian, OrbitalElements, StateVector};

/// p9-core's GM_SUN (AU^3/day^2) expressed in km^3/s^2 for starfield's
/// element APIs, so both sides describe the identical two-body problem.
fn gm_sun_km3_s2() -> f64 {
    GM_SUN * AU_KM * AU_KM * AU_KM / (DAY_S * DAY_S)
}

/// Propagate a state with starfield's prop2b port (universal variables +
/// Stumpff functions, SPICE parity) and return (pos AU, vel AU/day).
fn starfield_propagate(state: &StateVector, dt: f64) -> (Vector3<f64>, Vector3<f64>) {
    let ts = Timescale::default();
    let epoch = ts.tt_jd(2_451_545.0, None);
    let target = ts.tt_jd(2_451_545.0 + dt, None);
    let orbit = KeplerOrbit::new(state.pos, state.vel, &epoch, GM_SUN, None, None);
    let p = orbit.at(&target);
    (p.position, p.velocity)
}

// ===================================================================
// Work item 1: kepler_drift vs starfield::keplerlib::propagate
// ===================================================================

/// Compare kepler_drift to starfield over one (state, dt) case.
/// `tol` is a relative tolerance scaled by the orbit size (AU) and by the
/// number of orbits covered: the two independent universal-variable
/// formulations accumulate along-track phase differences of a few ulps per
/// revolution, so multi-orbit dt amplifies an irreducible rounding
/// disagreement linearly.
fn check_drift_against_starfield(elem: &OrbitalElements, dt: f64, tol: f64) {
    let state = elements_to_cartesian(elem, GM_SUN);
    let ours = kepler_drift(&state, dt, GM_SUN);
    let (sf_pos, sf_vel) = starfield_propagate(&state, dt);

    let period = TWO_PI * (elem.a.abs().powi(3) / GM_SUN).sqrt();
    let tol = tol * (1.0 + (dt / period).abs());

    let scale_r = state.pos.norm().max(sf_pos.norm());
    let scale_v = state.vel.norm().max(sf_vel.norm());

    let dpos = (ours.pos - sf_pos).norm();
    let dvel = (ours.vel - sf_vel).norm();
    assert!(
        dpos <= tol * scale_r,
        "position mismatch vs starfield: |dr| = {dpos:.3e} AU (scale {scale_r:.3e}), \
         a = {}, e = {}, M = {}, dt = {dt}",
        elem.a,
        elem.e,
        elem.mean_anomaly
    );
    assert!(
        dvel <= tol * scale_v,
        "velocity mismatch vs starfield: |dv| = {dvel:.3e} AU/d (scale {scale_v:.3e}), \
         a = {}, e = {}, M = {}, dt = {dt}",
        elem.a,
        elem.e,
        elem.mean_anomaly
    );
}

/// The documented failing corner of the old universal-variable iteration:
/// a = 800 AU, e = 0.95, M ≈ 5.97 (just before perihelion), dt = 4e5 days.
/// In debug builds the unfixed solver hit a convergence-plateau
/// debug_assert; the fixed solver must converge AND agree with starfield.
#[test]
fn test_kepler_drift_known_failing_corner() {
    let elem = OrbitalElements {
        a: 800.0,
        e: 0.95,
        i: 0.3,
        omega_big: 1.0,
        omega: 2.0,
        mean_anomaly: 5.97,
    };
    check_drift_against_starfield(&elem, 4.0e5, 1e-8);
}

/// Property grid: (a, e, M, dt) sweep spanning the integrator's operating
/// envelope (Kuiper belt through inner Oort cloud, fractions of an orbit
/// through many orbits, both time directions). Tolerance 1e-8 relative:
/// both solvers are independent universal-variable formulations, and
/// near-perihelion high-e states lose a few digits to cancellation in
/// either implementation.
#[test]
fn test_kepler_drift_grid_vs_starfield() {
    let a_grid = [5.0, 40.0, 250.0, 800.0, 3000.0];
    let e_grid = [0.0, 0.3, 0.7, 0.95, 0.99];
    let m_grid = [0.0, 1.0, 3.1, 5.97];
    // dt as a fraction of the orbital period, plus fixed offsets that put
    // long-period orbits in the "fraction of an orbit across perihelion"
    // regime where the plateau was observed.
    let dt_fracs = [0.05, 0.5, 1.0, 3.7, -0.5];

    for &a in &a_grid {
        let period = TWO_PI * (a * a * a / GM_SUN).sqrt();
        for &e in &e_grid {
            for &m in &m_grid {
                let elem = OrbitalElements {
                    a,
                    e,
                    i: 0.4,
                    omega_big: 0.7,
                    omega: 2.1,
                    mean_anomaly: m,
                };
                for &frac in &dt_fracs {
                    check_drift_against_starfield(&elem, frac * period, 1e-8);
                }
                // Fixed dt = 4e5 d (the oort-cloud integrator's step size).
                check_drift_against_starfield(&elem, 4.0e5, 1e-8);
            }
        }
    }
}

/// Hyperbolic states: kepler_drift vs starfield for e > 1.
#[test]
fn test_kepler_drift_hyperbolic_vs_starfield() {
    for &e in &[1.0 + 1e-6, 1.1, 1.5, 3.0] {
        for &q in &[5.0, 100.0, 800.0] {
            // Build a perihelion state directly: r = q x̂, v = v_peri ŷ.
            let v_peri = (GM_SUN * (1.0 + e) / q).sqrt();
            let state = StateVector::new(Vector3::new(q, 0.0, 0.0), Vector3::new(0.0, v_peri, 0.0));
            for &dt in &[-2.0e5, -1.0e3, 1.0e3, 2.0e5] {
                let ours = kepler_drift(&state, dt, GM_SUN);
                let (sf_pos, sf_vel) = starfield_propagate(&state, dt);
                let scale_r = sf_pos.norm();
                assert!(
                    (ours.pos - sf_pos).norm() <= 1e-9 * scale_r,
                    "hyperbolic position mismatch: e = {e}, q = {q}, dt = {dt}, \
                     |dr| = {:.3e} (scale {scale_r:.3e})",
                    (ours.pos - sf_pos).norm()
                );
                assert!(
                    (ours.vel - sf_vel).norm() <= 1e-9 * sf_vel.norm(),
                    "hyperbolic velocity mismatch: e = {e}, q = {q}, dt = {dt}"
                );
            }
        }
    }
}

// ===================================================================
// Work item 4: element round trips vs starfield::elementslib
// ===================================================================

/// Wrap an angle difference to (-pi, pi].
fn angle_diff(x: f64, y: f64) -> f64 {
    let mut d = (x - y) % TWO_PI;
    if d > std::f64::consts::PI {
        d -= TWO_PI;
    } else if d < -std::f64::consts::PI {
        d += TWO_PI;
    }
    d
}

/// Compare p9-core's cartesian_to_elements against starfield's
/// OsculatingElements for the same state vector.
///
/// Tolerances by regime:
/// - a (or q for near-parabolic), e, i: 1e-9 relative — pure algebra on
///   the state vector in both implementations.
/// - Ω, ω, M: 1e-8 rad away from degeneracies. Near-degenerate cases
///   (e ≤ 1e-6 makes ω/M individually ill-defined; i ≤ 1e-6 or i ≥ π−1e-6
///   makes Ω/ω individually ill-defined) are compared through the
///   well-defined combinations instead: ϖ = Ω ± ω (equatorial),
///   ω + ν ~ argument of latitude (circular), with 2e-5 rad slack since
///   the two codes resolve the degenerate angle from differently-rounded
///   tiny vectors.
fn check_elements_against_starfield(state: &StateVector) {
    let ours = cartesian_to_elements(state, GM_SUN);
    let sf =
        starfield::elementslib::osculating_elements_of(&state.pos, &state.vel, gm_sun_km3_s2());

    let e = sf.eccentricity();
    let i = sf.inclination_rad();

    // Shape: a for bound/hyperbolic, q always. Skip a near e = 1, where
    // a = -gm/(2*energy) is ill-conditioned (the ~1e-16 relative noise of
    // the energy is amplified by 1/|e-1|); q remains well-conditioned and
    // is always compared.
    if (e - 1.0).abs() > 1e-5 {
        let a_sf = sf.semi_major_axis_au();
        assert!(
            (ours.a - a_sf).abs() <= 1e-9 * a_sf.abs(),
            "a mismatch: {} vs {a_sf}",
            ours.a
        );
    }
    let q_sf = sf.periapsis_distance_au();
    let q_ours = ours.a * (1.0 - ours.e);
    assert!(
        (q_ours - q_sf).abs() <= 1e-8 * q_sf,
        "q mismatch: {q_ours} vs {q_sf} (e = {e})"
    );
    assert!(
        (ours.e - e).abs() <= 1e-9 * e.max(1.0),
        "e mismatch: {} vs {e}",
        ours.e
    );
    assert!((ours.i - i).abs() <= 1e-9, "i mismatch: {} vs {i}", ours.i);

    let near_circular = e <= 1e-6;
    let near_equatorial = !(1e-6..=std::f64::consts::PI - 1e-6).contains(&i);

    if !near_equatorial {
        assert!(
            angle_diff(ours.omega_big, sf.longitude_of_ascending_node_rad()).abs() <= 1e-8,
            "Omega mismatch: {} vs {}",
            ours.omega_big,
            sf.longitude_of_ascending_node_rad()
        );
    }
    if !near_circular && !near_equatorial {
        assert!(
            angle_diff(ours.omega, sf.argument_of_periapsis_rad()).abs() <= 1e-8,
            "omega mismatch: {} vs {}",
            ours.omega,
            sf.argument_of_periapsis_rad()
        );
    }
    if !near_circular && e < 1.0 {
        assert!(
            angle_diff(ours.mean_anomaly, sf.mean_anomaly_rad()).abs() <= 1e-8,
            "M mismatch: {} vs {} (e = {e})",
            ours.mean_anomaly,
            sf.mean_anomaly_rad()
        );
    }

    // Degenerate regimes: compare the non-degenerate combinations.
    // ϖ-like sum Ω + ω + M is well-defined for prograde orbits whenever the
    // state itself is (starfield folds the same convention into its
    // longitude accessors only for i < π/2, so restrict to prograde).
    if (near_circular || near_equatorial) && e < 1.0 && i < std::f64::consts::FRAC_PI_2 {
        let ours_l = ours.omega_big + ours.omega + ours.mean_anomaly;
        let sf_l = sf.longitude_of_ascending_node_rad()
            + sf.argument_of_periapsis_rad()
            + sf.mean_anomaly_rad();
        assert!(
            angle_diff(ours_l, sf_l).abs() <= 2e-5,
            "mean longitude mismatch in degenerate regime: {ours_l} vs {sf_l} \
             (e = {e}, i = {i})"
        );
    }
}

/// p9-core round trip: elements -> cartesian -> elements must reproduce the
/// state vector (the conversions are mutual inverses), checked through the
/// state because angle representations differ in degenerate regimes.
fn check_p9_round_trip(elem: &OrbitalElements) {
    let state = elements_to_cartesian(elem, GM_SUN);
    let back = cartesian_to_elements(&state, GM_SUN);
    let state2 = elements_to_cartesian(&back, GM_SUN);

    let scale_r = state.pos.norm();
    let scale_v = state.vel.norm();
    assert!(
        (state.pos - state2.pos).norm() <= 1e-9 * scale_r,
        "round-trip position drift: a = {}, e = {}, i = {}, O = {}, w = {}, M = {}: |dr| = {:.3e}",
        elem.a,
        elem.e,
        elem.i,
        elem.omega_big,
        elem.omega,
        elem.mean_anomaly,
        (state.pos - state2.pos).norm()
    );
    assert!(
        (state.vel - state2.vel).norm() <= 1e-9 * scale_v,
        "round-trip velocity drift: a = {}, e = {}, i = {}, O = {}, w = {}, M = {}",
        elem.a,
        elem.e,
        elem.i,
        elem.omega_big,
        elem.omega,
        elem.mean_anomaly
    );
}

#[test]
fn test_element_round_trips_vs_starfield() {
    let e_grid = [0.0, 1e-6, 0.3, 0.95, 0.999, 1.0 + 1e-6, 1.5];
    let i_grid_deg = [0.0, 1e-6_f64.to_degrees(), 30.0, 90.0, 150.0, 180.0];
    // All quadrants of Omega, omega, nu (via mean anomaly proxies).
    let angles = [0.5, 2.0, 3.8, 5.5];

    for &e in &e_grid {
        for &i_deg in &i_grid_deg {
            for &omega_big in &angles {
                for &omega in &angles {
                    for &m in &angles {
                        // Hyperbolic: a < 0 with q = 30 AU. Pick the anomaly
                        // through the true anomaly as a fraction of the
                        // asymptote angle so r stays O(q) — a raw M grid
                        // puts near-parabolic states (|a| = q/(e-1), huge)
                        // at absurd distances where every element is
                        // ill-conditioned. The fractions still cover
                        // inbound and outbound on both sides of perihelion.
                        let (a, mean_anomaly) = if e < 1.0 {
                            (450.0, m)
                        } else {
                            let nu_inf = (-1.0_f64 / e).acos();
                            let frac = 0.8 * (m - std::f64::consts::PI) / std::f64::consts::PI;
                            let nu = frac * nu_inf;
                            let h: f64 =
                                2.0 * (((e - 1.0) / (e + 1.0)).sqrt() * (nu / 2.0).tan()).atanh();
                            (-30.0 / (e - 1.0), e * h.sinh() - h)
                        };
                        let elem = OrbitalElements {
                            a,
                            e,
                            i: i_deg.to_radians(),
                            omega_big,
                            omega,
                            mean_anomaly,
                        };
                        let state = elements_to_cartesian(&elem, GM_SUN);
                        check_elements_against_starfield(&state);
                        check_p9_round_trip(&elem);
                    }
                }
            }
        }
    }
}
