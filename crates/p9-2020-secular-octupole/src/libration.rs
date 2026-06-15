//! The headline result: including the octupole term produces an apsidal
//! (`dvarpi`) LIBRATION island where the quadrupole-only theory only
//! CIRCULATES.
//!
//! We integrate the one-degree-of-freedom secular equations of motion in the
//! canonical pair (Gamma, dvarpi),
//!
//!   d(dvarpi)/dt = + dH/dGamma,
//!   dGamma/dt    = - dH/d(dvarpi),
//!
//! with Gamma = sqrt(GM_sun a)(1 - sqrt(1-e^2)) and `a` constant (secular
//! theory conserves semi-major axis). We evolve a test orbit and classify the
//! trajectory:
//!
//! - **Circulation**: `dvarpi` advances monotonically through a full 2*pi
//!   (the apsidal line sweeps all orientations).
//! - **Libration**: `dvarpi` oscillates within a bounded arc and reverses;
//!   the orbit is APSIDALLY CONFINED.
//!
//! At quadrupole order H = H(e) only, so dGamma/dt = -dH/d(dvarpi) = 0: Gamma
//! (hence e) is frozen and dvarpi advances at a constant rate -> pure
//! circulation, NEVER libration. The octupole term supplies the
//! d(dvarpi)-dependence that lets a libration island exist.
//!
//! # Fixed point and the aligned/anti-aligned sign
//!
//! The coplanar octupole has an elliptic fixed point (the libration center) at
//! dvarpi = 0 with e* = epsilon_oct/3, found by setting dGamma/dt = 0
//! (sin dvarpi = 0) and dvarpi/dt = 0 (dH/de = 0): the negative quadrupole
//! slope dH_quad/de = -3 C_quad e balances the positive octupole slope
//! C_oct e_p cos dvarpi only at dvarpi = 0. So the coplanar libration center is
//! APSIDALLY ALIGNED with the perturber. The observed ETNOs are
//! ANTI-aligned with Planet Nine; that sign comes from the inclined / resonant
//! dynamics (and the perihelion-detached, high-i regime), not from the
//! coplanar secular fixed point reproduced here. This crate reproduces the
//! octupole-driven libration MECHANISM and the e*, epsilon_oct scalings; the
//! observed anti-alignment sign is documented as out of reach of the coplanar
//! reduction (cf. the analogous note in p9-2019-selfgrav-disk).
//!
//! Reference: Kohne & Batygin (2020).

use crate::octupole::{hamiltonian_octupole, hamiltonian_quadrupole};
use p9_core::constants::GM_SUN;
use std::f64::consts::PI;

/// Eccentricity from the canonical action Gamma at fixed a.
fn e_of_gamma(a: f64, gamma: f64) -> f64 {
    // Gamma = sqrt(GM a)(1 - sqrt(1-e^2)) => sqrt(1-e^2) = 1 - Gamma/sqrt(GM a)
    let s = 1.0 - gamma / (GM_SUN * a).sqrt();
    let s = s.clamp(0.0, 1.0);
    (1.0 - s * s).sqrt()
}

/// Canonical action Gamma from eccentricity at fixed a.
pub fn gamma_of_e(a: f64, e: f64) -> f64 {
    (GM_SUN * a).sqrt() * (1.0 - (1.0 - e * e).sqrt())
}

/// Outcome of classifying a secular trajectory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ApsidalMotion {
    /// dvarpi sweeps a full circle: free precession.
    Circulation,
    /// dvarpi oscillates within a bounded arc: apsidal confinement.
    Libration,
}

/// Whether the Hamiltonian includes the octupole term.
#[derive(Debug, Clone, Copy)]
pub enum Order {
    Quadrupole,
    Octupole,
}

fn hamiltonian(order: Order, a: f64, e: f64, dvarpi: f64, a_p: f64, e_p: f64, gm_p: f64) -> f64 {
    match order {
        Order::Quadrupole => hamiltonian_quadrupole(a, e, a_p, e_p, gm_p),
        Order::Octupole => hamiltonian_octupole(a, e, dvarpi, a_p, e_p, gm_p),
    }
}

/// Numerical partials of H w.r.t. (Gamma, dvarpi) by central differences in
/// (e, dvarpi), with dGamma -> de via de/dGamma.
fn derivs(order: Order, a: f64, e: f64, dvarpi: f64, a_p: f64, e_p: f64, gm_p: f64) -> (f64, f64) {
    let de = 1e-6;
    let dd = 1e-6;
    let dh_de = (hamiltonian(order, a, e + de, dvarpi, a_p, e_p, gm_p)
        - hamiltonian(order, a, e - de, dvarpi, a_p, e_p, gm_p))
        / (2.0 * de);
    let dh_dd = (hamiltonian(order, a, e, dvarpi + dd, a_p, e_p, gm_p)
        - hamiltonian(order, a, e, dvarpi - dd, a_p, e_p, gm_p))
        / (2.0 * dd);
    // dGamma/de = sqrt(GM a) e / sqrt(1-e^2)  =>  dH/dGamma = dH/de * de/dGamma.
    let dgamma_de = (GM_SUN * a).sqrt() * e / (1.0 - e * e).sqrt();
    let dh_dgamma = dh_de / dgamma_de;
    (dh_dgamma, dh_dd)
}

/// Integrate the secular EOM from initial (e0, dvarpi0) and classify the
/// apsidal motion. `n_steps` * `dt` should cover at least one secular period;
/// we adapt by running until dvarpi either completes a full turn (circulation)
/// or reverses direction twice about a center (libration), or the step budget
/// is exhausted.
///
/// Returns the classification and the (min, max) of dvarpi visited (radians,
/// unwrapped relative to the start).
#[allow(clippy::too_many_arguments)]
pub fn classify_trajectory(
    order: Order,
    a: f64,
    e0: f64,
    dvarpi0: f64,
    a_p: f64,
    e_p: f64,
    gm_p: f64,
    dt: f64,
    n_steps: usize,
) -> (ApsidalMotion, f64, f64) {
    let mut gamma = gamma_of_e(a, e0);
    let mut dvarpi = dvarpi0;
    let mut dvarpi_unwrapped = dvarpi0;
    let mut prev = dvarpi0;

    let mut total_advance = 0.0_f64; // net unwrapped change
    let mut sign_changes = 0u32;
    let mut last_rate_sign = 0.0_f64;
    let mut min_dv = dvarpi0;
    let mut max_dv = dvarpi0;

    for _ in 0..n_steps {
        let e = e_of_gamma(a, gamma);
        if !(1e-6..=0.999).contains(&e) {
            break;
        }
        // RK4 in (Gamma, dvarpi): dGamma/dt = -dH/ddvarpi, ddvarpi/dt = dH/dGamma.
        let f = |g: f64, d: f64| {
            let ee = e_of_gamma(a, g);
            let (dh_dg, dh_dd) = derivs(order, a, ee, d, a_p, e_p, gm_p);
            (-dh_dd, dh_dg) // (dGamma/dt, ddvarpi/dt)
        };
        let (k1g, k1d) = f(gamma, dvarpi);
        let (k2g, k2d) = f(gamma + 0.5 * dt * k1g, dvarpi + 0.5 * dt * k1d);
        let (k3g, k3d) = f(gamma + 0.5 * dt * k2g, dvarpi + 0.5 * dt * k2d);
        let (k4g, k4d) = f(gamma + dt * k3g, dvarpi + dt * k3d);
        gamma += dt / 6.0 * (k1g + 2.0 * k2g + 2.0 * k3g + k4g);
        let ddv = dt / 6.0 * (k1d + 2.0 * k2d + 2.0 * k3d + k4d);
        dvarpi += ddv;
        dvarpi_unwrapped += ddv;

        let rate_sign = ddv.signum();
        if last_rate_sign != 0.0 && rate_sign != last_rate_sign {
            sign_changes += 1;
        }
        last_rate_sign = rate_sign;

        total_advance += dvarpi_unwrapped - prev;
        prev = dvarpi_unwrapped;
        min_dv = min_dv.min(dvarpi_unwrapped);
        max_dv = max_dv.max(dvarpi_unwrapped);

        // Circulation: net monotone advance exceeds a full turn.
        if total_advance.abs() > 2.0 * PI {
            return (ApsidalMotion::Circulation, min_dv, max_dv);
        }
        // Libration: the precession rate has reversed (turned around) at least
        // twice (out and back) within a bounded arc.
        if sign_changes >= 2 && (max_dv - min_dv) < 2.0 * PI {
            return (ApsidalMotion::Libration, min_dv, max_dv);
        }
    }
    // Step budget exhausted without a full turn and with a bounded excursion:
    // treat a clearly sub-2pi span with a reversal as libration, else
    // circulation (it was advancing but slowly).
    if sign_changes >= 1 && (max_dv - min_dv) < 2.0 * PI {
        (ApsidalMotion::Libration, min_dv, max_dv)
    } else {
        (ApsidalMotion::Circulation, min_dv, max_dv)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN};

    fn gm_p9(mass_earth: f64) -> f64 {
        mass_earth * EARTH_MASS_SOLAR * GM_SUN
    }

    #[test]
    fn gamma_e_roundtrip() {
        let a = 250.0;
        for &e in &[0.1, 0.3, 0.5, 0.7, 0.9] {
            let g = gamma_of_e(a, e);
            let e2 = e_of_gamma(a, g);
            assert!((e - e2).abs() < 1e-10, "e {e} -> {e2}");
        }
    }

    #[test]
    fn quadrupole_only_circulates() {
        // Quadrupole H = H(e): dGamma/dt = 0, dvarpi advances forever.
        // Whatever the apsidal start, the classification is circulation.
        let gm = gm_p9(10.0);
        let (a, a_p, e_p) = (250.0, 700.0, 0.6);
        // dt sized to a fraction of the secular period.
        let dt = 2.0e6;
        for &dv0 in &[0.0, 0.5, PI, 2.5] {
            let (motion, _, _) =
                classify_trajectory(Order::Quadrupole, a, 0.4, dv0, a_p, e_p, gm, dt, 200_000);
            assert_eq!(
                motion,
                ApsidalMotion::Circulation,
                "quadrupole should circulate at dv0={dv0}"
            );
        }
    }

    #[test]
    fn octupole_produces_libration_island() {
        // The headline result: with the octupole on, an orbit started near the
        // secular fixed point is trapped in a libration island, while the same
        // orbit circulates under quadrupole-only.
        //
        // The coplanar fixed point sits at dvarpi = 0 (the apse aligned with
        // the perturber) and e* = epsilon_oct/3 (where dH_quad/de balances the
        // octupole term). See the module note on the aligned-vs-anti-aligned
        // sign caveat. We start slightly off-center so the island has a finite
        // libration amplitude.
        let gm = gm_p9(10.0);
        let (a, a_p, e_p) = (250.0, 700.0, 0.6);
        let e0 = 0.3; // below the e* ~ 0.47 fixed point -> librates
        let dt = 2.0e6;

        let (quad_motion, _, _) =
            classify_trajectory(Order::Quadrupole, a, e0, 0.3, a_p, e_p, gm, dt, 400_000);
        let (oct_motion, lo, hi) =
            classify_trajectory(Order::Octupole, a, e0, 0.3, a_p, e_p, gm, dt, 400_000);

        assert_eq!(quad_motion, ApsidalMotion::Circulation);
        assert_eq!(
            oct_motion,
            ApsidalMotion::Libration,
            "octupole should librate; span {lo}..{hi}"
        );
        // The librating arc must be strictly less than a full turn.
        assert!((hi - lo) < 2.0 * PI, "libration span {} >= 2pi", hi - lo);
    }

    #[test]
    fn libration_center_is_aligned_apse() {
        // An orbit started exactly at the fixed point (dvarpi = 0,
        // e* = epsilon_oct/3) barely moves in dvarpi: it IS the libration
        // center. Confirms the coplanar octupole fixed point location.
        let gm = gm_p9(10.0);
        let (a, a_p, e_p) = (250.0, 700.0, 0.6);
        let eps = 4.0 * crate::octupole::K_OCT * (a / a_p) * e_p / (1.0 - e_p * e_p);
        let e_star = eps / 3.0;
        let dt = 1.0e6;
        let (motion, lo, hi) =
            classify_trajectory(Order::Octupole, a, e_star, 0.0, a_p, e_p, gm, dt, 2_000_000);
        assert_eq!(motion, ApsidalMotion::Libration);
        // Excursion around the center is tiny.
        assert!(
            (hi - lo) < 0.1,
            "fixed-point excursion {} too large",
            hi - lo
        );
    }

    #[test]
    fn quadrupole_circulates_from_same_island_seed() {
        // The very initial condition that LIBRATES under octupole must
        // CIRCULATE under quadrupole-only -- the contrast is the whole point.
        let gm = gm_p9(10.0);
        let (a, a_p, e_p) = (250.0, 700.0, 0.6);
        let dt = 2.0e6;
        let (motion, _, _) =
            classify_trajectory(Order::Quadrupole, a, 0.3, 0.6, a_p, e_p, gm, dt, 400_000);
        assert_eq!(motion, ApsidalMotion::Circulation);
    }
}
