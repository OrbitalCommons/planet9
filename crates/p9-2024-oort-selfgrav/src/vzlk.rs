//! Von Zeipel-Lidov-Kozai secular dynamics for IOC objects.
//!
//! The secular problem at fixed semi-major axis is one degree of freedom in
//! the canonically conjugate pair (omega, G) with
//!
//!   G = sqrt(GM_sun * a) * eta,   eta = sqrt(1 - e^2),
//!
//! at fixed vertical Delaunay momentum H_z = G cos i = sqrt(GM_sun a) * J_z
//! (conserved by the axisymmetry of the J2 + Miyamoto-Nagai terms). The
//! trajectory follows Hamilton's equations
//!
//!   d(omega)/dt = dH/dG,   dG/dt = -dH/d(omega),
//!
//! whose solutions are level sets of H(omega, G). The overall sign convention
//! of H only reverses the direction of motion along a level set; the level
//! sets themselves — and hence the perihelion excursion and the period
//! magnitude — are unaffected.

use serde::{Deserialize, Serialize};

use p9_core::constants::{GM_SUN, GYR_DAYS};
use p9_core::units::{au, days, radians, Angle, Length, Time};

use crate::hamiltonian::{secular_hamiltonian, HamiltonianParams};

/// A point on the vZLK phase portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VzlkPoint {
    /// Argument of perihelion (rad)
    pub omega: f64,
    /// Perihelion distance (AU)
    pub q: f64,
    /// Hamiltonian value
    pub h_value: f64,
}

impl VzlkPoint {
    /// Argument of perihelion as a typed [`Angle`].
    pub fn omega(&self) -> Angle {
        radians(self.omega)
    }
    /// Perihelion distance as a typed [`Length`].
    pub fn perihelion(&self) -> Length {
        au(self.q)
    }
}

/// A point along an integrated secular trajectory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrajectoryPoint {
    /// Time (days)
    pub t_days: f64,
    /// Argument of perihelion (rad)
    pub omega: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination (rad)
    pub i: f64,
    /// Perihelion distance (AU)
    pub q: f64,
    /// Hamiltonian value (conserved along the trajectory)
    pub h_value: f64,
}

impl TrajectoryPoint {
    /// Time as a typed [`Time`].
    pub fn time(&self) -> Time {
        days(self.t_days)
    }
    /// Argument of perihelion as a typed [`Angle`].
    pub fn omega(&self) -> Angle {
        radians(self.omega)
    }
    /// Inclination as a typed [`Angle`].
    pub fn inclination(&self) -> Angle {
        radians(self.i)
    }
    /// Perihelion distance as a typed [`Length`].
    pub fn perihelion(&self) -> Length {
        au(self.q)
    }
}

/// Configuration for vZLK phase portrait.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VzlkConfig {
    /// Fixed semi-major axis (AU)
    pub a: f64,
    /// Fixed vertical angular momentum J = sqrt(1-e^2)*cos(i)
    pub j_z: f64,
    /// Number of grid points in omega
    pub n_omega: usize,
    /// Number of grid points in q (perihelion)
    pub n_q: usize,
    /// Minimum perihelion to plot (AU)
    pub q_min: f64,
    /// Maximum perihelion to plot (AU)
    pub q_max: f64,
}

impl VzlkConfig {
    /// Fixed semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }
    /// Minimum perihelion as a typed [`Length`].
    pub fn q_min(&self) -> Length {
        au(self.q_min)
    }
    /// Maximum perihelion as a typed [`Length`].
    pub fn q_max(&self) -> Length {
        au(self.q_max)
    }

    /// Configuration for a typical IOC object at a=1000 AU.
    pub fn default_paper() -> Self {
        Self {
            a: 1000.0,
            j_z: 0.3, // Moderate inclination
            n_omega: 200,
            n_q: 100,
            q_min: 50.0,
            q_max: 500.0,
        }
    }
}

/// Generate a vZLK phase portrait on the (omega, q) plane.
///
/// At fixed (a, J_z), we scan omega and q, computing e and i from:
///   q = a*(1-e), J_z = sqrt(1-e^2)*cos(i)
pub fn compute_vzlk_portrait(config: &VzlkConfig, params: &HamiltonianParams) -> Vec<VzlkPoint> {
    let n_total = config.n_omega * config.n_q;
    let mut points = Vec::with_capacity(n_total);

    let domega = std::f64::consts::TAU / config.n_omega as f64;
    let dq = (config.q_max - config.q_min) / config.n_q as f64;

    for io in 0..config.n_omega {
        let omega = io as f64 * domega;
        for iq in 0..config.n_q {
            let q = config.q_min + iq as f64 * dq;
            let e = 1.0 - q / config.a;
            if e <= 0.0 || e >= 1.0 {
                continue;
            }

            let eta = (1.0 - e * e).sqrt();
            let cos_i = config.j_z / eta;
            if cos_i.abs() > 1.0 {
                continue;
            }
            let i = cos_i.acos();

            let h = secular_hamiltonian(config.a, e, i, omega, params);
            points.push(VzlkPoint {
                omega,
                q,
                h_value: h,
            });
        }
    }

    points
}

/// H expressed in the canonical momentum G at fixed (a, H_z, omega).
/// Returns None where (omega, G) is outside the physical region.
fn h_of_canonical(a: f64, h_z: f64, omega: f64, g: f64, params: &HamiltonianParams) -> Option<f64> {
    let l = (GM_SUN * a).sqrt();
    if g <= h_z.abs() || g >= l {
        return None;
    }
    let eta = g / l;
    let e = (1.0 - eta * eta).max(0.0).sqrt();
    let cos_i = h_z / g;
    if cos_i.abs() > 1.0 {
        return None;
    }
    Some(secular_hamiltonian(a, e, cos_i.acos(), omega, params))
}

/// Right-hand side of Hamilton's equations in (omega, G), via central
/// differences of H. `h_z` is the conserved vertical Delaunay momentum.
fn hamilton_rhs(a: f64, h_z: f64, omega: f64, g: f64, params: &HamiltonianParams) -> (f64, f64) {
    let l = (GM_SUN * a).sqrt();
    let dg = 1e-6 * l;
    let domega = 1e-6;

    let h_gp = h_of_canonical(a, h_z, omega, g + dg, params);
    let h_gm = h_of_canonical(a, h_z, omega, g - dg, params);
    let dh_dg = match (h_gp, h_gm) {
        (Some(p), Some(m)) => (p - m) / (2.0 * dg),
        _ => 0.0,
    };

    let h_op = h_of_canonical(a, h_z, omega + domega, g, params);
    let h_om = h_of_canonical(a, h_z, omega - domega, g, params);
    let dh_domega = match (h_op, h_om) {
        (Some(p), Some(m)) => (p - m) / (2.0 * domega),
        _ => 0.0,
    };

    (dh_dg, -dh_domega)
}

/// Integrate the secular equations of motion at fixed (a, J_z) from
/// (omega0, e0), using RK4 on the canonical pair (omega, G).
///
/// `j_z` is the dimensionless vertical momentum eta*cos(i); the initial
/// inclination follows from cos(i0) = j_z / eta0. Returns the trajectory
/// sampled at every step (n_steps + 1 points). H is recorded at each point;
/// its conservation is the primary invariant check (see tests).
pub fn integrate_vzlk(
    a: f64,
    j_z: f64,
    e0: f64,
    omega0: f64,
    params: &HamiltonianParams,
    dt_days: f64,
    n_steps: usize,
) -> Vec<TrajectoryPoint> {
    let l = (GM_SUN * a).sqrt();
    let h_z = l * j_z;

    let mut omega = omega0;
    let mut g = l * (1.0 - e0 * e0).sqrt();
    let mut out = Vec::with_capacity(n_steps + 1);

    let record = |t: f64, omega: f64, g: f64, out: &mut Vec<TrajectoryPoint>| {
        let eta = (g / l).clamp(0.0, 1.0);
        let e = (1.0 - eta * eta).max(0.0).sqrt();
        let cos_i = (h_z / g).clamp(-1.0, 1.0);
        let i = cos_i.acos();
        let h = secular_hamiltonian(a, e, i, omega, params);
        out.push(TrajectoryPoint {
            t_days: t,
            omega,
            e,
            i,
            q: a * (1.0 - e),
            h_value: h,
        });
    };

    record(0.0, omega, g, &mut out);

    for step in 0..n_steps {
        let (k1w, k1g) = hamilton_rhs(a, h_z, omega, g, params);
        let (k2w, k2g) = hamilton_rhs(
            a,
            h_z,
            omega + 0.5 * dt_days * k1w,
            g + 0.5 * dt_days * k1g,
            params,
        );
        let (k3w, k3g) = hamilton_rhs(
            a,
            h_z,
            omega + 0.5 * dt_days * k2w,
            g + 0.5 * dt_days * k2g,
            params,
        );
        let (k4w, k4g) = hamilton_rhs(a, h_z, omega + dt_days * k3w, g + dt_days * k3g, params);

        omega += dt_days / 6.0 * (k1w + 2.0 * k2w + 2.0 * k3w + k4w);
        g += dt_days / 6.0 * (k1g + 2.0 * k2g + 2.0 * k3g + k4g);

        // Keep G inside the physical strip (|H_z|, L).
        g = g.clamp(h_z.abs() * (1.0 + 1e-12) + 1e-300, l * (1.0 - 1e-12));

        record((step + 1) as f64 * dt_days, omega, g, &mut out);
    }

    out
}

/// Secular evolutionary timescale (Gyr): tau = 2*pi / |d(omega)/dt| with
/// d(omega)/dt = dH/dG evaluated at fixed H_z. Since G = sqrt(GM_sun a) eta,
/// this is the dimensionally consistent form of the chain-rule estimate
/// tau ~ 2*pi*sqrt(GM_sun a) / |dH/d(eta)| (H in AU^2/day^2, G in AU^2/day,
/// so dH/dG is a frequency in 1/day). This replaces the previous version
/// that divided dH/de — not a frequency — into 2*pi.
pub fn evolutionary_timescale_gyr(
    a: f64,
    e: f64,
    i: f64,
    omega: f64,
    params: &HamiltonianParams,
) -> f64 {
    let l = (GM_SUN * a).sqrt();
    let eta = (1.0 - e * e).sqrt();
    let g = l * eta;
    let h_z = g * i.cos();

    let (dh_dg, _) = hamilton_rhs(a, h_z, omega, g, params);
    if dh_dg.abs() < 1e-300 {
        return f64::INFINITY;
    }
    let period_days = std::f64::consts::TAU / dh_dg.abs();
    period_days / GYR_DAYS
}

/// Minimum perihelion reached along the secular trajectory (level set of H)
/// through (e0, omega0) at fixed (a, J_z).
///
/// Integrates Hamilton's equations for `n_cycles` estimated secular periods
/// and returns the smallest q = a(1 - e) encountered. The hard kinematic
/// floor is q >= a(1 - sqrt(1 - J_z^2)); the dynamical minimum returned here
/// is generally well above it because H conservation confines the motion.
pub fn minimum_perihelion(
    a: f64,
    j_z: f64,
    e0: f64,
    omega0: f64,
    params: &HamiltonianParams,
) -> f64 {
    let eta0 = (1.0 - e0 * e0).sqrt();
    let cos_i0 = (j_z / eta0).clamp(-1.0, 1.0);
    let tau_gyr = evolutionary_timescale_gyr(a, e0, cos_i0.acos(), omega0, params);

    // Two estimated secular cycles, 2000 steps per cycle. The timescale can
    // be enormous (>> Hubble time); that only sets the step in days.
    let n_steps = 4000;
    let dt_days = if tau_gyr.is_finite() {
        2.0 * tau_gyr * GYR_DAYS / n_steps as f64
    } else {
        // Static point: any step works, trajectory stays put.
        1e9
    };

    let traj = integrate_vzlk(a, j_z, e0, omega0, params, dt_days, n_steps);
    traj.iter().map(|p| p.q).fold(f64::INFINITY, f64::min)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamiltonian::HamiltonianParams;
    use p9_core::constants::DEG2RAD;

    fn cheap_params() -> HamiltonianParams {
        HamiltonianParams {
            n_quadrature: 64,
            ..HamiltonianParams::default_paper()
        }
    }

    #[test]
    fn typed_accessors_match_f64() {
        use approx::assert_relative_eq;
        use uom::si::angle::radian;
        use uom::si::length::astronomical_unit;
        use uom::si::time::day;
        let cfg = VzlkConfig::default_paper();
        assert_relative_eq!(
            cfg.semi_major_axis().get::<astronomical_unit>(),
            cfg.a,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            cfg.q_min().get::<astronomical_unit>(),
            cfg.q_min,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            cfg.q_max().get::<astronomical_unit>(),
            cfg.q_max,
            epsilon = 1e-12
        );
        let p = VzlkPoint {
            omega: 1.2,
            q: 200.0,
            h_value: 0.0,
        };
        assert_relative_eq!(p.omega().get::<radian>(), p.omega, epsilon = 1e-12);
        assert_relative_eq!(
            p.perihelion().get::<astronomical_unit>(),
            p.q,
            epsilon = 1e-12
        );
        let tp = TrajectoryPoint {
            t_days: 1.0e6,
            omega: 0.5,
            e: 0.9,
            i: 0.7,
            q: 150.0,
            h_value: 0.0,
        };
        assert_relative_eq!(tp.time().get::<day>(), tp.t_days, max_relative = 1e-12);
        assert_relative_eq!(tp.omega().get::<radian>(), tp.omega, epsilon = 1e-12);
        assert_relative_eq!(tp.inclination().get::<radian>(), tp.i, epsilon = 1e-12);
        assert_relative_eq!(
            tp.perihelion().get::<astronomical_unit>(),
            tp.q,
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_vzlk_portrait_nonempty() {
        let config = VzlkConfig {
            n_omega: 40,
            n_q: 20,
            ..VzlkConfig::default_paper()
        };
        let params = cheap_params();
        let points = compute_vzlk_portrait(&config, &params);
        assert!(!points.is_empty(), "portrait should have points");
        for p in &points {
            assert!(p.h_value.is_finite());
        }
    }

    #[test]
    fn test_hamiltonian_conserved_along_trajectory() {
        // Primary invariant check: H must be constant along an integrated
        // secular trajectory (it generates the flow).
        let params = cheap_params();
        let (a, j_z, e0, omega0) = (1000.0, 0.3, 0.85, 80.0 * DEG2RAD);
        let eta0 = (1.0_f64 - e0 * e0).sqrt();
        let tau = evolutionary_timescale_gyr(a, e0, (j_z / eta0).acos(), omega0, &params);
        let n_steps = 800;
        let dt = tau * p9_core::constants::GYR_DAYS / n_steps as f64;

        let traj = integrate_vzlk(a, j_z, e0, omega0, &params, dt, n_steps);
        let h0 = traj[0].h_value;
        let h_scale = h0.abs();
        let max_dev = traj
            .iter()
            .map(|p| (p.h_value - h0).abs())
            .fold(0.0, f64::max);
        assert!(
            max_dev / h_scale < 1e-6,
            "relative H drift = {:.3e}",
            max_dev / h_scale
        );
    }

    #[test]
    fn test_trajectory_respects_jz_bound() {
        let params = cheap_params();
        let j_z = 0.3;
        let traj = integrate_vzlk(1000.0, j_z, 0.85, 80.0 * DEG2RAD, &params, 1e10, 500);
        for p in &traj {
            let eta = (1.0 - p.e * p.e).sqrt();
            assert!(
                eta >= j_z.abs() - 1e-9,
                "eta = {eta} violates kinematic bound at e = {}",
                p.e
            );
        }
    }

    #[test]
    fn test_minimum_perihelion_bounds() {
        let params = cheap_params();
        let (a, j_z, e0) = (1000.0, 0.3, 0.5);
        let q_min = minimum_perihelion(a, j_z, e0, 90.0 * DEG2RAD, &params);
        // Kinematic floor: q >= a(1 - sqrt(1 - J_z^2)).
        let q_floor = a * (1.0 - (1.0 - j_z * j_z).sqrt());
        assert!(
            q_min >= q_floor - 1.0,
            "q_min = {q_min} below floor {q_floor}"
        );
        // Cannot exceed the initial perihelion by definition of a minimum
        // along a trajectory through (e0, omega0).
        assert!(q_min <= a * (1.0 - e0) + 1e-6, "q_min = {q_min}");
    }

    #[test]
    fn test_evolutionary_timescale_exceeds_solar_system_age() {
        // The paper's central conclusion: for nominal IOC parameters
        // (M_IOC = 3 M_earth, a = 1000 AU, q = 100 AU, i = 30 deg) the
        // vZLK evolutionary timescale far exceeds 4.5 Gyr.
        let params = cheap_params();
        let tau = evolutionary_timescale_gyr(1000.0, 0.9, 30.0 * DEG2RAD, 0.0, &params);
        assert!(
            tau > 4.5,
            "tau = {tau} Gyr should exceed the solar system age"
        );
    }
}
