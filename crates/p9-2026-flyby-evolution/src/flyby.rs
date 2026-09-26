//! The Pfalzner, Govind & Portegies Zwart (2024) flyby (their model A, the
//! starting point of the 2026 paper): a 0.8 M☉ star on a parabolic orbit
//! with periastron 110 AU, inclination 70° and argument of periastron 80°
//! passes a thin, initially circular disc of massless tracers of constant
//! surface density out to 150 AU (model A1).
//!
//! The encounter is the restricted three-body problem Sun + perturber +
//! tracer. The Sun–perturber relative motion is the exact parabola (Barker's
//! equation); tracers are integrated in the Sun–perturber barycentric frame
//! with an adaptive Dormand–Prince 5(4) scheme from the perturber's inbound
//! [`FlybyConfig::start_distance_au`] to the same distance outbound, then
//! converted back to heliocentric elements.

use nalgebra::{Matrix3, Vector3};
use p9_core::constants::{DEG2RAD, GM_SUN, TWO_PI};
use p9_core::types::{cartesian_to_elements, OrbitalElements, StateVector};
use rand::{Rng, SeedableRng};
use rayon::prelude::*;

/// Perturber orbit relative to the Sun.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FlybyParams {
    /// Perturber mass (solar masses).
    pub mass_solar: f64,
    /// Periastron distance (AU).
    pub q_au: f64,
    /// Inclination of the perturber orbit to the disc plane (deg).
    pub inclination_deg: f64,
    /// Argument of periastron (deg).
    pub arg_periastron_deg: f64,
    /// Longitude of the ascending node (deg); irrelevant for an
    /// axisymmetric disc and set to zero.
    pub node_deg: f64,
}

impl FlybyParams {
    /// Pfalzner et al. (2024) model A: 0.8 M☉, 110 AU, i = 70°, ω = 80°.
    pub fn model_a() -> Self {
        Self {
            mass_solar: 0.8,
            q_au: 110.0,
            inclination_deg: 70.0,
            arg_periastron_deg: 80.0,
            node_deg: 0.0,
        }
    }

    /// Gravitational parameter of the Sun–perturber relative orbit
    /// (AU³/day²).
    pub fn mu(&self) -> f64 {
        GM_SUN * (1.0 + self.mass_solar)
    }

    /// Perifocal → disc-frame rotation R_z(Ω) R_x(i) R_z(ω).
    fn rotation(&self) -> Matrix3<f64> {
        let rz = |t: f64| {
            let (s, c) = t.sin_cos();
            Matrix3::new(c, -s, 0.0, s, c, 0.0, 0.0, 0.0, 1.0)
        };
        let rx = |t: f64| {
            let (s, c) = t.sin_cos();
            Matrix3::new(1.0, 0.0, 0.0, 0.0, c, -s, 0.0, s, c)
        };
        rz(self.node_deg * DEG2RAD)
            * rx(self.inclination_deg * DEG2RAD)
            * rz(self.arg_periastron_deg * DEG2RAD)
    }

    /// Barker's equation: the parabolic anomaly D = tan(ν/2) at time `t`
    /// (days) from periastron, from D + D³/3 = t·√(μ/2q³).
    pub fn parabolic_anomaly(&self, t_days: f64) -> f64 {
        let tau = t_days * (self.mu() / (2.0 * self.q_au.powi(3))).sqrt();
        let a = 1.5 * tau;
        let b = (a + (a * a + 1.0).sqrt()).cbrt();
        b - 1.0 / b
    }

    /// Perturber position and velocity relative to the Sun at `t_days` from
    /// periastron (AU, AU/day), disc frame.
    pub fn relative_state(&self, t_days: f64) -> StateVector {
        let d = self.parabolic_anomaly(t_days);
        let q = self.q_au;
        let pos = Vector3::new(q * (1.0 - d * d), 2.0 * q * d, 0.0);
        let one_p = 1.0 + d * d;
        let (sin_nu, cos_nu) = (2.0 * d / one_p, (1.0 - d * d) / one_p);
        let vs = (self.mu() / (2.0 * q)).sqrt();
        let vel = Vector3::new(-vs * sin_nu, vs * (1.0 + cos_nu), 0.0);
        let r = self.rotation();
        StateVector::new(r * pos, r * vel)
    }

    /// Time (days) after periastron at which the perturber is `r_au` from
    /// the Sun on the outbound leg (the inbound time is its negative).
    pub fn time_at_distance(&self, r_au: f64) -> f64 {
        assert!(r_au >= self.q_au);
        let d = (r_au / self.q_au - 1.0).sqrt();
        (d + d * d * d / 3.0) / (self.mu() / (2.0 * self.q_au.powi(3))).sqrt()
    }

    /// Sun and perturber barycentric states at `t_days`.
    pub fn barycentric_states(&self, t_days: f64) -> (StateVector, StateVector) {
        let rel = self.relative_state(t_days);
        let m = self.mass_solar;
        let f_sun = -m / (1.0 + m);
        let f_star = 1.0 / (1.0 + m);
        (
            StateVector::new(f_sun * rel.pos, f_sun * rel.vel),
            StateVector::new(f_star * rel.pos, f_star * rel.vel),
        )
    }
}

/// Pre-flyby disc.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct DiscParams {
    /// Inner edge (AU). Not stated by the papers; the flyby leaves the disc
    /// undisturbed inside ~33 AU (their r_d = 0.28 M_p^{-0.32} r_peri), so
    /// the tracers start at the inner edge of the region that matters.
    pub r_min_au: f64,
    /// Outer edge (AU): 150 AU for model A1, 300 AU for A2.
    pub r_max_au: f64,
    pub n_particles: usize,
}

impl DiscParams {
    /// Model A1: constant surface density from 30 to 150 AU.
    pub fn model_a1(n_particles: usize) -> Self {
        Self {
            r_min_au: 30.0,
            r_max_au: 150.0,
            n_particles,
        }
    }
}

/// Full encounter set-up.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct FlybyConfig {
    pub flyby: FlybyParams,
    pub disc: DiscParams,
    /// Sun–perturber distance at which the integration starts and ends
    /// (AU). At 2000 AU the perturber's tidal acceleration on the outer disc
    /// edge is below 10⁻³ of the solar attraction.
    pub start_distance_au: f64,
    /// Relative tolerance of the adaptive integrator.
    pub rel_tol: f64,
    pub seed: u64,
}

impl FlybyConfig {
    /// Model A / disc A1 with `n` tracers.
    pub fn model_a1(n: usize) -> Self {
        Self {
            flyby: FlybyParams::model_a(),
            disc: DiscParams::model_a1(n),
            start_distance_au: 2000.0,
            rel_tol: 1e-10,
            seed: 20260903, // arXiv:2609.03575 posting date
        }
    }
}

/// Post-encounter state of every tracer.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct FlybyOutcome {
    pub config: FlybyConfig,
    /// Pre-flyby heliocentric radius of each tracer (AU).
    pub initial_radii: Vec<f64>,
    /// Post-flyby heliocentric elements; `None` if unbound from the Sun.
    pub elements: Vec<Option<OrbitalElements>>,
}

impl FlybyOutcome {
    /// Tracers still bound to the Sun, as `(initial radius, elements)`.
    pub fn bound(&self) -> Vec<(f64, OrbitalElements)> {
        self.initial_radii
            .iter()
            .zip(&self.elements)
            .filter_map(|(&r, e)| e.map(|e| (r, e)))
            .collect()
    }

    pub fn n_bound(&self) -> usize {
        self.elements.iter().filter(|e| e.is_some()).count()
    }

    /// Σ ∝ 1/r-weighted fraction of the disc lost from the Sun (captured by
    /// the perturber or ejected).
    pub fn unbound_fraction(&self) -> f64 {
        let (mut lost, mut total) = (0.0, 0.0);
        for (&r, e) in self.initial_radii.iter().zip(&self.elements) {
            let w = crate::groups::surface_density_weight(r);
            total += w;
            if e.is_none() {
                lost += w;
            }
        }
        lost / total
    }

    /// Weighted Table 1 census of the bound tracers.
    pub fn census(&self) -> crate::groups::Census {
        crate::groups::Census::of(self.bound().iter().map(|(r, e)| (*r, e)))
    }
}

/// Integrate the encounter for every tracer (parallel).
pub fn run_flyby(config: &FlybyConfig) -> FlybyOutcome {
    let mut rng = rand::rngs::StdRng::seed_from_u64(config.seed);
    let d = config.disc;
    let tracers: Vec<(f64, f64)> = (0..d.n_particles)
        .map(|_| {
            let u: f64 = rng.gen();
            let r = (d.r_min_au.powi(2) + u * (d.r_max_au.powi(2) - d.r_min_au.powi(2))).sqrt();
            (r, rng.gen_range(0.0..TWO_PI))
        })
        .collect();
    let t_end = config.flyby.time_at_distance(config.start_distance_au);
    let t_start = -t_end;
    let (sun0, _) = config.flyby.barycentric_states(t_start);
    let (sun1, _) = config.flyby.barycentric_states(t_end);

    let elements: Vec<Option<OrbitalElements>> = tracers
        .par_iter()
        .map(|&(r, phi)| {
            let v = (GM_SUN / r).sqrt();
            let pos = Vector3::new(r * phi.cos(), r * phi.sin(), 0.0) + sun0.pos;
            let vel = Vector3::new(-v * phi.sin(), v * phi.cos(), 0.0) + sun0.vel;
            let (p, v) = integrate_tracer(&config.flyby, pos, vel, t_start, t_end, config.rel_tol);
            let helio = StateVector::new(p - sun1.pos, v - sun1.vel);
            let e = cartesian_to_elements(&helio, GM_SUN);
            (e.e < 1.0 && e.a > 0.0).then_some(e)
        })
        .collect();
    FlybyOutcome {
        config: *config,
        initial_radii: tracers.iter().map(|t| t.0).collect(),
        elements,
    }
}

/// Barycentric acceleration on a tracer.
fn acceleration(flyby: &FlybyParams, t: f64, pos: &Vector3<f64>) -> Vector3<f64> {
    let (sun, star) = flyby.barycentric_states(t);
    let ds = pos - sun.pos;
    let dp = pos - star.pos;
    -GM_SUN * ds / ds.norm().powi(3) - GM_SUN * flyby.mass_solar * dp / dp.norm().powi(3)
}

/// Dormand–Prince 5(4) with step-size control on the 6-vector (pos, vel).
fn integrate_tracer(
    flyby: &FlybyParams,
    mut pos: Vector3<f64>,
    mut vel: Vector3<f64>,
    t0: f64,
    t1: f64,
    rel_tol: f64,
) -> (Vector3<f64>, Vector3<f64>) {
    const A: [[f64; 6]; 6] = [
        [1.0 / 5.0, 0.0, 0.0, 0.0, 0.0, 0.0],
        [3.0 / 40.0, 9.0 / 40.0, 0.0, 0.0, 0.0, 0.0],
        [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0, 0.0, 0.0, 0.0],
        [
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
            0.0,
            0.0,
        ],
        [
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
            0.0,
        ],
        [
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
        ],
    ];
    const C: [f64; 6] = [1.0 / 5.0, 3.0 / 10.0, 4.0 / 5.0, 8.0 / 9.0, 1.0, 1.0];
    // 5th − 4th order weights (error estimate).
    const E: [f64; 7] = [
        71.0 / 57600.0,
        0.0,
        -71.0 / 16695.0,
        71.0 / 1920.0,
        -17253.0 / 339200.0,
        22.0 / 525.0,
        -1.0 / 40.0,
    ];
    let abs_tol_pos = 1e-9;
    let abs_tol_vel = 1e-12;

    let mut t = t0;
    let mut h = 0.005 * TWO_PI * (pos.norm().powi(3) / GM_SUN).sqrt();
    let mut k = [(Vector3::zeros(), Vector3::zeros()); 7];
    while t < t1 {
        if t + h > t1 {
            h = t1 - t;
        }
        k[0] = (vel, acceleration(flyby, t, &pos));
        for s in 0..6 {
            let mut p = pos;
            let mut v = vel;
            for (j, &a) in A[s].iter().enumerate().take(s + 1) {
                if a != 0.0 {
                    p += h * a * k[j].0;
                    v += h * a * k[j].1;
                }
            }
            k[s + 1] = (v, acceleration(flyby, t + C[s] * h, &p));
        }
        // k[6] was evaluated at the 5th-order solution (FSAL): its stage
        // inputs are the new state.
        let mut p_new = pos;
        let mut v_new = vel;
        for (j, &a) in A[5].iter().enumerate() {
            if a != 0.0 {
                p_new += h * a * k[j].0;
                v_new += h * a * k[j].1;
            }
        }
        let mut err_p = Vector3::zeros();
        let mut err_v = Vector3::zeros();
        for (j, &e) in E.iter().enumerate() {
            err_p += h * e * k[j].0;
            err_v += h * e * k[j].1;
        }
        let scale_p = abs_tol_pos + rel_tol * pos.norm().max(p_new.norm());
        let scale_v = abs_tol_vel + rel_tol * vel.norm().max(v_new.norm());
        let err = (err_p.norm() / scale_p).max(err_v.norm() / scale_v);
        if err <= 1.0 {
            t += h;
            pos = p_new;
            vel = v_new;
        }
        let factor = if err > 0.0 {
            (0.9 * err.powf(-0.2)).clamp(0.2, 5.0)
        } else {
            5.0
        };
        h *= factor;
    }
    (pos, vel)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::YEAR_DAYS;

    #[test]
    fn barker_periastron_and_energy() {
        let f = FlybyParams::model_a();
        let peri = f.relative_state(0.0);
        assert!((peri.pos.norm() - f.q_au).abs() < 1e-9);
        // Zero specific energy at every epoch.
        for t in [-5000.0, -100.0, 0.0, 250.0, 9000.0] {
            let s = f.relative_state(t * YEAR_DAYS);
            let energy = 0.5 * s.vel.norm_squared() - f.mu() / s.pos.norm();
            assert!(energy.abs() < 1e-12, "energy {energy} at t = {t} yr");
        }
        // Inclination of the relative orbit's angular momentum.
        let s = f.relative_state(1000.0);
        let h = s.pos.cross(&s.vel).normalize();
        assert!((h.z.acos() / DEG2RAD - 70.0).abs() < 1e-9);
    }

    #[test]
    fn time_at_distance_inverts_the_trajectory() {
        let f = FlybyParams::model_a();
        let t = f.time_at_distance(2000.0);
        assert!((f.relative_state(t).pos.norm() - 2000.0).abs() < 1e-6);
        // ~5 kyr each way for this encounter.
        assert!((t / YEAR_DAYS) > 4000.0 && (t / YEAR_DAYS) < 7000.0);
    }

    #[test]
    fn unperturbed_circular_orbit_is_preserved_by_the_integrator() {
        // A massless perturber leaves the tracer on its circular orbit.
        let mut f = FlybyParams::model_a();
        f.mass_solar = 0.0;
        let r = 40.0;
        let v = (GM_SUN / r).sqrt();
        let t1 = 20.0 * TWO_PI * (r.powi(3) / GM_SUN).sqrt();
        let (p, vv) = integrate_tracer(
            &f,
            Vector3::new(r, 0.0, 0.0),
            Vector3::new(0.0, v, 0.0),
            0.0,
            t1,
            1e-10,
        );
        assert!((p.norm() - r).abs() < 1e-6, "r = {}", p.norm());
        assert!((vv.norm() - v).abs() < 1e-10);
    }
}
