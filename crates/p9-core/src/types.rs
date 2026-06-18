use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::constants::*;
use crate::units::{self, GravitationalParameter, Length, Time, Velocity};

/// 3D position + velocity state vector in heliocentric ecliptic J2000.
/// Units: AU for position, AU/day for velocity.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct StateVector {
    pub pos: Vector3<f64>,
    pub vel: Vector3<f64>,
}

impl StateVector {
    pub fn new(pos: Vector3<f64>, vel: Vector3<f64>) -> Self {
        Self { pos, vel }
    }

    pub fn zero() -> Self {
        Self {
            pos: Vector3::zeros(),
            vel: Vector3::zeros(),
        }
    }

    pub fn speed(&self) -> f64 {
        self.vel.norm()
    }

    /// Heliocentric distance as a dimension-checked [`Length`]. (The vector
    /// storage stays `f64`/AU because `nalgebra`'s norm cannot carry units;
    /// this is the typed boundary accessor.)
    pub fn distance_typed(&self) -> Length {
        units::au(self.pos.norm())
    }

    /// Speed as a dimension-checked [`Velocity`] (AU/day storage).
    pub fn speed_typed(&self) -> Velocity {
        units::au(self.speed()) / units::days(1.0)
    }
}

/// Classical Keplerian orbital elements.
/// Angles in radians, distances in AU.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct OrbitalElements {
    /// Semi-major axis (AU). Negative for hyperbolic orbits.
    pub a: f64,
    /// Eccentricity (dimensionless)
    pub e: f64,
    /// Inclination (radians, [0, pi])
    pub i: f64,
    /// Longitude of ascending node (radians, [0, 2pi))
    pub omega_big: f64,
    /// Argument of perihelion (radians, [0, 2pi))
    pub omega: f64,
    /// Mean anomaly (radians, [0, 2pi))
    pub mean_anomaly: f64,
}

impl OrbitalElements {
    /// Semi-latus rectum p = a(1-e^2)
    pub fn semi_latus_rectum(&self) -> f64 {
        self.a * (1.0 - self.e * self.e)
    }

    /// Longitude of perihelion (varpi = Omega + omega)
    pub fn longitude_of_perihelion(&self) -> f64 {
        (self.omega_big + self.omega) % TWO_PI
    }

    /// Mean motion (radians/day)
    pub fn mean_motion(&self, gm: f64) -> f64 {
        (gm / self.a.powi(3)).sqrt()
    }

    /// Convert to Cartesian state vector.
    /// `gm` is the central body GM in AU^3/day^2.
    pub fn to_state_vector(&self, gm: f64) -> StateVector {
        elements_to_cartesian(self, gm)
    }

    // ---- typed boundary accessors (dimension-checked views of the f64 fields) ----

    /// Semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        units::au(self.a)
    }
    /// Inclination as an [`Angle`].
    pub fn inclination(&self) -> units::Angle {
        units::radians(self.i)
    }
    /// Perihelion distance `q = a(1-e)` as a [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        units::au(self.a * (1.0 - self.e))
    }
    /// Aphelion distance `Q = a(1+e)` as a [`Length`].
    pub fn aphelion_typed(&self) -> Length {
        units::au(self.a * (1.0 + self.e))
    }
    /// Orbital period as a [`Time`]. `gm` is the central-body gravitational
    /// parameter in AU³/day² (the workspace's native unit).
    pub fn period_typed(&self, gm: f64) -> Time {
        units::days(TWO_PI * (self.a.powi(3) / gm).sqrt())
    }
}

/// A massive body (star, planet, or Planet Nine).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MassiveBody {
    pub name: String,
    /// GM in AU^3/day^2
    pub gm: f64,
    /// Mass in solar masses
    pub mass: f64,
    pub state: StateVector,
    /// Equatorial radius in AU (for Hill sphere computation)
    pub radius_au: f64,
    /// J2 oblateness (optional)
    pub j2: Option<f64>,
    /// J4 oblateness (optional)
    pub j4: Option<f64>,
}

impl MassiveBody {
    /// Hill sphere radius: r_H = a * (m / 3*M_central)^(1/3)
    pub fn hill_radius(&self, distance_from_central: f64) -> f64 {
        distance_from_central * (self.mass / 3.0).cbrt()
    }
}

/// A massless test particle.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Particle {
    pub id: u64,
    pub state: StateVector,
    pub active: bool,
}

/// Planet Nine orbital parameter set.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct P9Params {
    /// Mass in Earth masses
    pub mass_earth: f64,
    /// Semi-major axis in AU
    pub a: f64,
    /// Eccentricity
    pub e: f64,
    /// Inclination in radians
    pub i: f64,
    /// Argument of perihelion in radians
    pub omega: f64,
    /// Longitude of ascending node in radians
    pub omega_big: f64,
    /// Mean anomaly in radians
    pub mean_anomaly: f64,
}

impl P9Params {
    /// Batygin & Brown (2016) nominal parameters
    pub fn nominal_2016() -> Self {
        Self {
            mass_earth: 10.0,
            a: 700.0,
            e: 0.6,
            i: 30.0 * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    /// Batygin & Brown (2016) variant used in inclined-TNOs paper
    pub fn inclined_tnos_2016() -> Self {
        Self {
            mass_earth: 10.0,
            a: 600.0,
            e: 0.5,
            i: 30.0 * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    /// Batygin et al. (2019) revised parameters (review paper best fit)
    pub fn revised_2019() -> Self {
        Self {
            mass_earth: 5.0,
            a: 500.0,
            e: 0.25,
            i: 20.0 * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    /// Brown & Batygin (2021) MCMC posterior median
    pub fn mcmc_2021() -> Self {
        Self {
            mass_earth: 6.2,
            a: 380.0,
            e: 0.3,
            i: 16.0 * DEG2RAD,
            omega: 150.0 * DEG2RAD,
            omega_big: 100.0 * DEG2RAD,
            mean_anomaly: 0.0,
        }
    }

    /// Mass in solar masses
    pub fn mass_solar(&self) -> f64 {
        self.mass_earth * EARTH_MASS_SOLAR
    }

    /// GM in AU^3/day^2
    pub fn gm(&self) -> f64 {
        self.mass_solar() * GM_SUN
    }

    // ---- typed boundary accessors (dimension-checked views of the f64 fields) ----

    /// Planet Nine mass as a [`Mass`] (from the Earth-mass field).
    pub fn mass(&self) -> units::Mass {
        units::earth_masses(self.mass_earth)
    }
    /// Semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        units::au(self.a)
    }
    /// Inclination as an [`Angle`].
    pub fn inclination(&self) -> units::Angle {
        units::radians(self.i)
    }
    /// Perihelion distance as a [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        units::au(self.a * (1.0 - self.e))
    }
    /// Gravitational parameter `GM` as a typed [`GravitationalParameter`].
    pub fn gm_typed(&self) -> GravitationalParameter {
        units::gm_from_au3_day2(self.gm())
    }

    /// Convert to a MassiveBody at a given mean anomaly
    pub fn to_body(&self) -> MassiveBody {
        let elements = OrbitalElements {
            a: self.a,
            e: self.e,
            i: self.i,
            omega_big: self.omega_big,
            omega: self.omega,
            mean_anomaly: self.mean_anomaly,
        };
        let state = elements.to_state_vector(GM_SUN);
        // Neptune-anchored mass-radius relation, R ≈ 3.4 R⊕ at 10 M⊕
        // (see analysis::photometry::mass_radius_neptunian for the source).
        let radius_km =
            crate::analysis::photometry::mass_radius_neptunian(self.mass_earth) * EARTH_RADIUS_KM;
        let radius_au = radius_km / AU_KM;
        MassiveBody {
            name: "Planet Nine".to_string(),
            gm: self.gm(),
            mass: self.mass_solar(),
            state,
            radius_au,
            j2: None,
            j4: None,
        }
    }
}

/// Configuration for a simulation run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    /// Timestep in days
    pub dt: f64,
    /// Start time in days from epoch
    pub t_start: f64,
    /// End time in days from epoch
    pub t_end: f64,
    /// Remove particles closer than this to the Sun (AU)
    pub removal_inner_au: f64,
    /// Remove particles farther than this from the Sun (AU)
    pub removal_outer_au: f64,
    /// Record snapshot every this many days
    pub snapshot_interval_days: f64,
    /// Hill sphere multiplier for hybrid integrator switching
    pub hybrid_changeover_hill: f64,
    /// Bulirsch-Stoer accuracy parameter
    pub bs_epsilon: f64,
}

impl SimConfig {
    /// Standard config for 4 Gyr with 300 day timestep (Batygin & Brown 2016)
    pub fn standard_4gyr() -> Self {
        Self {
            dt: 300.0,
            t_start: 0.0,
            t_end: 4.0 * GYR_DAYS,
            removal_inner_au: 5.0,
            removal_outer_au: 10_000.0,
            snapshot_interval_days: 1_000_000.0 * YEAR_DAYS,
            hybrid_changeover_hill: 3.0,
            bs_epsilon: 1e-11,
        }
    }

    /// Pair this config with giant-planet starting bodies taken from the DE
    /// ephemeris at `epoch`, so an integration's t = `t_start` corresponds
    /// to a real survey/discovery epoch instead of J2000 only.
    ///
    /// Requires a cached DE kernel (never downloads); on a hermetic machine
    /// this returns `Err` and callers fall back to
    /// `initial_conditions::planets::giant_planets_j2000()`.
    pub fn with_bodies_at_epoch(
        self,
        epoch: &starfield::time::Time,
    ) -> Result<(Self, Vec<MassiveBody>), String> {
        let bodies = crate::initial_conditions::planets::giant_planets_at(epoch)?;
        Ok((self, bodies))
    }
}

// ============================================================
// Orbital element <-> Cartesian conversion
// ============================================================

/// Convert Keplerian elements to Cartesian state vector.
/// Internal implementation matching the standard Murray & Dermott formulation.
pub fn elements_to_cartesian(elem: &OrbitalElements, gm: f64) -> StateVector {
    let p = elem.semi_latus_rectum();

    // Solve Kepler's equation for eccentric anomaly
    let ea = kepler_equation(elem.e, elem.mean_anomaly);

    // True anomaly from eccentric anomaly
    let nu = if elem.e < 1.0 {
        let sin_nu = (1.0 - elem.e * elem.e).sqrt() * ea.sin() / (1.0 - elem.e * ea.cos());
        let cos_nu = (ea.cos() - elem.e) / (1.0 - elem.e * ea.cos());
        sin_nu.atan2(cos_nu)
    } else {
        // Hyperbolic case: sin ν = +√(e²−1)·sinh H / (e·cosh H − 1), so that
        // M > 0 (outbound, r·v > 0) maps to ν > 0. A negated sin ν mirrors the
        // state about perihelion (inbound/outbound swapped).
        let sin_nu = (elem.e * elem.e - 1.0).sqrt() * ea.sinh() / (elem.e * ea.cosh() - 1.0);
        let cos_nu = (elem.e - ea.cosh()) / (elem.e * ea.cosh() - 1.0);
        sin_nu.atan2(cos_nu)
    };

    // Distance
    let r = p / (1.0 + elem.e * nu.cos());

    // Position and velocity in orbital plane
    let r_orb = Vector3::new(r * nu.cos(), r * nu.sin(), 0.0);
    let mu_p = (gm / p).sqrt();
    let v_orb = Vector3::new(-mu_p * nu.sin(), mu_p * (elem.e + nu.cos()), 0.0);

    // Rotation from orbital plane to ecliptic
    let cos_o = elem.omega.cos();
    let sin_o = elem.omega.sin();
    let cos_i = elem.i.cos();
    let sin_i = elem.i.sin();
    let cos_big_o = elem.omega_big.cos();
    let sin_big_o = elem.omega_big.sin();

    let px = cos_big_o * cos_o - sin_big_o * sin_o * cos_i;
    let py = sin_big_o * cos_o + cos_big_o * sin_o * cos_i;
    let pz = sin_o * sin_i;

    let qx = -cos_big_o * sin_o - sin_big_o * cos_o * cos_i;
    let qy = -sin_big_o * sin_o + cos_big_o * cos_o * cos_i;
    let qz = cos_o * sin_i;

    let pos = Vector3::new(
        r_orb.x * px + r_orb.y * qx,
        r_orb.x * py + r_orb.y * qy,
        r_orb.x * pz + r_orb.y * qz,
    );

    let vel = Vector3::new(
        v_orb.x * px + v_orb.y * qx,
        v_orb.x * py + v_orb.y * qy,
        v_orb.x * pz + v_orb.y * qz,
    );

    StateVector { pos, vel }
}

// ============================================================
// True-anomaly geometry helper
// ============================================================

/// Heliocentric geometry of an orbit evaluated at one true anomaly.
#[derive(Debug, Clone, Copy)]
pub struct OrbitGeometry {
    /// True anomaly along the orbit (radians).
    pub true_anomaly: f64,
    /// Heliocentric position (AU, ecliptic J2000).
    pub position: Vector3<f64>,
    /// Heliocentric distance `r(ν) = a(1-e²)/(1 + e cos ν)` (AU).
    pub distance: f64,
}

/// Convert a true anomaly to the mean anomaly for an elliptic orbit
/// (`0 ≤ e < 1`). Inverse of the Kepler chain `M → E → ν`.
pub fn true_to_mean_anomaly(e: f64, nu: f64) -> f64 {
    // Eccentric anomaly E from true anomaly ν.
    let ea = ((1.0 - e) / (1.0 + e)).sqrt() * (nu / 2.0).tan();
    let ea = 2.0 * ea.atan();
    // Mean anomaly M = E - e sin E, wrapped to [0, 2π).
    let m = ea - e * ea.sin();
    m.rem_euclid(TWO_PI)
}

/// Heliocentric position of an orbit at true anomaly `nu` (radians) on the
/// orbit described by `params`. The `mean_anomaly` field of `params` is
/// overridden: the orbit's orientation (a, e, i, ω, Ω) is what matters, and the
/// caller sweeps ν. Built on [`elements_to_cartesian`].
pub fn position_at_true_anomaly(params: &P9Params, nu: f64) -> OrbitGeometry {
    let elements = OrbitalElements {
        a: params.a,
        e: params.e,
        i: params.i,
        omega_big: params.omega_big,
        omega: params.omega,
        mean_anomaly: true_to_mean_anomaly(params.e, nu),
    };
    let state: StateVector = elements_to_cartesian(&elements, GM_SUN);
    let position = state.pos;
    OrbitGeometry {
        true_anomaly: nu,
        position,
        distance: position.norm(),
    }
}

/// Analytic heliocentric distance `r(ν)` for the conic (AU), independent of the
/// Cartesian path; used as a cross-check.
pub fn helio_distance_at_true_anomaly(params: &P9Params, nu: f64) -> f64 {
    let p = params.a * (1.0 - params.e * params.e);
    p / (1.0 + params.e * nu.cos())
}

#[cfg(test)]
mod true_anomaly_geometry_tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;

    fn bb_orbit() -> P9Params {
        crate::data::ephemeris_constraint::brown_batygin_orbit()
    }

    #[test]
    fn true_mean_roundtrip_through_core_kepler() {
        let e = 0.6;
        for deg in [10.0, 73.0, 108.0, 200.0, 333.0] {
            let nu = deg * PI / 180.0;
            let m = true_to_mean_anomaly(e, nu);
            // Re-solve M → E → ν via the shared core solver and compare.
            let ea = solve_kepler(e, m);
            let sin_nu = (1.0 - e * e).sqrt() * ea.sin() / (1.0 - e * ea.cos());
            let cos_nu = (ea.cos() - e) / (1.0 - e * ea.cos());
            let nu_back = sin_nu.atan2(cos_nu).rem_euclid(TWO_PI);
            assert_relative_eq!(nu_back, nu.rem_euclid(TWO_PI), epsilon = 1e-9);
        }
    }

    #[test]
    fn perihelion_and_aphelion_distances() {
        let p = bb_orbit();
        let peri = position_at_true_anomaly(&p, 0.0);
        let apo = position_at_true_anomaly(&p, PI);
        assert_relative_eq!(peri.distance, p.a * (1.0 - p.e), epsilon = 1e-9);
        assert_relative_eq!(apo.distance, p.a * (1.0 + p.e), epsilon = 1e-9);
    }

    #[test]
    fn cartesian_distance_matches_conic_formula() {
        let p = bb_orbit();
        for deg in [0.0, 45.0, 108.0, 180.0, 270.0] {
            let nu = deg * PI / 180.0;
            let g = position_at_true_anomaly(&p, nu);
            let r_conic = helio_distance_at_true_anomaly(&p, nu);
            assert_relative_eq!(g.distance, r_conic, epsilon = 1e-7);
        }
    }
}

/// Convert Cartesian state vector to Keplerian elements.
pub fn cartesian_to_elements(state: &StateVector, gm: f64) -> OrbitalElements {
    let r = &state.pos;
    let v = &state.vel;
    let r_mag = r.norm();
    let v_mag = v.norm();

    // Specific angular momentum
    let h = r.cross(v);
    let h_mag = h.norm();

    // Node vector (z cross h)
    let n = Vector3::new(-h.y, h.x, 0.0);
    let n_mag = n.norm();

    // Eccentricity vector
    let e_vec = ((v_mag * v_mag - gm / r_mag) * r - r.dot(v) * v) / gm;
    let e = e_vec.norm();

    // Semi-major axis
    let energy = 0.5 * v_mag * v_mag - gm / r_mag;
    let a = if energy.abs() < 1e-30 {
        f64::INFINITY
    } else {
        -gm / (2.0 * energy)
    };

    // Inclination
    let i = (h.z / h_mag).clamp(-1.0, 1.0).acos();

    // Longitude of ascending node
    let omega_big = if n_mag < 1e-30 {
        0.0
    } else {
        let mut o = (n.x / n_mag).clamp(-1.0, 1.0).acos();
        if n.y < 0.0 {
            o = crate::constants::TWO_PI - o;
        }
        o
    };

    // Argument of perihelion
    let omega = if n_mag < 1e-30 || e < 1e-15 {
        if e >= 1e-15 {
            // Equatorial (no node): perihelion longitude measured along the
            // direction of motion; h.z's sign distinguishes prograde from
            // retrograde (a v.x or e_vec.y test alone breaks for i = π).
            let mut w = (e_vec.x / e).clamp(-1.0, 1.0).acos();
            if e_vec.y * h.z.signum() < 0.0 {
                w = TWO_PI - w;
            }
            w
        } else {
            0.0
        }
    } else {
        let cos_w = (n.dot(&e_vec) / (n_mag * e)).clamp(-1.0, 1.0);
        let mut w = cos_w.acos();
        if e_vec.z < 0.0 {
            w = TWO_PI - w;
        }
        w
    };

    // True anomaly
    let nu = if e < 1e-15 {
        if n_mag < 1e-30 {
            // Circular equatorial: true longitude along the direction of
            // motion, prograde/retrograde resolved by sign(h.z).
            let mut nu = (r.x / r_mag).clamp(-1.0, 1.0).acos();
            if r.y * h.z.signum() < 0.0 {
                nu = TWO_PI - nu;
            }
            nu
        } else {
            let cos_nu = (n.dot(r) / (n_mag * r_mag)).clamp(-1.0, 1.0);
            let mut nu = cos_nu.acos();
            if r.z < 0.0 {
                nu = TWO_PI - nu;
            }
            nu
        }
    } else {
        let cos_nu = (e_vec.dot(r) / (e * r_mag)).clamp(-1.0, 1.0);
        let mut nu = cos_nu.acos();
        if r.dot(v) < 0.0 {
            nu = TWO_PI - nu;
        }
        nu
    };

    // Mean anomaly from true anomaly
    let mean_anomaly = if e < 1.0 {
        let ea = ((1.0 - e) / (1.0 + e)).sqrt() * (nu / 2.0).tan();
        let ea = 2.0 * ea.atan();
        let ma = ea - e * ea.sin();
        if ma < 0.0 {
            ma + TWO_PI
        } else {
            ma
        }
    } else {
        // Hyperbolic
        let ha = ((e - 1.0) / (e + 1.0)).sqrt() * (nu / 2.0).tan();
        let ha = 2.0 * ha.atanh();
        e * ha.sinh() - ha
    };

    OrbitalElements {
        a,
        e,
        i,
        omega_big,
        omega,
        mean_anomaly,
    }
}

/// Solve Kepler's equation for the eccentric (e < 1) or hyperbolic (e > 1) anomaly.
///
/// Elliptic: M = E − e·sin E. Hyperbolic: M = e·sinh H − H.
///
/// Newton-Raphson with a bisection safeguard: the Kepler function is strictly
/// monotone in the anomaly, so the iterate is confined to a sign-changing
/// bracket and cannot diverge, even for e → 1. Starters: E₀ = π for e > 0.8
/// (elliptic, avoids the e·sin E ≈ M degeneracy near perihelion), Danby's
/// H₀ = ln(2M/e + 1.8) (hyperbolic).
///
/// This is the single Kepler solver for the workspace (several crates
/// previously carried unguarded Newton copies).
pub fn solve_kepler(e: f64, mean_anomaly: f64) -> f64 {
    if e < 1.0 {
        // Reduce M to (-pi, pi] and exploit the odd symmetry E(-M) = -E(M).
        let mut m = mean_anomaly % TWO_PI;
        if m > std::f64::consts::PI {
            m -= TWO_PI;
        } else if m < -std::f64::consts::PI {
            m += TWO_PI;
        }
        let sign = if m < 0.0 { -1.0 } else { 1.0 };
        let m = m.abs();

        // For m in [0, pi] the root lies in [m, m + e] (since 0 <= e sin E <= e).
        let mut lo = m;
        let mut hi = (m + e).min(std::f64::consts::PI + e);
        let mut ea = if e > 0.8 {
            std::f64::consts::PI
        } else {
            m + 0.85 * e
        };
        ea = ea.clamp(lo, hi);

        let mut converged = false;
        for _ in 0..100 {
            let f = ea - e * ea.sin() - m;
            if f > 0.0 {
                hi = ea;
            } else {
                lo = ea;
            }
            let fp = 1.0 - e * ea.cos();
            let delta = f / fp;
            let mut next = ea - delta;
            if !(lo..=hi).contains(&next) {
                // Newton left the bracket: bisect instead.
                next = 0.5 * (lo + hi);
            }
            let step = (next - ea).abs();
            ea = next;
            if step < 1e-15 * ea.abs().max(1.0) {
                converged = true;
                break;
            }
        }
        debug_assert!(
            converged || (ea - e * ea.sin() - m).abs() < 1e-12,
            "Kepler solver failed to converge: e = {e}, M = {mean_anomaly}"
        );
        sign * ea
    } else {
        // Hyperbolic, odd symmetry H(-M) = -H(M).
        let sign = if mean_anomaly < 0.0 { -1.0 } else { 1.0 };
        let m = mean_anomaly.abs();

        // Danby (1988) starter.
        let mut ha = (2.0 * m / e + 1.8).ln().max(1e-12);

        // Bracket: f(0) = -m <= 0; expand hi until f(hi) >= 0
        // (f' = e cosh H - 1 > 0, so the root is unique).
        let f = |h: f64| e * h.sinh() - h - m;
        let mut lo = 0.0_f64;
        let mut hi = ha.max(1.0);
        while f(hi) < 0.0 {
            hi *= 2.0;
        }
        ha = ha.clamp(lo, hi);

        let mut converged = false;
        for _ in 0..100 {
            let fv = f(ha);
            if fv > 0.0 {
                hi = ha;
            } else {
                lo = ha;
            }
            let fp = e * ha.cosh() - 1.0;
            let mut next = ha - fv / fp;
            if !(lo..=hi).contains(&next) {
                next = 0.5 * (lo + hi);
            }
            let step = (next - ha).abs();
            ha = next;
            if step < 1e-15 * ha.abs().max(1.0) {
                converged = true;
                break;
            }
        }
        debug_assert!(
            converged || f(ha).abs() < 1e-12 * m.max(1.0),
            "hyperbolic Kepler solver failed to converge: e = {e}, M = {mean_anomaly}"
        );
        sign * ha
    }
}

/// Backwards-compatible internal alias.
fn kepler_equation(e: f64, mean_anomaly: f64) -> f64 {
    solve_kepler(e, mean_anomaly)
}

#[cfg(test)]
mod typed_accessor_tests {
    use super::*;
    use approx::assert_relative_eq;
    use uom::si::length::astronomical_unit;
    use uom::si::mass::kilogram;

    #[test]
    fn p9params_typed_views_match_f64_fields() {
        let p = P9Params::nominal_2016();
        assert_relative_eq!(
            p.semi_major_axis().get::<astronomical_unit>(),
            p.a,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            p.perihelion_typed().get::<astronomical_unit>(),
            p.a * (1.0 - p.e),
            epsilon = 1e-12
        );
        // 10 M⊕ in kilograms.
        assert_relative_eq!(
            p.mass().get::<kilogram>(),
            10.0 * units::EARTH_MASS_KG,
            epsilon = 1.0
        );
        // Typed GM equals the f64 GM converted to SI.
        assert_relative_eq!(
            p.gm_typed().value,
            units::gm_from_au3_day2(p.gm()).value,
            epsilon = 1e6
        );
    }

    #[test]
    fn statevector_typed_views_match_norms() {
        let sv = OrbitalElements {
            a: 30.0,
            e: 0.1,
            i: 0.2,
            omega_big: 1.0,
            omega: 2.0,
            mean_anomaly: 0.5,
        }
        .to_state_vector(GM_SUN);
        assert_relative_eq!(
            sv.distance_typed().get::<astronomical_unit>(),
            sv.pos.norm(),
            epsilon = 1e-12
        );
    }
}
