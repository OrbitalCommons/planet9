use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use crate::constants::*;

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

    pub fn distance(&self) -> f64 {
        self.pos.norm()
    }

    pub fn speed(&self) -> f64 {
        self.vel.norm()
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
    /// Perihelion distance q = a(1-e) in AU
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Aphelion distance Q = a(1+e) in AU
    pub fn aphelion(&self) -> f64 {
        self.a * (1.0 + self.e)
    }

    /// Semi-latus rectum p = a(1-e^2)
    pub fn semi_latus_rectum(&self) -> f64 {
        self.a * (1.0 - self.e * self.e)
    }

    /// Orbital period in days (only meaningful for bound orbits)
    pub fn period(&self, gm: f64) -> f64 {
        TWO_PI * (self.a.powi(3) / gm).sqrt()
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

    /// Perihelion distance in AU
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
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
        // TODO: Estimate radius from mass-radius relation (r = (m/3M_earth) * R_earth)
        let radius_au = 3.0 * EARTH_MASS_SOLAR * self.mass_earth * 6_371.0 / AU_KM;
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
        // Hyperbolic case
        let sin_nu = -(elem.e * elem.e - 1.0).sqrt() * ea.sinh() / (elem.e * ea.cosh() - 1.0);
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
            // No inclination but has eccentricity
            let mut w = (e_vec.x / e).clamp(-1.0, 1.0).acos();
            if e_vec.y < 0.0 {
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
            let mut nu = (r.x / r_mag).clamp(-1.0, 1.0).acos();
            if v.x > 0.0 {
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

/// Solve Kepler's equation M = E - e*sin(E) for eccentric anomaly.
/// Uses Newton-Raphson iteration.
fn kepler_equation(e: f64, mean_anomaly: f64) -> f64 {
    if e < 1.0 {
        // Elliptic case
        let m = mean_anomaly % TWO_PI;
        let mut ea = if e > 0.8 { std::f64::consts::PI } else { m };

        for _ in 0..50 {
            let delta = (ea - e * ea.sin() - m) / (1.0 - e * ea.cos());
            ea -= delta;
            if delta.abs() < 1e-15 {
                break;
            }
        }
        ea
    } else {
        // Hyperbolic case
        let mut ha = mean_anomaly.signum() * mean_anomaly.abs().ln().max(1.0);
        for _ in 0..50 {
            let delta = (e * ha.sinh() - ha - mean_anomaly) / (e * ha.cosh() - 1.0);
            ha -= delta;
            if delta.abs() < 1e-15 {
                break;
            }
        }
        ha
    }
}
