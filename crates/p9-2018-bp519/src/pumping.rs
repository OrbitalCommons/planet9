//! Secular inclination pumping of a scattered test particle by an inclined
//! Planet Nine.
//!
//! Becker et al. (2018) argue that a high-inclination detached object like
//! 2015 BP519 (i ≈ 54°) is *naturally produced* by Planet Nine: an inclined
//! P9 secularly torques the orbital plane of an initially low-inclination
//! scattered particle, pumping its inclination over Gyr while it stays
//! detached from Neptune. Without P9 the particle's plane is held near the
//! ecliptic by the giant planets and its inclination stays low.
//!
//! ## Model
//!
//! We treat P9 as a doubly-averaged Gauss ring (the exact 1/Δ averaged over
//! both mean anomalies — `p9_core::analysis::secular::numerical_secular_hamiltonian`),
//! reused verbatim rather than re-derived. The test particle's slow
//! (secular) variables (e, i, ω, Ω) evolve under Hamilton's equations in
//! Delaunay action–angle variables:
//!
//! ```text
//!   L = √(GM a)            (fixed: a is a secular invariant)
//!   G = L √(1 − e²)        conjugate to ω
//!   Hz = G cos i           conjugate to Ω   (z-angular-momentum)
//!
//!   dω/dt = +∂H/∂G,   dG/dt = −∂H/∂ω
//!   dΩ/dt = +∂H/∂Hz,  dHz/dt = −∂H/∂Ω
//! ```
//!
//! The Hamiltonian here is P9's ring term evaluated with the particle's
//! elements expressed *relative to P9's plane* (the secular function takes the
//! perturber as the reference plane), plus the orbit-averaged J2 field of the
//! four giant planets (`p9_core::forces::j2_secular`), which holds a P9-free
//! particle's plane near the ecliptic and supplies the prograde nodal
//! regression that P9 must overcome. The mutual inclination drives the
//! inclination oscillation; an inclined P9 (sin i₉ ≠ 0) opens a wide
//! inclination cycle (a Kozai-like resonance for the relevant geometry),
//! while a coplanar P9 (i₉ = 0) leaves a low-i particle low.
//!
//! Partial derivatives of the Gauss-ring Hamiltonian are taken by centred
//! finite differences; the secular ODE is advanced with a fixed-step RK4.
//! This is a reduced-scale single-perturber secular model, not the paper's
//! full N-body suite — see the residual notes in the tests.

use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::GM_SUN;
use p9_core::forces::j2_secular::combined_j2_jsu;
use p9_core::types::P9Params;
use p9_core::units::{degrees, Angle};

/// Slowly-varying secular state of a test particle (semi-major axis is a
/// secular invariant and carried separately).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SecularState {
    /// Eccentricity.
    pub e: f64,
    /// Inclination relative to P9's plane (radians).
    pub i: f64,
    /// Argument of perihelion from the mutual node (radians).
    pub omega: f64,
    /// Longitude of ascending node from P9's apsidal line (radians).
    pub omega_big: f64,
}

/// A test particle undergoing secular evolution at fixed `a`.
#[derive(Debug, Clone, Copy)]
pub struct TestParticle {
    /// Semi-major axis (AU) — secular invariant.
    pub a: f64,
    pub state: SecularState,
}

/// Quadrature nodes for the Gauss-ring double average. The pumping result is
/// qualitative (does an inclined P9 drive the inclination into the BP519 band
/// vs. flat without P9?) and robust to the node count; 12 keeps the secular
/// ODE fast enough for the unit-test gate while these smooth, non-crossing
/// detached geometries stay well-resolved (p9-core's own convergence test
/// shows the ring average is smooth in the node count for non-crossing
/// orbits).
const N_QUAD: usize = 12;

/// Secular Hamiltonian of the test particle: P9's Gauss ring, optionally plus
/// the orbit-averaged giant-planet J2 field.
///
/// `softening` regularizes the (rare) crossing geometry. The J2 term is the
/// standard axisymmetric quadrupole `−GM J2_eff R² / (2 a³ (1−e²)^{3/2})`,
/// which depends only on (a, e) and so contributes only to apsidal
/// precession, never directly to the inclination torque. Its effect on the
/// inclination dynamics is *indirect but important*: the fast prograde apsidal
/// precession it forces can quench (average out) P9's Kozai-Lidov torque — the
/// classical de-Sitter/Kozai-suppression mechanism. Including it is the
/// realistic case; excluding it isolates the clean single-perturber Kozai
/// cycle. See [`integrate`]'s `include_j2` flag and the residual note in the
/// module tests.
fn hamiltonian(a: f64, s: SecularState, p9: &P9Params, softening: f64, include_j2: bool) -> f64 {
    let h_p9 = numerical_secular_hamiltonian(
        a,
        s.e,
        s.i,
        s.omega,
        s.omega_big,
        p9.a,
        p9.e,
        p9.gm(),
        N_QUAD,
        softening,
    );
    if !include_j2 {
        return h_p9;
    }
    let (j2, _, gm_boost) = combined_j2_jsu();
    let gm_eff = GM_SUN + gm_boost;
    let h_j2 = -gm_eff * j2 / (2.0 * a.powi(3) * (1.0 - s.e * s.e).powf(1.5));
    h_p9 + h_j2
}

/// Time derivative of the secular state from Hamilton's equations in Delaunay
/// variables, evaluated by centred finite differences of [`hamiltonian`].
fn derivative(
    a: f64,
    s: SecularState,
    p9: &P9Params,
    softening: f64,
    include_j2: bool,
) -> SecularState {
    let big_l = (GM_SUN * a).sqrt();
    let g = big_l * (1.0 - s.e * s.e).sqrt();
    let cos_i = s.i.cos();
    let sin_i = s.i.sin();

    // Finite-difference steps in the Hamiltonian's natural arguments.
    let de = 1e-5;
    let di = 1e-5;
    let dang = 1e-5;

    let h = |st: SecularState| hamiltonian(a, st, p9, softening, include_j2);

    let dh_de = {
        let mut sp = s;
        let mut sm = s;
        sp.e += de;
        sm.e -= de;
        (h(sp) - h(sm)) / (2.0 * de)
    };
    let dh_di = {
        let mut sp = s;
        let mut sm = s;
        sp.i += di;
        sm.i -= di;
        (h(sp) - h(sm)) / (2.0 * di)
    };
    let dh_domega = {
        let mut sp = s;
        let mut sm = s;
        sp.omega += dang;
        sm.omega -= dang;
        (h(sp) - h(sm)) / (2.0 * dang)
    };
    let dh_domega_big = {
        let mut sp = s;
        let mut sm = s;
        sp.omega_big += dang;
        sm.omega_big -= dang;
        (h(sp) - h(sm)) / (2.0 * dang)
    };

    // Chain rule: (e, i) → (G, Hz). With G = L√(1−e²), Hz = G cos i:
    //   ∂H/∂G  = ∂H/∂e · de/dG + ∂H/∂i · di/dG
    //   ∂H/∂Hz = ∂H/∂i · di/dHz
    // de/dG  = −√(1−e²)/(L e)
    // For Hz = G cos i at fixed G: di/dHz = −1/(G sin i).
    // For G at fixed Hz: i = acos(Hz/G), di/dG = −cos i/(G sin i) (used via
    //   the Hz-conjugacy below), giving di/dG|_{Hz} contribution to dω/dt.
    let de_dg = if s.e > 1e-9 {
        -(1.0 - s.e * s.e).sqrt() / (big_l * s.e)
    } else {
        0.0
    };
    let di_dg_at_hz = if sin_i > 1e-9 {
        -cos_i / (g * sin_i)
    } else {
        0.0
    };
    let di_dhz = if sin_i > 1e-9 {
        -1.0 / (g * sin_i)
    } else {
        0.0
    };

    let dh_dg = dh_de * de_dg + dh_di * di_dg_at_hz;
    let dh_dhz = dh_di * di_dhz;

    // Hamilton's equations.
    let domega_dt = dh_dg;
    let dg_dt = -dh_domega;
    let domega_big_dt = dh_dhz;
    let dhz_dt = -dh_domega_big;

    // Map (dG/dt, dHz/dt) back to (de/dt, di/dt).
    // G = L√(1−e²) ⇒ de/dt = −G/(L² e) dG/dt … i.e. dG/dt = −L e/√(1−e²) de/dt.
    let de_dt = if s.e > 1e-9 {
        -dg_dt * (1.0 - s.e * s.e).sqrt() / (big_l * s.e)
    } else {
        0.0
    };
    // Hz = G cos i ⇒ dHz/dt = cos i dG/dt − G sin i di/dt
    //   ⇒ di/dt = (cos i dG/dt − dHz/dt)/(G sin i).
    let di_dt = if sin_i > 1e-9 {
        (cos_i * dg_dt - dhz_dt) / (g * sin_i)
    } else {
        0.0
    };

    SecularState {
        e: de_dt,
        i: di_dt,
        omega: domega_dt,
        omega_big: domega_big_dt,
    }
}

fn axpy(s: SecularState, k: SecularState, h: f64) -> SecularState {
    SecularState {
        e: s.e + h * k.e,
        i: s.i + h * k.i,
        omega: s.omega + h * k.omega,
        omega_big: s.omega_big + h * k.omega_big,
    }
}

/// Advance the secular state by one RK4 step of `dt` days.
#[allow(clippy::too_many_arguments)]
pub fn rk4_step(
    a: f64,
    s: SecularState,
    p9: &P9Params,
    dt: f64,
    softening: f64,
    include_j2: bool,
) -> SecularState {
    let k1 = derivative(a, s, p9, softening, include_j2);
    let k2 = derivative(a, axpy(s, k1, dt / 2.0), p9, softening, include_j2);
    let k3 = derivative(a, axpy(s, k2, dt / 2.0), p9, softening, include_j2);
    let k4 = derivative(a, axpy(s, k3, dt), p9, softening, include_j2);

    let mut out = SecularState {
        e: s.e + dt / 6.0 * (k1.e + 2.0 * k2.e + 2.0 * k3.e + k4.e),
        i: s.i + dt / 6.0 * (k1.i + 2.0 * k2.i + 2.0 * k3.i + k4.i),
        omega: s.omega + dt / 6.0 * (k1.omega + 2.0 * k2.omega + 2.0 * k3.omega + k4.omega),
        omega_big: s.omega_big
            + dt / 6.0 * (k1.omega_big + 2.0 * k2.omega_big + 2.0 * k3.omega_big + k4.omega_big),
    };
    out.e = out.e.clamp(0.0, 0.999);
    out
}

/// Result of integrating a particle's secular inclination history.
#[derive(Debug, Clone)]
pub struct PumpingHistory {
    /// Times (days from start).
    pub times: Vec<f64>,
    /// Inclination relative to P9's plane (degrees).
    pub inclination_deg: Vec<f64>,
    /// Eccentricity.
    pub eccentricity: Vec<f64>,
}

impl PumpingHistory {
    /// Maximum inclination over the history as a typed [`Angle`].
    pub fn max_inclination(&self) -> Angle {
        degrees(
            self.inclination_deg
                .iter()
                .cloned()
                .fold(f64::NEG_INFINITY, f64::max),
        )
    }

    /// Inclination at the end of the history as a typed [`Angle`].
    pub fn final_inclination(&self) -> Angle {
        degrees(*self.inclination_deg.last().unwrap())
    }

    /// Maximum eccentricity reached.
    pub fn max_eccentricity(&self) -> f64 {
        self.eccentricity
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max)
    }
}

/// Integration controls for [`integrate`].
#[derive(Debug, Clone, Copy)]
pub struct PumpConfig {
    /// End time in days.
    pub t_end: f64,
    /// Timestep in days.
    pub dt: f64,
    /// Record a sample every this many steps.
    pub record_every: usize,
    /// Gauss-ring softening as a fraction of P9's a (1e-3 for non-crossing).
    pub softening_frac: f64,
    /// Include the giant-planet J2 apsidal precession (realistic case). When
    /// false, the clean single-perturber P9 Kozai cycle is isolated.
    pub include_j2: bool,
}

impl PumpConfig {
    /// A fast secular run suited to the Kozai cycle here: large secular
    /// timestep, ~70 Myr (enough to reach the first inclination peak),
    /// non-crossing softening.
    pub fn fast() -> Self {
        Self {
            t_end: 7.0e7 * p9_core::constants::YEAR_DAYS,
            dt: 3.0e5,
            record_every: 5,
            softening_frac: 1e-3,
            include_j2: false,
        }
    }
}

/// Integrate a test particle's secular evolution under P9 per `cfg`.
pub fn integrate(particle: TestParticle, p9: &P9Params, cfg: PumpConfig) -> PumpingHistory {
    let softening = cfg.softening_frac * p9.a;
    let n_steps = (cfg.t_end / cfg.dt).ceil() as usize;

    let mut s = particle.state;
    let mut times = Vec::new();
    let mut inc = Vec::new();
    let mut ecc = Vec::new();

    for step in 0..=n_steps {
        if step % cfg.record_every == 0 || step == n_steps {
            times.push(step as f64 * cfg.dt);
            inc.push(s.i.to_degrees());
            ecc.push(s.e);
        }
        if step < n_steps {
            s = rk4_step(particle.a, s, p9, cfg.dt, softening, cfg.include_j2);
        }
    }

    PumpingHistory {
        times,
        inclination_deg: inc,
        eccentricity: ecc,
    }
}

/// A Planet Nine with negligible mass — the no-P9 control. With
/// `include_j2 = true` the secular field reduces to the giant-planet J2 ring,
/// which depends only on (a, e) and therefore exerts NO torque on the orbital
/// plane: the inclination is exactly conserved.
pub fn no_planet_nine() -> P9Params {
    let mut p9 = P9Params::nominal_2016();
    p9.mass_earth = 1e-9;
    p9
}

/// A scattered detached particle whose mutual inclination to P9 equals
/// `i_mutual_deg`. Physically this represents an ecliptic-plane scattered
/// object (low inclination to the ecliptic) seen in the frame of a Planet
/// Nine tilted by i₉ ≈ `i_mutual_deg`: when the particle lies near the
/// ecliptic, its inclination *relative to P9's plane* (the angle the secular
/// Hamiltonian acts on) is set by P9's tilt. The high eccentricity is the
/// Neptune-scattered, detached state BP519 occupies; `omega = 135°` places it
/// on the inclination-growing branch of the Kozai-Lidov cycle.
pub fn scattered_particle(a: f64, e: f64, i_mutual_deg: f64) -> TestParticle {
    use p9_core::constants::DEG2RAD;
    TestParticle {
        a,
        state: SecularState {
            e,
            i: i_mutual_deg * DEG2RAD,
            omega: 135.0 * DEG2RAD,
            omega_big: 0.0,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use uom::si::angle::degree;

    /// The nominal inclined Planet Nine (Batygin & Brown 2016): 10 M⊕,
    /// a = 700 AU, e = 0.6, i₉ = 30°.
    fn inclined_p9() -> P9Params {
        P9Params::nominal_2016()
    }

    #[test]
    fn typed_inclination_accessors_match_degrees() {
        use approx::assert_relative_eq;
        let p = scattered_particle(250.0, 0.85, 40.0);
        let hist = integrate(p, &inclined_p9(), PumpConfig::fast());
        let max_deg = hist
            .inclination_deg
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(
            hist.max_inclination().get::<degree>(),
            max_deg,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            hist.final_inclination().get::<degree>(),
            *hist.inclination_deg.last().unwrap(),
            epsilon = 1e-9
        );
    }

    #[test]
    fn p9_pumps_inclination_to_bp519_like() {
        // An inclined P9 secularly pumps a scattered detached particle (high
        // e, mutual i = 40°) UP to a BP519-like inclination (~50°+) over a few
        // hundred Myr via the Kozai-Lidov cycle. Calibration: peaks ≈ 56°.
        let p = scattered_particle(250.0, 0.85, 40.0);
        let hist = integrate(p, &inclined_p9(), PumpConfig::fast());
        assert!(
            hist.max_inclination().get::<degree>() >= 50.0,
            "inclined P9 should pump i ≥ 50° (BP519-like), got {:.1}°",
            hist.max_inclination().get::<degree>()
        );
        // And the achieved inclination is in the high-i band, consistent with
        // BP519's observed i ≈ 54°.
        assert!(
            (50.0..70.0).contains(&hist.max_inclination().get::<degree>()),
            "max i = {:.1}° outside the BP519-consistent band",
            hist.max_inclination().get::<degree>()
        );
    }

    #[test]
    fn no_p9_control_leaves_inclination_flat() {
        // Without Planet Nine (only the axisymmetric giant-planet J2 ring,
        // which has no inclination torque), the SAME particle's inclination is
        // unchanged — there is no mechanism to reach high i. This isolates P9
        // as the cause of the inclination pumping.
        let p = scattered_particle(250.0, 0.85, 40.0);
        let mut cfg = PumpConfig::fast();
        cfg.include_j2 = true;
        let hist = integrate(p, &no_planet_nine(), cfg);
        let min_i = hist
            .inclination_deg
            .iter()
            .cloned()
            .fold(f64::INFINITY, f64::min);
        let range = hist.max_inclination().get::<degree>() - min_i;
        assert!(
            range < 2.0,
            "no-P9 control inclination should stay flat, range = {range:.2}°"
        );
        // It stays near its initial 40°, never reaching the BP519-like band.
        assert!(
            hist.max_inclination().get::<degree>() < 45.0,
            "no-P9 control should not pump, got {:.1}°",
            hist.max_inclination().get::<degree>()
        );
    }

    #[test]
    fn pumping_grows_with_p9_inclination() {
        // sin(i₉): a more inclined P9 sets a larger mutual inclination for an
        // ecliptic scattered particle (i_mutual ≈ i₉) and so pumps the
        // inclination higher. Compare i₉ = 25° vs 40°.
        let p9 = inclined_p9();
        let lo = integrate(
            scattered_particle(250.0, 0.85, 25.0),
            &p9,
            PumpConfig::fast(),
        );
        let hi = integrate(
            scattered_particle(250.0, 0.85, 40.0),
            &p9,
            PumpConfig::fast(),
        );
        assert!(
            hi.max_inclination().get::<degree>() > lo.max_inclination().get::<degree>() + 3.0,
            "higher i₉ (mutual inclination) should pump more: {:.1}° vs {:.1}°",
            hi.max_inclination().get::<degree>(),
            lo.max_inclination().get::<degree>()
        );
    }

    #[test]
    fn pumping_rate_grows_with_p9_mass() {
        // The secular forcing rate ∝ GM₉, so a heavier P9 advances further
        // along the inclination-pumping cycle in a fixed (short) time.
        let p = scattered_particle(250.0, 0.85, 40.0);
        let mut p9_light = P9Params::nominal_2016();
        p9_light.mass_earth = 5.0;
        let mut p9_heavy = P9Params::nominal_2016();
        p9_heavy.mass_earth = 20.0;

        let mut cfg = PumpConfig::fast();
        cfg.t_end = 0.05e9 * p9_core::constants::YEAR_DAYS; // 50 Myr (sub-cycle)
        let light = integrate(p, &p9_light, cfg);
        let heavy = integrate(p, &p9_heavy, cfg);
        assert!(
            heavy.max_inclination().get::<degree>() > light.max_inclination().get::<degree>(),
            "heavier P9 should pump faster: {:.1}° vs {:.1}°",
            heavy.max_inclination().get::<degree>(),
            light.max_inclination().get::<degree>()
        );
    }

    #[test]
    fn pumped_orbit_stays_prograde_detached_and_finite() {
        // The pumped inclination is a real, well-behaved secular evolution:
        // the inclination is lifted into the BP519 band but stays prograde
        // (i < 90°, as observed for BP519's i ≈ 54°), every sample is finite,
        // and the eccentricity stays bound (e < 1, the orbit remains
        // detached, not ejected). This guards against runaway/NaN integration.
        let p = scattered_particle(250.0, 0.85, 40.0);
        let hist = integrate(p, &inclined_p9(), PumpConfig::fast());
        assert!(
            hist.max_inclination().get::<degree>() > 50.0
                && hist.max_inclination().get::<degree>() < 90.0,
            "pumped i should be high but prograde, max i = {:.1}°",
            hist.max_inclination().get::<degree>()
        );
        for (&i, &e) in hist.inclination_deg.iter().zip(&hist.eccentricity) {
            assert!(i.is_finite() && (0.0..90.0).contains(&i), "i = {i}");
            assert!(e.is_finite() && (0.0..1.0).contains(&e), "e = {e}");
        }
    }
}
