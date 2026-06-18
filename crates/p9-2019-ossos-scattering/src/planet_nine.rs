//! Planet Nine's secular perihelion lifting.
//!
//! A distant, massive, eccentric P9 exerts a secular (orbit-averaged) torque
//! on a scattering TNO. Because the TNO's semi-major axis is (between Neptune
//! encounters) fixed, the secular dynamics conserve the doubly-averaged
//! Hamiltonian H(e, Δϖ); the angular momentum and hence the perihelion
//! q = a(1−e) oscillate. The relevant effect for OSSOS XV is that P9 *lifts*
//! perihelion at large a, pushing objects across the scattering/detached
//! divide and depleting the high-a scattering population.
//!
//! We do not invent the secular field: the Hamiltonian is the numerical
//! Gauss-ring average `p9_core::analysis::secular::numerical_secular_hamiltonian`
//! (the multipole expansion does not converge for the α = a/a₉ of this
//! problem), plus the giants' axisymmetric J2 precession term from
//! `p9_core::forces::j2_secular`.
//!
//! Two quantities are exposed:
//!  - `max_perihelion`: the cycle *amplitude* — the highest q reachable on the
//!    conserved level curve H(e, Δϖ) at fixed a (a Δϖ scan, restricted to the
//!    high-e branch where scattering TNOs live).
//!  - `forcing_rate` / `perihelion_lift`: the secular eccentricity-forcing
//!    rate |de/dt| ∝ GM₉ (from ∂H/∂ϖ), integrated over a sculpting time and
//!    capped at the amplitude. Because the single-perturber Gauss-ring
//!    Hamiltonian scales as GM₉, its level-curve amplitude is P9-mass-
//!    independent (a real feature of linear secular theory); the OSSOS-relevant
//!    mass dependence enters through the *rate*, hence the rate-limited lift.
//!
//! This is a coplanar, test-particle, secular model: it captures apsidal-driven
//! perihelion raising (the mechanism Kaib et al. invoke for the
//! scattering→detached transition) but not Neptune's simultaneous random walk
//! in a, nor the Kozai inclination coupling. Documented as such (see
//! `REPRODUCTION_NOTES.md`).

use p9_core::analysis::secular::numerical_secular_hamiltonian;
use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN, TWO_PI};
use p9_core::forces::j2_secular::combined_j2_jsu;
use p9_core::units::{au, earth_masses, gm_from_au3_day2, GravitationalParameter, Length, Mass};

/// Planet Nine orbital/mass parameters. Defaults are the Batygin & Brown
/// (2021) nominal P9: 6.2 M⊕, a₉ = 380 AU, e₉ = 0.2.
#[derive(Debug, Clone, Copy)]
pub struct PlanetNine {
    /// Mass in Earth masses.
    pub mass_earth: f64,
    /// Semi-major axis (AU).
    pub a_au: f64,
    /// Eccentricity.
    pub e: f64,
}

impl Default for PlanetNine {
    fn default() -> Self {
        // Batygin & Brown (2021) nominal P9.
        PlanetNine {
            mass_earth: 6.2,
            a_au: 380.0,
            e: 0.2,
        }
    }
}

impl PlanetNine {
    /// Heliocentric gravitational parameter GM₉ (AU³/day²).
    pub fn gm(&self) -> f64 {
        self.mass_earth * EARTH_MASS_SOLAR * GM_SUN
    }

    /// Mass as a typed [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }

    /// Semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }

    /// Gravitational parameter GM₉ as a typed [`GravitationalParameter`].
    pub fn gm_typed(&self) -> GravitationalParameter {
        gm_from_au3_day2(self.gm())
    }
}

/// Number of apsidal-angle samples used to trace a secular level curve.
const N_DVARPI: usize = 48;
/// Quadrature nodes per mean anomaly in the Gauss-ring average. 48 is enough
/// for the smooth (non-crossing) interior-particle geometry here and keeps the
/// population sweep affordable (N_DVARPI · n² evaluations per object).
const N_QUAD: usize = 48;

/// Coplanar secular (orbit-averaged J2) energy of a test particle in the
/// combined Jupiter+Saturn+Uranus quadrupole field (AU²/day²):
///
///   H_J2(a, e) = −(GM_sun / 2a) · J2_eff · (R_ref / a)² · (1 − e²)^{−3/2}
///
/// with `J2_eff` and `R_ref = 1 AU` the effective coefficients from
/// `p9_core::forces::j2_secular::combined_j2_jsu`. This axisymmetric term
/// drives the inner giants' apsidal precession; it depends on (a, e) only (no
/// Δϖ), so it does *not* by itself change e — but added to P9's
/// Δϖ-dependent torque it sets the competition that makes the P9 lift depend
/// on P9's mass and distance. Without it the pure single-perturber secular
/// Hamiltonian scales as GM₉ and its level curves (hence the e-excursion) are
/// independent of P9's mass, which is unphysical for OSSOS XV.
fn giants_secular_energy(a_au: f64) -> impl Fn(f64) -> f64 {
    let (j2_eff, _j4, _gm) = combined_j2_jsu();
    let r_ref = 1.0_f64;
    let coef = -(GM_SUN / (2.0 * a_au)) * j2_eff * (r_ref / a_au).powi(2);
    move |e: f64| coef * (1.0 - e * e).powf(-1.5)
}

/// E below which a scattering TNO would sit on the low-e branch of the
/// Gauss-ring Hamiltonian; the scattering population is always above this.
const E_BRANCH_MIN: f64 = 0.7;

/// Maximum perihelion (AU) the secular cycle can reach: the largest q on the
/// conserved level curve H(e, Δϖ) = H(e₀, Δϖ₀) at fixed a, the *amplitude* of
/// the perihelion oscillation regardless of how long the cycle takes.
///
/// The conserved Hamiltonian is P9's Δϖ-dependent torque
/// (`p9_core::analysis::secular::numerical_secular_hamiltonian`, softened by
/// 1% of a₉ as `phase_portrait` does) plus the giants' axisymmetric J2 term
/// (`giants_secular_energy`). The Gauss-ring H(e) is non-monotone (it peaks in
/// magnitude near e ≈ 0.5); scattering TNOs live on the high-e branch
/// (e ≳ 0.9), where H is monotone, so the level-set search is restricted to
/// e ≥ `E_BRANCH_MIN` — both the physics and a guard against spurious low-e
/// roots past the peak.
pub fn max_perihelion(p9: &PlanetNine, a_au: f64, q0_au: f64, dvarpi0: f64) -> f64 {
    let e0 = 1.0 - q0_au / a_au;
    let gm = p9.gm();
    let soft = 0.01 * p9.a_au;
    if e0 <= E_BRANCH_MIN {
        return q0_au;
    }
    let h_giants = giants_secular_energy(a_au);
    let h_total = |e: f64, dv: f64| {
        numerical_secular_hamiltonian(a_au, e, 0.0, dv, 0.0, p9.a_au, p9.e, gm, N_QUAD, soft)
            + h_giants(e)
    };
    let h0 = h_total(e0, dvarpi0);

    let mut e_min = e0;
    for j in 0..N_DVARPI {
        let dv = j as f64 / N_DVARPI as f64 * TWO_PI;
        let h_at = |e: f64| h_total(e, dv);
        let h_lo = h_at(E_BRANCH_MIN);
        let h_hi = h_at(e0);
        if (h0 - h_lo) * (h0 - h_hi) <= 0.0 {
            let (mut lo, mut hi) = (E_BRANCH_MIN, e0);
            for _ in 0..40 {
                let mid = 0.5 * (lo + hi);
                let hm = h_at(mid);
                if (hm - h0) * (h_hi - h0) <= 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            e_min = e_min.min(0.5 * (lo + hi));
        }
    }
    a_au * (1.0 - e_min)
}

/// Magnitude of P9's secular eccentricity-forcing rate |de/dt| (per day) for a
/// coplanar TNO at (a, q), evaluated at the apsidal phase of strongest torque.
///
/// From the orbit-averaged Lagrange planetary equation
///   de/dt = −(√(1−e²) / (n a² e)) ∂H/∂ϖ ,
/// with n the TNO mean motion and ∂H/∂ϖ taken from a centred finite difference
/// of P9's Gauss-ring Hamiltonian (the giants' term is axisymmetric and drops
/// out). This rate is ∝ GM₉, so a heavier / closer P9 drives the perihelion up
/// faster — the lever the OSSOS XV constraint pulls on.
pub fn forcing_rate(p9: &PlanetNine, a_au: f64, q0_au: f64) -> f64 {
    let e0 = 1.0 - q0_au / a_au;
    if e0 <= E_BRANCH_MIN || e0 <= 0.0 {
        return 0.0;
    }
    let gm = p9.gm();
    let soft = 0.01 * p9.a_au;
    let n = (GM_SUN / a_au.powi(3)).sqrt(); // mean motion (rad/day)
    let dphi = 1e-3;
    // Maximum |∂H/∂ϖ| over apsidal phase (the strongest-torque phase).
    let mut max_dh = 0.0_f64;
    for j in 0..N_DVARPI {
        let dv = j as f64 / N_DVARPI as f64 * TWO_PI;
        let hp = numerical_secular_hamiltonian(
            a_au,
            e0,
            0.0,
            dv + dphi,
            0.0,
            p9.a_au,
            p9.e,
            gm,
            N_QUAD,
            soft,
        );
        let hm = numerical_secular_hamiltonian(
            a_au,
            e0,
            0.0,
            dv - dphi,
            0.0,
            p9.a_au,
            p9.e,
            gm,
            N_QUAD,
            soft,
        );
        let dh = (hp - hm) / (2.0 * dphi);
        max_dh = max_dh.max(dh.abs());
    }
    (1.0 - e0 * e0).sqrt() / (n * a_au * a_au * e0) * max_dh
}

/// The perihelion lift Δq (AU) a TNO accumulates over `time_days` of secular
/// evolution under P9, capped at the secular-cycle amplitude.
///
/// The eccentricity drifts at `forcing_rate` (∝ GM₉) toward lower e (higher q);
/// over a time T the lift is a·|de/dt|·T, but it cannot exceed the conserved
/// cycle's reachable maximum (`max_perihelion`, averaged over starting apsidal
/// phase). The two together give the correct limits: a heavier P9 reaches the
/// cap sooner (so detaches more within a fixed cluster age), while no P9 gives
/// zero rate and zero lift.
pub fn perihelion_lift(p9: &PlanetNine, a_au: f64, q0_au: f64, time_days: f64) -> f64 {
    let rate = forcing_rate(p9, a_au, q0_au); // |de/dt|, per day
    let dq_rate = a_au * rate * time_days;

    // Cycle-amplitude cap, averaged over starting apsidal phase.
    const N_PHASE: usize = 8;
    let mut cap = 0.0;
    for k in 0..N_PHASE {
        let dv0 = k as f64 / N_PHASE as f64 * TWO_PI;
        cap += (max_perihelion(p9, a_au, q0_au, dv0) - q0_au).max(0.0);
    }
    cap /= N_PHASE as f64;

    dq_rate.min(cap).max(0.0)
}

/// Baseline secular-sculpting time used by the population detachment test
/// (50 Myr). This is the duration of a scattering episode over which P9's
/// secular torque acts coherently before a Neptune encounter resets the orbit;
/// it is short enough that the perihelion lift is in the rate-limited regime
/// (∝ GM₉) rather than saturated at the cycle-amplitude cap, so the detached
/// fraction responds monotonically to P9's mass. Documented as a model choice,
/// not a derived quantity (see `REPRODUCTION_NOTES.md`).
pub const SECULAR_SCULPTING_DAYS: f64 = 50.0e6 * 365.25;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn gm_scales_with_mass() {
        let small = PlanetNine {
            mass_earth: 5.0,
            ..Default::default()
        };
        let big = PlanetNine {
            mass_earth: 10.0,
            ..Default::default()
        };
        assert_relative_eq!(big.gm() / small.gm(), 2.0, epsilon = 1e-12);
        // Sanity: 6.2 M⊕ in solar masses times GM_SUN.
        let nominal = PlanetNine::default();
        assert_relative_eq!(
            nominal.gm(),
            6.2 * EARTH_MASS_SOLAR * GM_SUN,
            epsilon = 1e-18
        );
    }

    #[test]
    fn planet_nine_typed_accessors_match_f64() {
        use uom::si::length::astronomical_unit;
        use uom::si::mass::kilogram;
        let p9 = PlanetNine::default();
        assert_relative_eq!(
            p9.semi_major_axis().get::<astronomical_unit>(),
            p9.a_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            p9.mass().get::<kilogram>(),
            p9.mass_earth * p9_core::units::EARTH_MASS_KG,
            epsilon = 1e12
        );
        // GM round-trips from AU³/day² to SI base via the core converter.
        assert_relative_eq!(
            p9.gm_typed().value,
            gm_from_au3_day2(p9.gm()).value,
            epsilon = 1e-30
        );
    }

    #[test]
    fn lift_is_non_negative_and_bounded() {
        let p9 = PlanetNine::default();
        let dq = perihelion_lift(&p9, 600.0, 35.0, SECULAR_SCULPTING_DAYS);
        assert!(dq >= 0.0, "lift {dq}");
        // A secular cycle cannot raise q above a (e cannot go negative), and
        // the lift never exceeds the cycle-amplitude cap.
        assert!(35.0 + dq <= 600.0);
    }

    #[test]
    fn forcing_rate_scales_with_p9_mass() {
        // The secular eccentricity-forcing rate is linear in GM₉ -> doubling
        // the mass doubles the rate (the lever the OSSOS XV constraint pulls).
        let a = 500.0;
        let q0 = 35.0;
        let m1 = PlanetNine {
            mass_earth: 5.0,
            ..Default::default()
        };
        let m2 = PlanetNine {
            mass_earth: 10.0,
            ..Default::default()
        };
        let r1 = forcing_rate(&m1, a, q0);
        let r2 = forcing_rate(&m2, a, q0);
        assert!(r1 > 0.0);
        assert_relative_eq!(r2 / r1, 2.0, epsilon = 1e-6);
    }

    #[test]
    fn massive_p9_lifts_more_than_no_p9() {
        // With essentially no perturber the secular torque -> 0 -> no lift.
        let none = PlanetNine {
            mass_earth: 1e-6,
            ..Default::default()
        };
        let real = PlanetNine::default();
        let q0 = 35.0;
        let a = 600.0;
        let t = SECULAR_SCULPTING_DAYS;
        let dq_none = perihelion_lift(&none, a, q0, t);
        let dq_real = perihelion_lift(&real, a, q0, t);
        assert!(dq_none < 0.5, "negligible perturber lift {dq_none}");
        assert!(
            dq_real > dq_none + 1.0,
            "real P9 lift {dq_real} should exceed null {dq_none}"
        );
    }

    #[test]
    fn heavier_p9_lifts_at_least_as_much() {
        // In the rate-limited regime a heavier P9 raises perihelion further
        // within the same sculpting time.
        let a = 500.0;
        let q0 = 36.0;
        let t = SECULAR_SCULPTING_DAYS;
        let light = PlanetNine {
            mass_earth: 2.0,
            ..Default::default()
        };
        let heavy = PlanetNine {
            mass_earth: 8.0,
            ..Default::default()
        };
        let dl = perihelion_lift(&light, a, q0, t);
        let dh = perihelion_lift(&heavy, a, q0, t);
        assert!(dh >= dl, "heavier P9 lift {dh} >= lighter {dl}");
        assert!(dl > 0.0, "even a light P9 should lift somewhat: {dl}");
    }
}
