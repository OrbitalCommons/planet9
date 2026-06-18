//! Secular Hamiltonian for spin-orbit coupling from Bailey+ (2016).
//!
//! Models the interaction between three angular momentum vectors using
//! the vector formulation of quadrupole secular theory:
//!   1. Giant planet angular momentum L_gp (massive wire)
//!   2. Planet Nine orbital angular momentum L_9 (massive wire)
//!   3. Solar spin angular momentum L_sun (test ring)
//!
//! The total angular momentum L_total = L_gp + L_9 is conserved
//! (L_sun is negligible). Each pair interacts through the quadrupole
//! secular Hamiltonian H_ij = C_ij * P₂(L̂_i · L̂_j).
//!
//! The equations of motion in vector form are:
//!   dL̂_i/dt = (3/L_i) Σ_j C_ij (L̂_i · L̂_j) (L̂_i × L̂_j)

use p9_core::analysis::circular::wrap_to_pi;
use p9_core::constants::*;
use p9_core::initial_conditions::giant_planets;
use p9_core::units::{au, days, radians, solar_masses, Angle, Length, Mass, Time};

use crate::solar_model;

/// State of the spin-orbit system as three unit vectors.
///
/// All vectors are in the frame where the total angular momentum
/// L_gp + L_9 initially defines the z-axis.
#[derive(Debug, Clone, Copy, serde::Serialize, serde::Deserialize)]
pub struct SpinOrbitState {
    /// Giant planet angular momentum unit vector
    pub l_gp: [f64; 3],
    /// Planet Nine orbital angular momentum unit vector
    pub l_9: [f64; 3],
    /// Solar spin axis unit vector
    pub l_sun: [f64; 3],
}

impl SpinOrbitState {
    /// Create initial state from inclinations relative to the total angular momentum.
    ///
    /// `i_9` is the inclination of Planet Nine relative to the invariable plane.
    /// The giant planet plane tilts oppositely to conserve total angular momentum.
    /// The solar spin starts aligned with the giant planet plane.
    pub fn from_inclinations(i_9: f64, omega_big_9: f64, l_gp_mag: f64, l_9_mag: f64) -> Self {
        // Total L defines z-axis.
        // L_gp + L_9 = L_total * ẑ
        // In the total-L frame, L_gp tilts by angle i_gp such that:
        //   L_gp * sin(i_gp) = L_9 * sin(i_9) (perpendicular balance)
        let sin_i_gp = (l_9_mag / l_gp_mag) * i_9.sin();
        let i_gp = sin_i_gp.asin();

        // L_gp and L_9 are on opposite sides of the z-axis (in the same plane)
        let l_gp = [
            -i_gp.sin() * omega_big_9.cos(),
            -i_gp.sin() * omega_big_9.sin(),
            i_gp.cos(),
        ];
        let l_9 = [
            i_9.sin() * omega_big_9.cos(),
            i_9.sin() * omega_big_9.sin(),
            i_9.cos(),
        ];

        // Solar spin initially aligned with giant planet plane
        let l_sun = l_gp;

        SpinOrbitState { l_gp, l_9, l_sun }
    }

    /// Angle between solar spin and giant planet angular momentum (solar
    /// obliquity), as a typed [`Angle`].
    pub fn obliquity_typed(&self) -> Angle {
        let dot = self.l_sun[0] * self.l_gp[0]
            + self.l_sun[1] * self.l_gp[1]
            + self.l_sun[2] * self.l_gp[2];
        radians(dot.clamp(-1.0, 1.0).acos())
    }

    /// Mutual inclination between giant planets and Planet Nine, as a typed
    /// [`Angle`].
    pub fn mutual_inclination_typed(&self) -> Angle {
        let dot =
            self.l_gp[0] * self.l_9[0] + self.l_gp[1] * self.l_9[1] + self.l_gp[2] * self.l_9[2];
        radians(dot.clamp(-1.0, 1.0).acos())
    }
}

/// Parameters for the secular integration.
#[derive(Debug, Clone)]
pub struct SecularParams {
    /// Planet Nine mass (solar masses)
    pub m9_solar: f64,
    /// Planet Nine semi-major axis (AU)
    pub a9: f64,
    /// Planet Nine eccentricity
    pub e9: f64,
    /// Total integration time (days)
    pub t_total: f64,
    /// Time step (days)
    pub dt: f64,
}

impl SecularParams {
    pub fn epsilon_9(&self) -> f64 {
        (1.0 - self.e9 * self.e9).sqrt()
    }

    // ---- typed (uom) boundary: dimension-checked views of the f64 fields ----

    /// Planet Nine mass as a [`Mass`].
    pub fn mass(&self) -> Mass {
        solar_masses(self.m9_solar)
    }
    /// Planet Nine semi-major axis as a [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a9)
    }
    /// Total integration time as a [`Time`].
    pub fn total_time(&self) -> Time {
        days(self.t_total)
    }
    /// Integration time step as a [`Time`].
    pub fn timestep(&self) -> Time {
        days(self.dt)
    }
}

/// Quadrupole coupling constant for two wires.
///
///   C = G * m₁ * m₂ / (4 * a_outer) * (a_inner/a_outer)² / ε_outer³
///
/// For circular outer orbit, ε = 1.
fn coupling_constant(m1: f64, m2: f64, a_inner: f64, a_outer: f64, epsilon_outer: f64) -> f64 {
    let g = G_AU3_MSUN_DAY2; // G in AU³/(M_sun * day²)
    let ratio = a_inner / a_outer;
    g * m1 * m2 / (4.0 * a_outer) * ratio * ratio / epsilon_outer.powi(3)
}

/// Giant-planet ↔ Planet Nine coupling, summed planet by planet
/// (∝ Σ mᵢaᵢ²). Collapsing the giants into one L-conserving ring (as a
/// previous version did) underestimates this sum by ~45%, a near-2×
/// systematic in the spin-orbit coupling that shifts the required-i₉
/// contours.
fn coupling_gp_9(m9_solar: f64, a9: f64, epsilon_9: f64) -> f64 {
    giant_planets::GIANT_PLANETS
        .iter()
        .map(|&(m_i, a_i)| coupling_constant(m_i, m9_solar, a_i, a9, epsilon_9))
        .sum()
}

/// Solar-spin-ring ↔ giant-planet coupling, summed planet by planet
/// (∝ Σ mᵢ/aᵢ³; the collapsed ring underestimates this by ~42%).
fn coupling_sun_gp(m_sun_eff: f64, a_tilde: f64) -> f64 {
    giant_planets::GIANT_PLANETS
        .iter()
        .map(|&(m_i, a_i)| coupling_constant(m_sun_eff, m_i, a_tilde, a_i, 1.0))
        .sum()
}

/// Compute the cross product a × b.
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Compute the dot product a · b.
fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Normalize a vector in place.
fn normalize(v: &mut [f64; 3]) {
    let mag = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if mag > 0.0 {
        v[0] /= mag;
        v[1] /= mag;
        v[2] /= mag;
    }
}

/// Compute the derivative of the state vector.
///
/// dL̂_i/dt = (3/L_i) * Σ_j C_ij * (L̂_i · L̂_j) * (L̂_i × L̂_j)
fn derivatives(
    state: &SpinOrbitState,
    l_gp_mag: f64,
    l_9_mag: f64,
    l_sun_mag: f64,
    c_gp_9: f64,
    c_sun_gp: f64,
    c_sun_9: f64,
) -> ([f64; 3], [f64; 3], [f64; 3]) {
    // Giant planet evolution: torqued by P9 and Sun-spin
    let dot_gp_9 = dot(&state.l_gp, &state.l_9);
    let cross_gp_9 = cross(&state.l_gp, &state.l_9);
    let dot_gp_sun = dot(&state.l_gp, &state.l_sun);
    let cross_gp_sun = cross(&state.l_gp, &state.l_sun);

    let mut dl_gp = [0.0; 3];
    for k in 0..3 {
        dl_gp[k] = (3.0 / l_gp_mag)
            * (c_gp_9 * dot_gp_9 * cross_gp_9[k] + c_sun_gp * dot_gp_sun * cross_gp_sun[k]);
    }

    // Planet Nine evolution: torqued by giant planets and Sun-spin
    let cross_9_gp = cross(&state.l_9, &state.l_gp);
    let dot_9_sun = dot(&state.l_9, &state.l_sun);
    let cross_9_sun = cross(&state.l_9, &state.l_sun);

    let mut dl_9 = [0.0; 3];
    for k in 0..3 {
        dl_9[k] = (3.0 / l_9_mag)
            * (c_gp_9 * dot_gp_9 * cross_9_gp[k] + c_sun_9 * dot_9_sun * cross_9_sun[k]);
    }

    // Solar spin evolution: torqued by giant planets and P9
    let cross_sun_gp = cross(&state.l_sun, &state.l_gp);
    let dot_sun_9 = dot(&state.l_sun, &state.l_9);
    let cross_sun_9 = cross(&state.l_sun, &state.l_9);

    let mut dl_sun = [0.0; 3];
    for k in 0..3 {
        dl_sun[k] = (3.0 / l_sun_mag)
            * (c_sun_gp * dot_gp_sun * cross_sun_gp[k] + c_sun_9 * dot_sun_9 * cross_sun_9[k]);
    }

    (dl_gp, dl_9, dl_sun)
}

/// Add scaled derivative: dst = src + scale * deriv.
///
/// RK4 stages are *not* renormalized: projecting the intermediate stages
/// onto the unit sphere violates the Butcher conditions and degrades the
/// method to 2nd order. Only the completed step is renormalized (an O(dt⁵)
/// projection, since the exact flow preserves the norms), which keeps the
/// scheme 4th order — see the dt-halving convergence test.
fn add_scaled(
    src: &SpinOrbitState,
    scale: f64,
    dl_gp: &[f64; 3],
    dl_9: &[f64; 3],
    dl_sun: &[f64; 3],
) -> SpinOrbitState {
    SpinOrbitState {
        l_gp: [
            src.l_gp[0] + scale * dl_gp[0],
            src.l_gp[1] + scale * dl_gp[1],
            src.l_gp[2] + scale * dl_gp[2],
        ],
        l_9: [
            src.l_9[0] + scale * dl_9[0],
            src.l_9[1] + scale * dl_9[1],
            src.l_9[2] + scale * dl_9[2],
        ],
        l_sun: [
            src.l_sun[0] + scale * dl_sun[0],
            src.l_sun[1] + scale * dl_sun[1],
            src.l_sun[2] + scale * dl_sun[2],
        ],
    }
}

/// Snapshot of the integration at a given time.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ObliquitySnapshot {
    /// Time in days
    pub t: f64,
    /// Solar obliquity (angle between spin and giant planet plane) in radians
    pub obliquity: f64,
    /// Mutual inclination between giant planets and P9 (rad)
    pub mutual_inclination: f64,
    /// Solar spin node longitude (rad)
    pub omega_big_sun: f64,
    /// Planet Nine inclination (rad)
    pub i_9: f64,
    /// Planet Nine node longitude (rad)
    pub omega_big_9: f64,
    /// ΔΩ = Ω_sun - Ω_9 (rad)
    pub delta_omega_big: f64,
}

impl ObliquitySnapshot {
    // ---- typed (uom) boundary: dimension-checked views of the f64 fields ----

    /// Snapshot time as a [`Time`].
    pub fn time(&self) -> Time {
        days(self.t)
    }
    /// Solar obliquity as an [`Angle`].
    pub fn obliquity_typed(&self) -> Angle {
        radians(self.obliquity)
    }
    /// Mutual inclination (giant planets vs P9) as an [`Angle`].
    pub fn mutual_inclination_typed(&self) -> Angle {
        radians(self.mutual_inclination)
    }
    /// Solar spin node longitude as an [`Angle`].
    pub fn omega_big_sun_typed(&self) -> Angle {
        radians(self.omega_big_sun)
    }
    /// Planet Nine inclination as an [`Angle`].
    pub fn i9_typed(&self) -> Angle {
        radians(self.i_9)
    }
    /// Planet Nine node longitude as an [`Angle`].
    pub fn omega_big_9_typed(&self) -> Angle {
        radians(self.omega_big_9)
    }
    /// ΔΩ = Ω_sun − Ω_9 as an [`Angle`].
    pub fn delta_omega_typed(&self) -> Angle {
        radians(self.delta_omega_big)
    }
}

/// Integrate the secular spin-orbit equations over the solar system lifetime.
///
/// Uses RK4 with time-dependent coupling (Skumanich spin-down).
pub fn integrate_obliquity(
    initial: SpinOrbitState,
    params: &SecularParams,
    snapshot_interval: f64,
) -> Vec<ObliquitySnapshot> {
    let t_age = params.t_total;
    let dt = params.dt;
    let n_steps = (t_age / dt).ceil() as usize;

    let l_gp_mag = giant_planets::giant_planet_angular_momentum();
    let epsilon_9 = params.epsilon_9();
    let l_9_mag = giant_planets::p9_angular_momentum(params.m9_solar, params.a9, params.e9);

    // Giant planet ↔ Planet Nine coupling (constant, per-planet sum)
    let c_gp_9 = coupling_gp_9(params.m9_solar, params.a9, epsilon_9);

    let mut state = initial;
    let mut t = 0.0;
    let mut snapshots = Vec::new();
    let mut next_snapshot = 0.0;

    snapshots.push(make_snapshot(&state, t));

    for _ in 0..n_steps {
        // Time-dependent solar parameters
        let compute_solar_couplings = |t_now: f64| {
            // ω in rad/day, extracted version-safely from the typed quantity.
            let omega_rot = (solar_model::solar_omega_at_time_typed(t_now, t_age)
                / (radians(1.0) / days(1.0)))
            .value;
            let a_tilde = (solar_model::solar_ring_semimajor_axis_typed(omega_rot) / au(1.0)).value;
            let l_sun_mag = solar_model::solar_spin_angular_momentum(omega_rot);

            // Effective mass of solar spin ring
            let m_sun_eff = l_sun_mag / (G_AU3_MSUN_DAY2 * a_tilde).sqrt();

            let c_sun_gp = coupling_sun_gp(m_sun_eff, a_tilde);
            let c_sun_9 =
                coupling_constant(m_sun_eff, params.m9_solar, a_tilde, params.a9, epsilon_9);

            (l_sun_mag, c_sun_gp, c_sun_9)
        };

        // RK4 with time-dependent coefficients
        let (l_sun_1, c_sg_1, c_s9_1) = compute_solar_couplings(t);
        let (k1_gp, k1_9, k1_sun) =
            derivatives(&state, l_gp_mag, l_9_mag, l_sun_1, c_gp_9, c_sg_1, c_s9_1);

        let s2 = add_scaled(&state, 0.5 * dt, &k1_gp, &k1_9, &k1_sun);
        let (l_sun_2, c_sg_2, c_s9_2) = compute_solar_couplings(t + 0.5 * dt);
        let (k2_gp, k2_9, k2_sun) =
            derivatives(&s2, l_gp_mag, l_9_mag, l_sun_2, c_gp_9, c_sg_2, c_s9_2);

        let s3 = add_scaled(&state, 0.5 * dt, &k2_gp, &k2_9, &k2_sun);
        let (l_sun_3, c_sg_3, c_s9_3) = compute_solar_couplings(t + 0.5 * dt);
        let (k3_gp, k3_9, k3_sun) =
            derivatives(&s3, l_gp_mag, l_9_mag, l_sun_3, c_gp_9, c_sg_3, c_s9_3);

        let s4 = add_scaled(&state, dt, &k3_gp, &k3_9, &k3_sun);
        let (l_sun_4, c_sg_4, c_s9_4) = compute_solar_couplings(t + dt);
        let (k4_gp, k4_9, k4_sun) =
            derivatives(&s4, l_gp_mag, l_9_mag, l_sun_4, c_gp_9, c_sg_4, c_s9_4);

        // RK4 update
        for k in 0..3 {
            state.l_gp[k] += dt / 6.0 * (k1_gp[k] + 2.0 * k2_gp[k] + 2.0 * k3_gp[k] + k4_gp[k]);
            state.l_9[k] += dt / 6.0 * (k1_9[k] + 2.0 * k2_9[k] + 2.0 * k3_9[k] + k4_9[k]);
            state.l_sun[k] +=
                dt / 6.0 * (k1_sun[k] + 2.0 * k2_sun[k] + 2.0 * k3_sun[k] + k4_sun[k]);
        }
        normalize(&mut state.l_gp);
        normalize(&mut state.l_9);
        normalize(&mut state.l_sun);

        t += dt;

        if t >= next_snapshot + snapshot_interval {
            snapshots.push(make_snapshot(&state, t));
            next_snapshot += snapshot_interval;
        }
    }

    snapshots.push(make_snapshot(&state, t));
    snapshots
}

fn make_snapshot(state: &SpinOrbitState, t: f64) -> ObliquitySnapshot {
    let obliquity = (state.obliquity_typed() / radians(1.0)).value;
    let mutual_inclination = (state.mutual_inclination_typed() / radians(1.0)).value;

    // Extract angles from unit vectors for diagnostics
    let i_9 = state.l_9[2].clamp(-1.0, 1.0).acos();
    let omega_big_9 = state.l_9[1].atan2(state.l_9[0]);

    let _i_sun = state.l_sun[2].clamp(-1.0, 1.0).acos();
    let omega_big_sun = state.l_sun[1].atan2(state.l_sun[0]);

    let delta = wrap_to_pi(omega_big_sun - omega_big_9);

    ObliquitySnapshot {
        t,
        obliquity,
        mutual_inclination,
        omega_big_sun,
        i_9,
        omega_big_9,
        delta_omega_big: delta,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use std::f64::consts::PI;
    use uom::si::angle::radian;
    use uom::si::length::astronomical_unit;
    use uom::si::time::day;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let params = SecularParams {
            m9_solar: 10.0 * EARTH_MASS_SOLAR,
            a9: 700.0,
            e9: 0.6,
            t_total: 4.5 * GYR_DAYS,
            dt: 1e5 * YEAR_DAYS,
        };
        assert_relative_eq!(
            params.semi_major_axis().get::<astronomical_unit>(),
            params.a9
        );
        assert_relative_eq!(params.total_time().get::<day>(), params.t_total);
        assert_relative_eq!(params.timestep().get::<day>(), params.dt);

        let state = SpinOrbitState::from_inclinations(0.5, 1.0, 1.0, 0.3);
        // obliquity = acos(L̂_sun · L̂_gp).
        let dot_sun_gp = state.l_sun[0] * state.l_gp[0]
            + state.l_sun[1] * state.l_gp[1]
            + state.l_sun[2] * state.l_gp[2];
        assert_relative_eq!(
            state.obliquity_typed().get::<radian>(),
            dot_sun_gp.clamp(-1.0, 1.0).acos()
        );
        // mutual inclination = acos(L̂_gp · L̂_9).
        let dot_gp_9 = state.l_gp[0] * state.l_9[0]
            + state.l_gp[1] * state.l_9[1]
            + state.l_gp[2] * state.l_9[2];
        assert_relative_eq!(
            state.mutual_inclination_typed().get::<radian>(),
            dot_gp_9.clamp(-1.0, 1.0).acos()
        );

        let snap = ObliquitySnapshot {
            t: 1.0e9,
            obliquity: 0.1,
            mutual_inclination: 0.2,
            omega_big_sun: 0.3,
            i_9: 0.4,
            omega_big_9: 0.5,
            delta_omega_big: -0.2,
        };
        assert_relative_eq!(snap.time().get::<day>(), snap.t);
        assert_relative_eq!(snap.obliquity_typed().get::<radian>(), snap.obliquity);
        assert_relative_eq!(
            snap.delta_omega_typed().get::<radian>(),
            snap.delta_omega_big
        );
    }

    #[test]
    fn test_coupling_constant_positive() {
        let c = coupling_constant(1e-3, 3e-5, 10.0, 700.0, 0.8);
        assert!(c > 0.0, "Coupling constant should be positive");
        assert!(c.is_finite());
    }

    #[test]
    fn test_per_planet_coupling_exceeds_collapsed_ring() {
        // H6 regression: the collapsed L-conserving ring (m_total at a_eff
        // with L = m_total √(GM a_eff)) underestimates both coupling sums
        // by ~40-45%; the per-planet sums must be substantially larger.
        let m_total: f64 = giant_planets::GIANT_PLANETS.iter().map(|&(m, _)| m).sum();
        let l_total = giant_planets::giant_planet_angular_momentum();
        let a_eff = (l_total / m_total).powi(2) / G_AU3_MSUN_DAY2;

        let m9 = 10.0 * EARTH_MASS_SOLAR;
        let (a9, eps9) = (700.0, 0.8);
        let c_collapsed_gp9 = coupling_constant(m_total, m9, a_eff, a9, eps9);
        let c_per_planet_gp9 = coupling_gp_9(m9, a9, eps9);
        let ratio_gp9 = c_per_planet_gp9 / c_collapsed_gp9;
        assert!(
            ratio_gp9 > 1.5 && ratio_gp9 < 2.5,
            "Σmᵢaᵢ² per-planet / collapsed = {ratio_gp9:.2} (expected ~1.8)"
        );

        let a_tilde = 1e-4; // a_tilde << a_planet: ratio independent of it
        let m_sun_eff = 1e-2;
        let c_collapsed_sgp = coupling_constant(m_sun_eff, m_total, a_tilde, a_eff, 1.0);
        let c_per_planet_sgp = coupling_sun_gp(m_sun_eff, a_tilde);
        let ratio_sgp = c_per_planet_sgp / c_collapsed_sgp;
        assert!(
            ratio_sgp > 1.4 && ratio_sgp < 2.5,
            "Σmᵢ/aᵢ³ per-planet / collapsed = {ratio_sgp:.2} (expected ~1.7)"
        );
    }

    #[test]
    fn test_node_regression_matches_two_vector_frequency() {
        // L4: in the two-vector (gp ↔ 9) problem both L̂'s precess rigidly
        // about the conserved total L with frequency
        //   Ω̇ = −3 C cos θ |L_tot| / (L_gp L_9)
        // (regression for prograde mutual inclination). The solar spin ring
        // contributes only at the ~1e-10 level, so the integrated node rate
        // must match the analytic frequency.
        let m9_solar = 10.0 * EARTH_MASS_SOLAR;
        let (a9, e9) = (700.0, 0.6);
        let l_gp_mag = giant_planets::giant_planet_angular_momentum();
        let l_9_mag = giant_planets::p9_angular_momentum(m9_solar, a9, e9);

        let initial = SpinOrbitState::from_inclinations(20.0 * DEG2RAD, PI, l_gp_mag, l_9_mag);
        let theta = (initial.mutual_inclination_typed() / radians(1.0)).value;
        let epsilon_9 = (1.0_f64 - e9 * e9).sqrt();
        let c_gp_9 = coupling_gp_9(m9_solar, a9, epsilon_9);

        // |L_tot| with the two vectors θ apart
        let l_tot =
            (l_gp_mag * l_gp_mag + l_9_mag * l_9_mag + 2.0 * l_gp_mag * l_9_mag * theta.cos())
                .sqrt();
        let omega_node_analytic = -3.0 * c_gp_9 * theta.cos() * l_tot / (l_gp_mag * l_9_mag);

        let params = SecularParams {
            m9_solar,
            a9,
            e9,
            t_total: 1e8 * YEAR_DAYS,
            dt: 1e4 * YEAR_DAYS,
        };
        let snapshots = integrate_obliquity(initial, &params, 1e7 * YEAR_DAYS);

        let first = &snapshots[0];
        let last = snapshots.last().unwrap();
        // Total node motion is ~0.02 rad here, far from wrapping.
        let omega_node_measured = (last.omega_big_9 - first.omega_big_9) / (last.t - first.t);

        assert!(
            omega_node_measured < 0.0,
            "prograde mutual tilt must regress the node: {omega_node_measured:.3e}"
        );
        let rel = ((omega_node_measured - omega_node_analytic) / omega_node_analytic).abs();
        assert!(
            rel < 0.05,
            "node rate {omega_node_measured:.4e} vs analytic {omega_node_analytic:.4e} (rel {rel:.2e})"
        );
    }

    #[test]
    fn test_rk4_dt_halving_convergence() {
        // L5: with stage renormalization removed, the integrator must show
        // 4th-order convergence (error ratio ~16 under dt halving). Uses a
        // compact toy orbit (a9 = 70 AU) so the precession is fast enough
        // for truncation error to dominate roundoff.
        let m9_solar = 10.0 * EARTH_MASS_SOLAR;
        let (a9, e9) = (70.0, 0.6);
        let l_gp_mag = giant_planets::giant_planet_angular_momentum();
        let l_9_mag = giant_planets::p9_angular_momentum(m9_solar, a9, e9);
        let initial = SpinOrbitState::from_inclinations(30.0 * DEG2RAD, PI, l_gp_mag, l_9_mag);

        let t_total = 1e7 * YEAR_DAYS;
        let run = |dt: f64| -> [f64; 3] {
            let params = SecularParams {
                m9_solar,
                a9,
                e9,
                t_total,
                dt,
            };
            let snap = integrate_obliquity(initial, &params, t_total);
            let last = snap.last().unwrap();
            [
                last.i_9.sin() * last.omega_big_9.cos(),
                last.i_9.sin() * last.omega_big_9.sin(),
                last.i_9.cos(),
            ]
        };

        let dist = |a: [f64; 3], b: [f64; 3]| -> f64 {
            ((a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)).sqrt()
        };

        let reference = run(1.25e4 * YEAR_DAYS);
        let e1 = dist(run(1e5 * YEAR_DAYS), reference);
        let e2 = dist(run(5e4 * YEAR_DAYS), reference);

        assert!(
            e1 > 1e-12,
            "error too small to measure convergence: {e1:.2e}"
        );
        let ratio = e1 / e2;
        assert!(
            (8.0..32.0).contains(&ratio),
            "dt-halving error ratio = {ratio:.1} (4th order expects ~16; e1={e1:.2e}, e2={e2:.2e})"
        );
    }

    #[test]
    fn test_no_obliquity_without_p9_inclination() {
        // With i_9 = 0, all vectors are aligned → no obliquity excitation
        let l_gp_mag = giant_planets::giant_planet_angular_momentum();
        let m9_solar = 10.0 * EARTH_MASS_SOLAR;
        let l_9_mag = giant_planets::p9_angular_momentum(m9_solar, 700.0, 0.6);

        let initial = SpinOrbitState::from_inclinations(0.001, PI, l_gp_mag, l_9_mag);

        let params = SecularParams {
            m9_solar,
            a9: 700.0,
            e9: 0.6,
            t_total: 1e8 * YEAR_DAYS,
            dt: 1e5 * YEAR_DAYS,
        };

        let snapshots = integrate_obliquity(initial, &params, 1e7 * YEAR_DAYS);
        let final_snap = snapshots.last().unwrap();

        assert!(
            final_snap.obliquity < 0.1,
            "Obliquity grew to {:.4} rad without inclined P9",
            final_snap.obliquity
        );
    }

    #[test]
    fn test_obliquity_excited_by_inclined_p9() {
        let m9_solar = 15.0 * EARTH_MASS_SOLAR;
        let l_gp_mag = giant_planets::giant_planet_angular_momentum();
        let l_9_mag = giant_planets::p9_angular_momentum(m9_solar, 500.0, 0.5);

        let initial = SpinOrbitState::from_inclinations(30.0 * DEG2RAD, PI, l_gp_mag, l_9_mag);

        let params = SecularParams {
            m9_solar,
            a9: 500.0,
            e9: 0.5,
            t_total: 4.5 * GYR_DAYS,
            dt: 1e5 * YEAR_DAYS,
        };

        let snapshots = integrate_obliquity(initial, &params, 0.5 * GYR_DAYS);
        let final_snap = snapshots.last().unwrap();
        let obliquity_deg = final_snap.obliquity * RAD2DEG;

        // The paper shows ~6° for optimal parameters. We should get significant excitation.
        assert!(
            obliquity_deg > 1.0,
            "Obliquity should be excited above 1°, got {:.2}°",
            obliquity_deg
        );
    }

    #[test]
    fn test_mutual_inclination_roughly_conserved() {
        let m9_solar = 10.0 * EARTH_MASS_SOLAR;
        let l_gp_mag = giant_planets::giant_planet_angular_momentum();
        let l_9_mag = giant_planets::p9_angular_momentum(m9_solar, 700.0, 0.6);

        let initial = SpinOrbitState::from_inclinations(20.0 * DEG2RAD, PI, l_gp_mag, l_9_mag);

        let params = SecularParams {
            m9_solar,
            a9: 700.0,
            e9: 0.6,
            t_total: 1e8 * YEAR_DAYS,
            dt: 1e5 * YEAR_DAYS,
        };

        let snapshots = integrate_obliquity(initial, &params, 1e7 * YEAR_DAYS);

        let initial_mi = snapshots[0].mutual_inclination;
        for snap in &snapshots {
            assert!(snap.obliquity.is_finite());
            assert!(snap.mutual_inclination.is_finite());
            // Mutual inclination should be roughly conserved (P9↔gp dominates)
            let diff_deg = (snap.mutual_inclination - initial_mi).abs() * RAD2DEG;
            assert!(
                diff_deg < 5.0,
                "Mutual inclination changed by {:.2}° (too much)",
                diff_deg
            );
        }
    }

    #[test]
    fn test_integration_produces_snapshots() {
        let m9_solar = 10.0 * EARTH_MASS_SOLAR;
        let l_gp_mag = giant_planets::giant_planet_angular_momentum();
        let l_9_mag = giant_planets::p9_angular_momentum(m9_solar, 700.0, 0.6);

        let initial = SpinOrbitState::from_inclinations(20.0 * DEG2RAD, PI, l_gp_mag, l_9_mag);

        let params = SecularParams {
            m9_solar,
            a9: 700.0,
            e9: 0.6,
            t_total: 1e8 * YEAR_DAYS,
            dt: 1e5 * YEAR_DAYS,
        };

        let snapshots = integrate_obliquity(initial, &params, 1e7 * YEAR_DAYS);

        assert!(snapshots.len() >= 2, "Should have multiple snapshots");

        for snap in &snapshots {
            assert!(snap.obliquity.is_finite());
            assert!(snap.i_9.is_finite());
        }
    }
}
