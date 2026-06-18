//! Integrate the resonant secular pendulum and classify a trajectory as
//! **libration** (Δϖ confined — captured in the secular resonance) or
//! **circulation** (Δϖ sweeps the full circle — no apsidal confinement).
//!
//! The flow is Hamilton's equations for the resonant Hamiltonian
//! `K(Δϖ, Γ) = a₀(Γ) − g₉·Γ + c₁(Γ)·cos Δϖ` (see [`crate::portrait`]):
//!
//! ```text
//!     dΔϖ/dt = ∂K/∂Γ,      dΓ/dt = −∂K/∂Δϖ = c₁(Γ)·sin Δϖ,
//! ```
//!
//! integrated with RK4. `K` is conserved to integrator order, so the
//! trajectory follows a level curve of the resonant Hamiltonian — Beust's
//! secular phase-space portrait. A particle near the resonant action `Γ*` and
//! the stable apse librates; one launched far off-resonance circulates.

use crate::portrait::ResonantHamiltonian;
use p9_core::units::{radians, Angle};
use std::f64::consts::PI;

/// Outcome of following a resonant trajectory.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TrajectoryKind {
    /// Δϖ oscillates about a center without completing a full circle.
    Libration,
    /// Δϖ advances monotonically through 2π (no apsidal confinement).
    Circulation,
}

/// A sampled resonant trajectory in `(Δϖ, Γ)`.
#[derive(Debug, Clone)]
pub struct Trajectory {
    /// Unwrapped Δϖ samples (radians) — continuous, not folded to [0, 2π).
    pub dvarpi: Vec<f64>,
    /// Action Γ samples (AU²/day).
    pub gamma: Vec<f64>,
    /// Conserved resonant-Hamiltonian samples (diagnostics).
    pub k: Vec<f64>,
    /// Energy scale for the drift metric (the pendulum potential amplitude).
    scale: Option<f64>,
    /// Classification.
    pub kind: TrajectoryKind,
}

impl Trajectory {
    /// Libration center: midpoint of the Δϖ swing (radians, folded to [0, 2π)),
    /// as a dimension-checked [`Angle`].
    pub fn center_typed(&self) -> Angle {
        let (lo, hi) = self.bounds();
        radians((0.5 * (lo + hi)).rem_euclid(2.0 * PI))
    }

    /// Libration amplitude: half the peak-to-peak Δϖ swing, as a
    /// dimension-checked [`Angle`].
    pub fn amplitude_typed(&self) -> Angle {
        let (lo, hi) = self.bounds();
        radians(0.5 * (hi - lo))
    }

    fn bounds(&self) -> (f64, f64) {
        let lo = self.dvarpi.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = self
            .dvarpi
            .iter()
            .cloned()
            .fold(f64::NEG_INFINITY, f64::max);
        (lo, hi)
    }

    /// Peak-to-peak excursion of the conserved resonant Hamiltonian, scaled by
    /// the pendulum's potential amplitude (the cos-term swing along the orbit).
    /// A value ≪ 1 means the trajectory faithfully follows a single level curve
    /// of K; the natural scale is the angular forcing the libration rides on,
    /// not the (much larger, dynamically irrelevant) axisymmetric offset.
    pub fn energy_drift(&self) -> f64 {
        let lo = self.k.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = self.k.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        self.scale
            .map(|s| (hi - lo) / s)
            .unwrap_or_else(|| (hi - lo).abs())
    }
}

fn rhs(rh: &ResonantHamiltonian, dvarpi: f64, gamma: f64) -> (f64, f64) {
    (rh.ddvarpi_dt(dvarpi, gamma), rh.dgamma_dt(dvarpi, gamma))
}

/// Integrate the resonant trajectory from `(Δϖ₀, e₀)` for `steps` steps of
/// `dt` days and classify it.
pub fn classify(
    rh: &ResonantHamiltonian,
    dvarpi0: f64,
    e0: f64,
    dt: f64,
    steps: usize,
) -> Trajectory {
    let mut dvarpi = dvarpi0;
    let mut gamma = rh.model.action(e0);

    let gamma_min = rh.model.action(0.02);
    let gamma_max = rh.model.action(0.985);

    // The pendulum potential amplitude sets the energy scale for conservation.
    let scale = rh.c1(gamma).abs();

    let mut traj = Trajectory {
        dvarpi: Vec::with_capacity(steps + 1),
        gamma: Vec::with_capacity(steps + 1),
        k: Vec::with_capacity(steps + 1),
        scale: if scale > 0.0 { Some(scale) } else { None },
        kind: TrajectoryKind::Libration,
    };

    let mut unwrapped = dvarpi0;
    let mut prev_wrapped = dvarpi0.rem_euclid(2.0 * PI);

    traj.dvarpi.push(unwrapped);
    traj.gamma.push(gamma);
    traj.k.push(rh.k(unwrapped, gamma));

    for _ in 0..steps {
        let (k1a, k1g) = rhs(rh, dvarpi, gamma);
        let (k2a, k2g) = rhs(rh, dvarpi + 0.5 * dt * k1a, gamma + 0.5 * dt * k1g);
        let (k3a, k3g) = rhs(rh, dvarpi + 0.5 * dt * k2a, gamma + 0.5 * dt * k2g);
        let (k4a, k4g) = rhs(rh, dvarpi + dt * k3a, gamma + dt * k3g);

        dvarpi += dt / 6.0 * (k1a + 2.0 * k2a + 2.0 * k3a + k4a);
        gamma += dt / 6.0 * (k1g + 2.0 * k2g + 2.0 * k3g + k4g);
        gamma = gamma.clamp(gamma_min, gamma_max);

        let wrapped = dvarpi.rem_euclid(2.0 * PI);
        let mut d = wrapped - prev_wrapped;
        if d > PI {
            d -= 2.0 * PI;
        } else if d < -PI {
            d += 2.0 * PI;
        }
        unwrapped += d;
        prev_wrapped = wrapped;

        traj.dvarpi.push(unwrapped);
        traj.gamma.push(gamma);
        traj.k.push(rh.k(unwrapped, gamma));
    }

    let (lo, hi) = traj.bounds();
    traj.kind = if (hi - lo) >= 2.0 * PI {
        TrajectoryKind::Circulation
    } else {
        TrajectoryKind::Libration
    };
    traj
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hamiltonian::SecularModel;
    use crate::portrait::{libration_island, ResonantHamiltonian};
    use crate::published;
    use approx::assert_relative_eq;

    // A representative distant-ETNO semi-major axis with a robust, well-isolated
    // resonance (2015 RX245 sits at a ≈ 430 AU; 2007 TG422 at ≈ 500 AU).
    const A_TEST: f64 = 400.0;

    /// A resonant (commensurate) particle near the anti-aligned apse librates;
    /// the *same* Planet Nine forcing but with the apse rate detuned away from
    /// commensurability circulates — the apsidal clustering is the signature of
    /// the secular resonance, not of the forcing alone.
    #[test]
    fn resonant_librates_detuned_circulates() {
        let p9 = published::nominal_p9();
        let m = SecularModel::new(A_TEST, &p9);
        let island = libration_island(&m).expect("island should exist");
        let dt = island.suggested_dt;
        let steps = island.suggested_steps;

        // Commensurate: g₉ tuned to the free precession at the resonant e.
        let resonant = ResonantHamiltonian::new(m, island.e_center);
        let lib = classify(
            &resonant,
            island.center_dvarpi + 0.3,
            island.e_center,
            dt,
            steps,
        );
        // Typed angle accessors mirror the inline center/amplitude formulae.
        let lo = lib.dvarpi.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = lib.dvarpi.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        assert_relative_eq!(
            (lib.center_typed() / radians(1.0)).value,
            (0.5 * (lo + hi)).rem_euclid(2.0 * PI),
            max_relative = 1e-12
        );
        let amplitude_rad = (lib.amplitude_typed() / radians(1.0)).value;
        assert_relative_eq!(amplitude_rad, 0.5 * (hi - lo), max_relative = 1e-12);
        assert_eq!(
            lib.kind,
            TrajectoryKind::Libration,
            "commensurate start should librate; amplitude {:.1}°, K-drift {:.1e}",
            amplitude_rad.to_degrees(),
            lib.energy_drift()
        );
        // K is conserved to integrator precision (the pendulum is exact).
        assert!(
            lib.energy_drift() < 1e-3,
            "K not conserved: drift {:.2e}",
            lib.energy_drift()
        );

        // Detuned: push g₉ well above the whole free-precession band → no
        // commensurability → the relative apse circulates everywhere.
        let (_lo, hi) = resonant.free_precession_band();
        let detuned = resonant
            .clone()
            .with_external_rate(3.0 * hi.abs().max(1e-20));
        let circ = classify(&detuned, island.center_dvarpi, island.e_center, dt, steps);
        assert_eq!(
            circ.kind,
            TrajectoryKind::Circulation,
            "detuned (non-commensurate) start should circulate; net swing {:.0}°",
            (2.0 * (circ.amplitude_typed() / radians(1.0)).value).to_degrees()
        );
    }

    /// A large angular displacement (beyond the pendulum separatrix) circulates
    /// while a small one librates — at a *fixed*, commensurate Planet Nine. The
    /// separatrix in the resonant angle is the libration/circulation divide.
    #[test]
    fn small_amplitude_librates_large_circulates() {
        let p9 = published::nominal_p9();
        // Use the narrower resonance at a = 250 AU, where the separatrix sits
        // inside the angle range so a large launch can cross it. (At a = 400 AU
        // the resonance is so wide it spans the whole eccentricity grid — see
        // REPRODUCTION_NOTES; there the circulating regime is only reachable by
        // detuning, exercised in `resonant_librates_detuned_circulates`.)
        let m = SecularModel::new(250.0, &p9);
        let island = libration_island(&m).unwrap();
        let rh = ResonantHamiltonian::new(m, island.e_center);

        let small = classify(
            &rh,
            island.center_dvarpi + 0.2,
            island.e_center,
            island.suggested_dt,
            island.suggested_steps,
        );
        assert_eq!(small.kind, TrajectoryKind::Libration);

        let large = classify(
            &rh,
            island.center_dvarpi + 1.0,
            island.e_center,
            island.suggested_dt,
            island.suggested_steps,
        );
        assert_eq!(
            large.kind,
            TrajectoryKind::Circulation,
            "large launch should cross the separatrix and circulate; swing {:.0}°",
            (2.0 * (large.amplitude_typed() / radians(1.0)).value).to_degrees()
        );
    }

    /// The librating trajectory is centered on the anti-aligned apse (Δϖ ≈ π).
    #[test]
    fn libration_center_is_anti_aligned() {
        let p9 = published::nominal_p9();
        let m = SecularModel::new(A_TEST, &p9);
        let island = libration_island(&m).unwrap();
        let rh = ResonantHamiltonian::new(m, island.e_center);

        let lib = classify(
            &rh,
            island.center_dvarpi + 0.25,
            island.e_center,
            island.suggested_dt,
            island.suggested_steps,
        );
        assert_eq!(lib.kind, TrajectoryKind::Libration);

        let center_rad = (lib.center_typed() / radians(1.0)).value;
        let center_err = (center_rad - published::ANTI_ALIGNED_CENTER_RAD).abs();
        assert!(
            center_err < 0.15,
            "libration center {:.1}° not near 180° (err {:.1}°)",
            center_rad.to_degrees(),
            center_err.to_degrees()
        );
    }

    /// Libration amplitude grows with the launch displacement (pendulum
    /// hallmark), staying anti-aligned.
    #[test]
    fn amplitude_grows_with_displacement() {
        let p9 = published::nominal_p9();
        let m = SecularModel::new(A_TEST, &p9);
        let island = libration_island(&m).unwrap();
        let rh = ResonantHamiltonian::new(m, island.e_center);

        let small = classify(
            &rh,
            island.center_dvarpi + 0.15,
            island.e_center,
            island.suggested_dt,
            island.suggested_steps,
        );
        let large = classify(
            &rh,
            island.center_dvarpi + 0.45,
            island.e_center,
            island.suggested_dt,
            island.suggested_steps,
        );
        assert_eq!(small.kind, TrajectoryKind::Libration);
        assert_eq!(large.kind, TrajectoryKind::Libration);
        let small_amp = (small.amplitude_typed() / radians(1.0)).value;
        let large_amp = (large.amplitude_typed() / radians(1.0)).value;
        assert!(
            large_amp > small_amp,
            "amplitude should grow with displacement: {:.1}° vs {:.1}°",
            large_amp.to_degrees(),
            small_amp.to_degrees()
        );
    }
}
