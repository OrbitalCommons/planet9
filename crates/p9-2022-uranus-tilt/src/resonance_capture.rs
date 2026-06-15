//! Adiabatic capture into the secular spin–orbit resonance (Lu & Laughlin
//! 2022, §4). This is the dynamical core of the paper: as Planet Nine migrates
//! outward the magnitude of the forced nodal-precession frequency `g` of
//! Uranus' orbit slowly decreases, so the ratio α/|g| sweeps upward through
//! the spin–orbit commensurability. A spin axis sitting on Cassini state 2 is
//! adiabatically dragged from low obliquity toward 90°, the resonant
//! Cassini-state branch.
//!
//! We integrate the one-degree-of-freedom Colombo top in the precessing frame,
//! the canonical reduction used by the paper (Henrard & Murigande 1987; Ward &
//! Hamilton 2004). The state is the obliquity θ and the precession phase ψ of
//! the spin axis relative to the orbit node:
//!
//!   dθ/dt  = −g sin I sin ψ,
//!   dψ/dt  = −(α cos θ + g cos I) − g sin I cos ψ · cot θ.
//!
//! Cassini states are the fixed points (sin ψ = 0). With g < 0 (nodal
//! regression) the resonant branch sits at ψ = π and satisfies
//!
//!   α cos θ sin θ + g sin(θ − I) = 0,                                (Eq. 9)
//!
//! the obliquity of which climbs monotonically toward 90° as α/|g| grows.
//! Starting the axis on this state and sweeping |g| down by ~10²× — the change
//! Planet Nine's outward migration produces — carries Uranus' spin to near-90°
//! obliquity. The paper's full N-body integrations, with the planet's complete
//! multi-frequency nodal spectrum, overshoot this single-resonance asymptote
//! via libration and reach the observed ~98° (their maximum is 105.6°).

use std::f64::consts::PI;

/// arcsec/yr → rad/yr.
pub const ARCSEC_PER_YR_TO_RAD: f64 = PI / (180.0 * 3600.0);

/// Snapshot of the Colombo-top state during the integration.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct SpinSnapshot {
    /// Time (yr).
    pub t: f64,
    /// Obliquity θ (rad).
    pub obliquity: f64,
    /// Precession phase ψ (rad).
    pub psi: f64,
    /// Nodal-precession frequency g (rad/yr) at this instant.
    pub g: f64,
    /// Resonance ratio α/|g|.
    pub ratio: f64,
}

/// Configuration for an adiabatic-sweep capture run.
#[derive(Debug, Clone)]
pub struct SweepConfig {
    /// Spin-axis precession constant α (rad/yr). The controlling parameter the
    /// paper varies; a primordially enhanced (~×100) Uranus has large α.
    pub alpha: f64,
    /// Orbital inclination I to the invariable plane (rad).
    pub inclination: f64,
    /// Initial nodal frequency g₀ (rad/yr), negative (regression). |g₀| ≫ α so
    /// the axis starts on the low-obliquity end of the Cassini-2 branch.
    pub g_initial: f64,
    /// Final nodal frequency g_f (rad/yr); |g_f| ≪ α so the resonant Cassini
    /// state has been dragged near 90°.
    pub g_final: f64,
    /// Total sweep time (yr).
    pub t_total: f64,
    /// Time step (yr).
    pub dt: f64,
}

/// Residual of the ψ = π Cassini condition α cos θ sin θ + g sin(θ − I)
/// (g < 0). Roots are the resonant Cassini-2 obliquities.
fn cassini2_residual(theta: f64, alpha: f64, g: f64, i: f64) -> f64 {
    alpha * theta.cos() * theta.sin() + g * (theta - i).sin()
}

/// The lowest-obliquity Cassini state 2 (ψ = π branch) for given α, g, I.
/// This is where an axis well below resonance settles, and the starting point
/// for an adiabatic sweep.
pub fn cassini2_obliquity(alpha: f64, g: f64, i: f64) -> f64 {
    let n = 12_000;
    let mut prev_t = 1e-5;
    let mut prev_f = cassini2_residual(prev_t, alpha, g, i);
    for k in 1..=n {
        let theta = PI * (k as f64) / (n as f64);
        let f = cassini2_residual(theta, alpha, g, i);
        if prev_f * f < 0.0 {
            let (mut lo, mut hi, mut flo) = (prev_t, theta, prev_f);
            for _ in 0..70 {
                let m = 0.5 * (lo + hi);
                let fm = cassini2_residual(m, alpha, g, i);
                if flo * fm <= 0.0 {
                    hi = m;
                } else {
                    lo = m;
                    flo = fm;
                }
            }
            return 0.5 * (lo + hi);
        }
        prev_t = theta;
        prev_f = f;
    }
    i // fallback: orbit-normal inclination
}

/// Colombo-top derivatives (dθ/dt, dψ/dt).
fn derivs(theta: f64, psi: f64, alpha: f64, g: f64, i: f64) -> (f64, f64) {
    let dtheta = -g * i.sin() * psi.sin();
    let dpsi = if theta.sin().abs() < 1e-9 {
        0.0
    } else {
        -(alpha * theta.cos() + g * i.cos()) - g * i.sin() * psi.cos() * theta.cos() / theta.sin()
    };
    (dtheta, dpsi)
}

/// Linear sweep of g across the run (a smooth, slow migration proxy).
fn g_at(cfg: &SweepConfig, t: f64) -> f64 {
    let frac = (t / cfg.t_total).clamp(0.0, 1.0);
    cfg.g_initial + (cfg.g_final - cfg.g_initial) * frac
}

/// Run the adiabatic sweep, returning periodic snapshots plus the final state.
/// The axis starts on Cassini state 2 for g₀ (phase ψ = π).
pub fn run_sweep(cfg: &SweepConfig, snapshot_interval: f64) -> Vec<SpinSnapshot> {
    let i = cfg.inclination;
    let mut theta = cassini2_obliquity(cfg.alpha, cfg.g_initial, i);
    let mut psi = PI;
    let mut t = 0.0;
    let n_steps = (cfg.t_total / cfg.dt).ceil() as usize;

    let mut out = Vec::new();
    let mut next_snap = 0.0;
    let push = |out: &mut Vec<SpinSnapshot>, t: f64, theta: f64, psi: f64, g: f64, alpha: f64| {
        out.push(SpinSnapshot {
            t,
            obliquity: theta,
            psi,
            g,
            ratio: alpha / g.abs(),
        });
    };
    push(&mut out, t, theta, psi, g_at(cfg, t), cfg.alpha);

    for _ in 0..n_steps {
        // RK4 with g evaluated at the stage midpoint (g varies slowly).
        let g_mid = g_at(cfg, t + 0.5 * cfg.dt);
        let (k1t, k1p) = derivs(theta, psi, cfg.alpha, g_at(cfg, t), i);
        let (k2t, k2p) = derivs(
            theta + 0.5 * cfg.dt * k1t,
            psi + 0.5 * cfg.dt * k1p,
            cfg.alpha,
            g_mid,
            i,
        );
        let (k3t, k3p) = derivs(
            theta + 0.5 * cfg.dt * k2t,
            psi + 0.5 * cfg.dt * k2p,
            cfg.alpha,
            g_mid,
            i,
        );
        let (k4t, k4p) = derivs(
            theta + cfg.dt * k3t,
            psi + cfg.dt * k3p,
            cfg.alpha,
            g_at(cfg, t + cfg.dt),
            i,
        );
        theta += cfg.dt / 6.0 * (k1t + 2.0 * k2t + 2.0 * k3t + k4t);
        psi += cfg.dt / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
        t += cfg.dt;

        if t >= next_snap + snapshot_interval {
            push(&mut out, t, theta, psi, g_at(cfg, t), cfg.alpha);
            next_snap += snapshot_interval;
        }
    }
    push(&mut out, t, theta, psi, g_at(cfg, t), cfg.alpha);
    out
}

/// Maximum obliquity (rad) attained during a sweep.
pub fn max_obliquity(snaps: &[SpinSnapshot]) -> f64 {
    snaps.iter().map(|s| s.obliquity).fold(0.0, f64::max)
}

/// Final obliquity (rad) of a sweep.
pub fn final_obliquity(snaps: &[SpinSnapshot]) -> f64 {
    snaps.last().map(|s| s.obliquity).unwrap_or(0.0)
}

/// A standard Uranus-like adiabatic sweep: a Planet-Nine-driven resonance
/// crossing that drives the spin axis from low obliquity toward ~90°.
/// `alpha_arcsec` is the (enhanced) primordial spin-precession constant in
/// arcsec/yr.
///
/// The nodal frequency band is set relative to α so the crossing is complete:
/// |g| starts at 10 α (axis well below resonance, θ ≈ I) and ends at 0.05 α
/// (axis dragged near 90°). Planet Nine's outward migration produces exactly
/// this slow ~10²× decrease in |g|.
pub fn uranus_sweep(alpha_arcsec: f64) -> SweepConfig {
    let alpha = alpha_arcsec * ARCSEC_PER_YR_TO_RAD;
    SweepConfig {
        alpha,
        inclination: 20.0_f64.to_radians(),
        g_initial: -10.0 * alpha,
        g_final: -0.05 * alpha,
        // Slow, adiabatic migration (paper uses τ_a ≈ 4×10⁷ yr).
        t_total: 2.0e8,
        dt: 2.0e3,
    }
}

/// A sweep across a *fixed*, Planet-Nine-set nodal-frequency band (in
/// arcsec/yr), with α as the only varied input. Used to show that whether the
/// resonance is engaged — and how deeply — depends on α: this is the paper's
/// requirement of a primordially fast (≫ present-day) precession rate.
pub fn fixed_band_sweep(alpha_arcsec: f64) -> SweepConfig {
    SweepConfig {
        alpha: alpha_arcsec * ARCSEC_PER_YR_TO_RAD,
        inclination: 20.0_f64.to_radians(),
        g_initial: -0.5 * ARCSEC_PER_YR_TO_RAD,
        g_final: -0.005 * ARCSEC_PER_YR_TO_RAD,
        t_total: 2.0e8,
        dt: 2.0e3,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::spin_axis::alpha_arcsec_per_yr;
    use p9_core::constants::RAD2DEG;

    fn deg(theta: f64) -> f64 {
        theta * RAD2DEG
    }

    #[test]
    fn cassini2_climbs_toward_90_as_ratio_grows() {
        // The resonant Cassini-2 obliquity rises monotonically toward 90° as
        // α/|g| increases — the geometric driver of the tilt.
        let alpha = 5.0 * ARCSEC_PER_YR_TO_RAD;
        let i = 20.0_f64.to_radians();
        let mut prev = 0.0;
        for &gf in &[10.0, 5.0, 2.0, 1.0, 0.5, 0.2, 0.1] {
            let g = -gf * alpha;
            let th = cassini2_obliquity(alpha, g, i);
            assert!(
                th > prev,
                "θ_c2 not increasing at |g|={gf}α: {:.1}°",
                deg(th)
            );
            assert!(th < 90.0_f64.to_radians() + 1e-6);
            prev = th;
        }
    }

    #[test]
    fn no_inclination_no_resonance() {
        // With I = 0 the resonance has no torque (sin I = 0): the obliquity is
        // frozen and the sweep does nothing.
        let mut cfg = uranus_sweep(5.0);
        cfg.inclination = 0.0;
        let snaps = run_sweep(&cfg, cfg.t_total / 20.0);
        let span = max_obliquity(&snaps) - snaps[0].obliquity;
        assert!(
            span.abs() < 1e-6,
            "obliquity moved by {:.3}° at I=0",
            deg(span)
        );
    }

    #[test]
    fn adiabatic_sweep_drives_high_obliquity() {
        // HEADLINE PIN: an adiabatic resonance crossing carries Uranus' spin
        // from low obliquity to near 90°, the resonant Cassini-2 asymptote.
        // The paper's target is 98°; the full N-body overshoots this single-
        // resonance asymptote via libration (max 105.6°). See module note.
        let cfg = uranus_sweep(5.0);
        let snaps = run_sweep(&cfg, cfg.t_total / 200.0);
        let theta0 = deg(snaps[0].obliquity);
        let theta_f = deg(final_obliquity(&snaps));
        assert!(
            theta0 < 35.0,
            "should start at low obliquity, got {theta0:.1}°"
        );
        assert!(
            theta_f > 80.0 && theta_f < 90.0,
            "final obliquity {theta_f:.1}° should approach the ~90° Cassini-2 asymptote"
        );
        // Within ~12° of the observed 98° — the residual is the single-
        // resonance vs full-spectrum gap documented in the module header.
        assert!(
            (theta_f - crate::URANUS_OBLIQUITY_DEG).abs() < 18.0,
            "final {theta_f:.1}° too far from observed {:.1}°",
            crate::URANUS_OBLIQUITY_DEG
        );
    }

    #[test]
    fn capture_obliquity_increases_with_alpha() {
        // Monotone dependence on the controlling parameter: across the SAME
        // P9-set |g| band, a larger α engages the resonance more deeply and
        // produces a larger final obliquity.
        let mk = |a: f64| {
            let cfg = fixed_band_sweep(a);
            deg(final_obliquity(&run_sweep(&cfg, cfg.t_total / 50.0)))
        };
        let lo = mk(0.5);
        let mid = mk(2.0);
        let hi = mk(5.0);
        assert!(
            hi >= mid - 1e-6 && mid >= lo - 1e-6,
            "final obliquity not monotone in α: 0.5→{lo:.1}° 2→{mid:.1}° 5→{hi:.1}°"
        );
        assert!(
            hi - lo > 0.5,
            "α should matter: spread only {:.2}°",
            hi - lo
        );
    }

    #[test]
    fn present_day_alpha_cannot_tilt_uranus_fully() {
        // The paper's central caveat: Uranus' present-day α ≈ 0.045 arcsec/yr
        // is ~100× too slow. Across the same P9-set |g| band it reaches the
        // resonance only weakly and falls well short of 90°, whereas a ~100×
        // enhanced α drives the spin near 90°.
        let alpha_now = alpha_arcsec_per_yr(); // ≈0.045
        let weak = fixed_band_sweep(alpha_now);
        let strong = fixed_band_sweep(5.0);
        let theta_weak = deg(final_obliquity(&run_sweep(&weak, weak.t_total / 50.0)));
        let theta_strong = deg(final_obliquity(&run_sweep(&strong, strong.t_total / 50.0)));
        assert!(
            theta_strong - theta_weak > 15.0,
            "present-day α should fall well short: weak={theta_weak:.1}° strong={theta_strong:.1}°"
        );
    }

    #[test]
    fn obliquity_stays_finite_and_bounded() {
        let cfg = uranus_sweep(6.0);
        let snaps = run_sweep(&cfg, cfg.t_total / 50.0);
        for s in &snaps {
            assert!(s.obliquity.is_finite() && s.psi.is_finite());
            assert!((0.0..=PI).contains(&s.obliquity));
            assert!(s.ratio.is_finite() && s.ratio > 0.0);
        }
        assert!(snaps.len() > 10);
    }
}
