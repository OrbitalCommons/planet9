//! Colombo-top / Cassini-state geometry (Lu & Laughlin 2022, Eq. 8–9).
//!
//! In a frame co-precessing with the planet's orbit (whose pole precesses at
//! rate g < 0 about the invariable-plane normal, the two planes inclined by
//! I) the spin axis has up to four equilibria — the Cassini states. With
//! r ≡ α/|g| ≥ 0 they split across the two azimuths of the co-precessing
//! frame:
//!
//!   ψ = π (resonant family, state 2):  r cos θ sin θ − sin(θ − I) = 0,
//!   ψ = 0 (states 3, and 1/4 above the bifurcation):
//!                                      r cos θ sin θ + sin(θ − I) = 0.
//!
//! (This matches `resonance_capture`'s equation of motion
//! α cos θ sin θ + g sin(θ − I) with g < 0. A previous version used the
//! "+" form for state 2 — the θ → π − θ mirror of the wrong family — whose
//! middle root is the state-4 saddle: it DECREASED with r below the
//! bifurcation and jumped past 90°, contradicting the crate's own ≤ 90°
//! single-resonance asymptote.)
//!
//! States 1 and 4 appear when
//!
//!   r ≥ r_crit = (sin^{2/3} I + cos^{2/3} I)^{3/2}                    (Eq. 8)
//!
//! State 2 is the resonant branch relevant to adiabatic capture: it rises
//! from θ ≈ I at r ≪ 1 monotonically toward (but never past) 90° as r grows.
//! An outward-migrating Planet Nine slowly drives r upward (a primordially
//! fast spin precession, large α, is what the paper requires), dragging
//! Uranus' spin to high obliquity along this branch.

use std::f64::consts::PI;

/// Residual of the RESONANT (ψ = π) Cassini condition, hosting state 2:
///
///   f(θ) = r · cos θ · sin θ − sin(θ − i)
///
/// The single root on θ ∈ (0, π) is the state-2 obliquity.
pub fn cassini_residual(theta: f64, r: f64, i: f64) -> f64 {
    r * theta.cos() * theta.sin() - (theta - i).sin()
}

/// Residual of the ψ = 0 Cassini condition, hosting state 3 (always) and
/// states 1 and 4 (above the bifurcation r ≥ r_crit):
///
///   f(θ) = r · cos θ · sin θ + sin(θ − i)
pub fn cassini_residual_psi0(theta: f64, r: f64, i: f64) -> f64 {
    r * theta.cos() * theta.sin() + (theta - i).sin()
}

/// Critical r = α/|g| at which the second pair of Cassini states (1 and 4)
/// appears (Eq. 8): r_crit = (sin^{2/3} I + cos^{2/3} I)^{3/2}.
pub fn critical_ratio(i: f64) -> f64 {
    (i.sin().abs().powf(2.0 / 3.0) + i.cos().abs().powf(2.0 / 3.0)).powf(1.5)
}

/// Roots of a residual on θ ∈ [0, π] by sign-change bracketing + bisection.
fn roots_on_0_pi(f: impl Fn(f64) -> f64) -> Vec<f64> {
    let n = 4000;
    let mut roots = Vec::new();
    let mut prev_theta = 0.0;
    let mut prev_f = f(prev_theta);
    for k in 1..=n {
        let theta = PI * (k as f64) / (n as f64);
        let fv = f(theta);
        if prev_f == 0.0 {
            roots.push(prev_theta);
        } else if prev_f * fv < 0.0 {
            let (mut lo, mut hi) = (prev_theta, theta);
            let mut flo = prev_f;
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let fm = f(mid);
                if flo * fm <= 0.0 {
                    hi = mid;
                } else {
                    lo = mid;
                    flo = fm;
                }
            }
            roots.push(0.5 * (lo + hi));
        }
        prev_theta = theta;
        prev_f = fv;
    }
    roots
}

/// Solve for all Cassini-state obliquities at ratio `r = α/|g| ≥ 0` and
/// inclination `i`: the union of the ψ = π root (state 2) and the ψ = 0
/// roots (state 3, plus states 1 and 4 above the bifurcation). Returns
/// obliquities (rad) in ascending order — 2 states below r_crit, 4 above.
pub fn cassini_states(r: f64, i: f64) -> Vec<f64> {
    let mut roots = roots_on_0_pi(|t| cassini_residual(t, r, i));
    roots.extend(roots_on_0_pi(|t| cassini_residual_psi0(t, r, i)));
    roots.sort_by(|a, b| a.partial_cmp(b).unwrap());
    roots
}

/// Obliquity (rad) of Cassini state 2 — the resonant (ψ = π) branch followed
/// during adiabatic capture. Rises monotonically from θ ≈ I at r ≪ 1 toward
/// (never past) 90° as r → ∞; e.g. at I = 20°: 46.7° at r = 0.9, ~71° at
/// r = 3, → 90⁻.
pub fn cassini_state_2(r: f64, i: f64) -> f64 {
    let roots = roots_on_0_pi(|t| cassini_residual(t, r, i));
    *roots
        .first()
        .expect("the resonant psi=pi branch always has its state-2 root")
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::constants::DEG2RAD;

    #[test]
    fn critical_ratio_matches_paper_form() {
        // At I = 20°, (α/g)_crit ≈ 1.74 (Eq. 8 evaluated).
        let i = 20.0 * DEG2RAD;
        let c = critical_ratio(i);
        assert_relative_eq!(c, 1.7432, epsilon = 1e-3);
        // Bounds: between 1 (I→0 or 90°) and 2 (I=45°).
        assert!(c >= 1.0 && c <= 2.0 + 1e-9);
    }

    #[test]
    fn critical_ratio_peaks_near_45_degrees() {
        // r_crit = (sin^{2/3}I + cos^{2/3}I)^{3/2} peaks at I=45° with value 2.
        let mid = critical_ratio(45.0 * DEG2RAD);
        let lo = critical_ratio(5.0 * DEG2RAD);
        let hi = critical_ratio(85.0 * DEG2RAD);
        assert!(mid > lo && mid > hi, "crit should peak near 45°");
        assert_relative_eq!(mid, 2.0, epsilon = 1e-9);
    }

    #[test]
    fn two_states_below_critical_four_above() {
        // The bifurcation at r_crit: states 2 and 3 below it, states 1 and 4
        // appearing above (Eq. 8).
        let i = 20.0 * DEG2RAD;
        let crit = critical_ratio(i);
        assert_eq!(cassini_states(0.5 * crit, i).len(), 2);
        assert_eq!(cassini_states(2.0 * crit, i).len(), 4);
    }

    #[test]
    fn slow_precession_pins_obliquity_near_inclination() {
        // r = α/|g| ≪ 1 (slow spin precession): the single Cassini state sits
        // close to the orbit-normal inclination — low obliquity.
        let i = 20.0 * DEG2RAD;
        let theta2 = cassini_state_2(0.2, i);
        assert!(
            (theta2 - i).abs() < 10.0 * DEG2RAD,
            "θ₂ = {:.1}° should be near I=20° for slow precession",
            theta2 / DEG2RAD
        );
    }

    #[test]
    fn state_2_rises_monotonically_toward_90_with_r() {
        // The resonant state-2 obliquity rises monotonically with r from
        // θ ≈ I toward — but never past — 90°: the single-resonance asymptote
        // (REPRODUCTION_NOTES §19). The old mirrored-family solver DECREASED
        // below the bifurcation and jumped past 90°.
        let i = 20.0 * DEG2RAD;
        let mut prev = 0.0;
        for &r in &[0.3, 0.9, 1.5, 3.0, 10.0, 100.0] {
            let t = cassini_state_2(r, i);
            assert!(
                t > prev && t < 90.0 * DEG2RAD,
                "θ₂(r={r}) = {:.2}°, prev {:.2}°",
                t / DEG2RAD,
                prev / DEG2RAD
            );
            prev = t;
        }
        // Pinned reference value: θ₂(r = 0.9, I = 20°) ≈ 46.7°.
        let t09 = cassini_state_2(0.9, i);
        assert!(
            (t09 - 46.7 * DEG2RAD).abs() < 1.0 * DEG2RAD,
            "θ₂(0.9) = {:.2}°",
            t09 / DEG2RAD
        );
    }

    #[test]
    fn state_2_matches_resonance_capture_equilibrium() {
        // Cross-check against the EOM module's ψ = π condition
        // α·cosθ·sinθ + g·sin(θ−I) with g < 0: the state-2 root must satisfy
        // it to solver precision.
        let i = 20.0 * DEG2RAD;
        let (alpha, g): (f64, f64) = (2.0, -1.0);
        let t = cassini_state_2(alpha / g.abs(), i);
        let residual = alpha * t.cos() * t.sin() + g * (t - i).sin();
        assert!(residual.abs() < 1e-9, "EOM residual = {residual:.2e}");
    }

    #[test]
    fn residual_is_zero_at_returned_states() {
        // Every returned state is a root of one of the two azimuth branches.
        let i = 20.0 * DEG2RAD;
        for &r in &[50.0, 2.0, 1.0, 0.3] {
            for theta in cassini_states(r, i) {
                let f_pi = cassini_residual(theta, r, i).abs();
                let f_0 = cassini_residual_psi0(theta, r, i).abs();
                assert!(
                    f_pi.min(f_0) < 1e-6,
                    "θ = {theta:.4} not a root of either branch for r={r}"
                );
            }
        }
    }
}
