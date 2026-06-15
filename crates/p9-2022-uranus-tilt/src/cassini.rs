//! Colombo-top / Cassini-state geometry (Lu & Laughlin 2022, Eq. 8–9).
//!
//! In a frame co-precessing with the planet's orbit (whose pole precesses at
//! rate g about the invariable-plane normal, the two planes inclined by I) the
//! spin axis has up to four equilibria — the Cassini states. Writing the ratio
//! r ≡ α/|g| (positive; the sign of the nodal regression is absorbed into the
//! sin θ cos θ term as the Ward–Hamilton convention), their obliquities θ_c
//! measured from the orbit normal solve
//!
//!   r cos θ_c sin θ_c + sin(θ_c − I) = 0.                            (Eq. 9)
//!
//! There are four Cassini states when
//!
//!   r ≥ r_crit = (sin^{2/3} I + cos^{2/3} I)^{3/2}                    (Eq. 8)
//!
//! and only state 2 (and its mirror, state 3) otherwise. State 2 is the
//! resonant branch relevant to adiabatic capture: it migrates from θ_c ≈ I at
//! r ≪ 1 up through 90° and beyond as r grows. An outward-migrating Planet
//! Nine slowly drives r upward (a primordially fast spin precession, large α,
//! is what the paper requires), dragging Uranus' spin to high obliquity.

use std::f64::consts::PI;

/// Right-hand side of the Cassini equilibrium relation (Eq. 9), as a function
/// of obliquity θ for ratio `r = α/|g| ≥ 0` and orbital inclination `i`.
///
///   f(θ) = r · cos θ · sin θ + sin(θ − i)
///
/// Cassini states are the roots f(θ) = 0 on θ ∈ [0, π].
pub fn cassini_residual(theta: f64, r: f64, i: f64) -> f64 {
    r * theta.cos() * theta.sin() + (theta - i).sin()
}

/// Critical r = α/|g| at which the second pair of Cassini states (1 and 4)
/// appears (Eq. 8): r_crit = (sin^{2/3} I + cos^{2/3} I)^{3/2}.
pub fn critical_ratio(i: f64) -> f64 {
    (i.sin().abs().powf(2.0 / 3.0) + i.cos().abs().powf(2.0 / 3.0)).powf(1.5)
}

/// Solve for all Cassini-state obliquities at ratio `r = α/|g| ≥ 0` and
/// inclination `i`. Roots are found by sign-change bracketing on a fine grid
/// followed by bisection. Returns obliquities (rad) in ascending order.
pub fn cassini_states(r: f64, i: f64) -> Vec<f64> {
    let n = 4000;
    let mut roots = Vec::new();
    let mut prev_theta = 0.0;
    let mut prev_f = cassini_residual(prev_theta, r, i);
    for k in 1..=n {
        let theta = PI * (k as f64) / (n as f64);
        let f = cassini_residual(theta, r, i);
        if prev_f == 0.0 {
            roots.push(prev_theta);
        } else if prev_f * f < 0.0 {
            // Bisect on [prev_theta, theta].
            let (mut lo, mut hi) = (prev_theta, theta);
            let (mut flo, _fhi) = (prev_f, f);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                let fm = cassini_residual(mid, r, i);
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
        prev_f = f;
    }
    roots
}

/// Obliquity (rad) of Cassini state 2 — the resonant branch followed during
/// adiabatic capture. `r = α/|g| ≥ 0`. Below r_crit there is a single root
/// (state 2), which sits near θ ≈ I for small r and climbs toward 90° as r
/// approaches r_crit. Above r_crit three roots exist (states 1, 2, 4 by
/// ascending θ); state 2 is the middle root, which continues past 90° toward
/// the high-obliquity values the paper attributes to Uranus.
pub fn cassini_state_2(r: f64, i: f64) -> f64 {
    let roots = cassini_states(r, i);
    match roots.len() {
        0 => panic!("at least one Cassini state must exist"),
        1 => roots[0],
        _ => {
            // The middle of the (typically three) roots is state 2.
            roots[roots.len() / 2]
        }
    }
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
    fn one_state_below_critical_three_above() {
        // The bifurcation at r_crit: a single Cassini state below it, two new
        // ones (states 1 and 4) appearing above (Eq. 8).
        let i = 20.0 * DEG2RAD;
        let crit = critical_ratio(i);
        assert_eq!(cassini_states(0.5 * crit, i).len(), 1);
        assert_eq!(cassini_states(2.0 * crit, i).len(), 3);
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
    fn near_bifurcation_drives_high_obliquity() {
        // Just above r_crit the resonant Cassini state 2 sits past 90° — the
        // paper's high-obliquity branch that captures Uranus' spin.
        let i = 20.0 * DEG2RAD;
        let crit = critical_ratio(i);
        let theta2 = cassini_state_2(1.05 * crit, i);
        assert!(
            theta2 > 90.0 * DEG2RAD,
            "θ₂ = {:.1}° should exceed 90° just past the bifurcation",
            theta2 / DEG2RAD
        );
    }

    #[test]
    fn state_2_grows_monotonically_as_r_falls_to_crit() {
        // Above the bifurcation, the resonant state-2 obliquity rises
        // monotonically as r decreases toward r_crit (from ~90° at large r
        // toward ~125° near the bifurcation): a clean monotone dependence on
        // the controlling ratio.
        let i = 20.0 * DEG2RAD;
        let crit = critical_ratio(i);
        let ratios = [50.0, 10.0, 4.0, 2.0, 1.1 * crit];
        let mut prev = 0.0;
        for &r in &ratios {
            let t = cassini_state_2(r, i);
            assert!(
                t >= prev - 1e-6,
                "θ₂ not monotone: r={r}, θ₂={:.2}° prev={:.2}°",
                t / DEG2RAD,
                prev / DEG2RAD
            );
            prev = t;
        }
    }

    #[test]
    fn residual_is_zero_at_returned_states() {
        let i = 20.0 * DEG2RAD;
        for &r in &[50.0, 2.0, 1.0, 0.3] {
            for theta in cassini_states(r, i) {
                assert!(
                    cassini_residual(theta, r, i).abs() < 1e-6,
                    "f({theta:.4}) not a root for r={r}"
                );
            }
        }
    }
}
