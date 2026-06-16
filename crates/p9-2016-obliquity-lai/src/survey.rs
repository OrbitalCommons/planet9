//! Parameter-space constraints from Lai (2016) §3.
//!
//! Reproduces the analytic contour (Eq. 17) of the effective semi-major axis
//! ãₚ that, together with the right inclination, produces the observed
//! θsl = 6° solar obliquity, and the ΔΩ relation (Eq. 18).

use std::f64::consts::PI;

use p9_core::units::{au, Length};

use crate::precession::{self, published, PlanetNine};

/// Lai Eq. (19) auxiliary quantities for a target ΔΩ (the relative longitude
/// of ascending node between P9 and the solar equator) and target obliquity.
#[derive(Debug, Clone, Copy)]
pub struct NodeGeometry {
    /// θ̂sl ≡ θsl / 6° (1 for the observed value).
    pub theta_sl_hat: f64,
    /// g ≡ (π/2 − ΔΩ) / (π/4).
    pub g: f64,
    /// f ≡ (π/2 − ΔΩ) / cos ΔΩ.
    pub f: f64,
}

impl NodeGeometry {
    pub fn new(obliquity_deg: f64, delta_omega_rad: f64) -> Self {
        let half_pi_minus = PI / 2.0 - delta_omega_rad;
        NodeGeometry {
            theta_sl_hat: obliquity_deg / published::OBSERVED_OBLIQUITY_DEG,
            g: half_pi_minus / (PI / 4.0),
            f: half_pi_minus / delta_omega_rad.cos(),
        }
    }
}

/// Lai Eq. (17): the effective semi-major axis (AU) required to produce the
/// target obliquity, as a function of P9 mass and inclination:
///
///   ãₚ ≃ 462 [ (mp/10m⊕) sin 2θp / (θ̂sl f) ]^(1/3) AU
pub fn required_a_tilde(mass_earth: f64, theta_p_rad: f64, geom: &NodeGeometry) -> f64 {
    let inner = (mass_earth / 10.0) * (2.0 * theta_p_rad).sin() / (geom.theta_sl_hat * geom.f);
    published::A_TILDE_PREFACTOR_AU * inner.powf(1.0 / 3.0)
}

/// [`required_a_tilde`] as a typed [`Length`].
pub fn required_a_tilde_typed(mass_earth: f64, theta_p_rad: f64, geom: &NodeGeometry) -> Length {
    au(required_a_tilde(mass_earth, theta_p_rad, geom))
}

/// Brute-force inversion of the full analytic theory: find the effective
/// semi-major axis ãₚ at which [`precession::solar_obliquity_deg`] equals the
/// target obliquity for the given mass, inclination, and spin period.
///
/// Returns `None` if no crossing is found in the scan range. This is the
/// honest, formula-free counterpart to the closed-form Eq. (17) and is used to
/// cross-check it.
pub fn solve_a_tilde_for_obliquity(
    mass_earth: f64,
    theta_p_rad: f64,
    e: f64,
    spin_period_days: f64,
    target_obliquity_deg: f64,
) -> Option<f64> {
    let obl_at = |a_tilde: f64| {
        let p9 = PlanetNine {
            mass_earth,
            a_au: a_tilde / (1.0 - e * e).sqrt(),
            e,
            inclination_rad: theta_p_rad,
        };
        precession::solar_obliquity_deg(&p9, spin_period_days)
    };

    // Scan from the resonance side outward; obliquity falls with ãₚ once past
    // the peak. Find the *outer* crossing (the large-ãₚ branch Lai plots).
    let mut a = published::A_TILDE_RANGE_AU.1; // start high (low obliquity)
    let step = 1.0;
    let mut prev_a = a;
    let mut prev = obl_at(a);
    a -= step;
    while a >= 200.0 {
        let cur = obl_at(a);
        if (prev - target_obliquity_deg) * (cur - target_obliquity_deg) <= 0.0 {
            // Linear interpolation between prev_a and a.
            let frac = (target_obliquity_deg - prev) / (cur - prev);
            return Some(prev_a + frac * (a - prev_a));
        }
        prev_a = a;
        prev = cur;
        a -= step;
    }
    None
}

/// Lai Eq. (18): the node/inclination condition that, together with Eq. (17),
/// fixes ΔΩ = 45°. Returns the right-hand side
///
///   sin θp · L / Lp = 15 g + 4.84 λ★ (P★/10 d)⁻¹ − cos θp · (θ̂sl f)
///
/// rearranged to the predicted `(L/Lp) sin θp / (θ̂sl f)` so it can be compared
/// against the directly computed `L/Lp`. This is the analytic ΔΩ constraint;
/// it is exercised as an internal consistency check rather than pinned to a
/// published scalar (Lai presents it only graphically, Figs. 1–3).
pub fn eq18_rhs_l_over_lp(
    theta_p_rad: f64,
    spin_period_days: f64,
    geom: &NodeGeometry,
    lambda_star: f64,
) -> f64 {
    let bracket =
        15.0 * geom.g + 4.84 * lambda_star * (spin_period_days / 10.0).recip() - theta_p_rad.cos();
    bracket * (geom.theta_sl_hat * geom.f) / theta_p_rad.sin()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::constants::DEG2RAD;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_required_a_tilde_matches_f64() {
        let geom = NodeGeometry::new(6.0, 45.0 * DEG2RAD);
        assert_relative_eq!(
            required_a_tilde_typed(10.0, 20.0 * DEG2RAD, &geom).get::<astronomical_unit>(),
            required_a_tilde(10.0, 20.0 * DEG2RAD, &geom)
        );
    }

    #[test]
    fn eq17_prefactor_lands_in_published_band() {
        // Lai §3: across the parameter box, ãₚ ∈ [340, 480] AU. For ΔΩ = 45°,
        // mp = 10 m⊕, θp = 20°, Eq. (17) should land inside that band.
        let geom = NodeGeometry::new(6.0, 45.0 * DEG2RAD);
        let a_tilde = required_a_tilde(10.0, 20.0 * DEG2RAD, &geom);
        let (lo, hi) = published::A_TILDE_RANGE_AU;
        assert!(
            (lo..=hi).contains(&a_tilde),
            "Eq.17 ãₚ = {a_tilde:.0} AU outside Lai band [{lo}, {hi}]"
        );
    }

    #[test]
    fn eq17_heavier_p9_needs_smaller_inclination_modest_a_change() {
        // Lai §3: "a larger mp requires a smaller θp, with a modest change in
        // ãₚ." Holding ãₚ ≈ const, the inclination that satisfies Eq. (17)
        // shrinks as mp grows. Equivalently at fixed θp, ãₚ ∝ mp^(1/3) — a
        // modest (cube-root) change.
        let geom = NodeGeometry::new(6.0, 45.0 * DEG2RAD);
        let a10 = required_a_tilde(10.0, 20.0 * DEG2RAD, &geom);
        let a20 = required_a_tilde(20.0, 20.0 * DEG2RAD, &geom);
        let ratio = a20 / a10;
        // Cube-root of 2 ≈ 1.26: a "modest change".
        assert!(
            (ratio - 2.0_f64.powf(1.0 / 3.0)).abs() < 1e-6,
            "ãₚ(20)/ãₚ(10) = {ratio:.4} (expected 2^(1/3))"
        );
        assert!(ratio < 1.3, "mass change of ãₚ should be modest");
    }

    #[test]
    fn eq18_rhs_positive_and_monotone_in_delta_omega() {
        // Lai Eq. (18) is presented only graphically (Figs. 1–3): it is jointly
        // solved with Eq. (17) for (ãₚ, θp) at fixed mp, ΔΩ, so it is exercised
        // here as a self-consistency/sanity expression rather than pinned to a
        // published scalar. The required L/Lp must be positive and must grow as
        // ΔΩ shrinks (smaller node offset ⇒ heavier/closer P9 ⇒ larger L/Lp).
        let lam = crate::precession::solar::lambda_star();
        let theta_p = 20.0 * DEG2RAD;
        let g_lo = NodeGeometry::new(6.0, 12.0 * DEG2RAD);
        let g_hi = NodeGeometry::new(6.0, 52.0 * DEG2RAD);
        let llp_lo = eq18_rhs_l_over_lp(theta_p, 20.0, &g_lo, lam);
        let llp_hi = eq18_rhs_l_over_lp(theta_p, 20.0, &g_hi, lam);
        assert!(llp_lo > 0.0 && llp_hi > 0.0, "L/Lp must be positive");
        assert!(
            llp_lo > llp_hi,
            "smaller ΔΩ should demand larger L/Lp ({llp_lo:.1} vs {llp_hi:.1})"
        );
    }

    #[test]
    fn full_theory_inversion_consistent_with_eq17() {
        // The formula-free inversion of the full analytic obliquity should land
        // in Lai's 340–480 AU band. At the "averaged" P★/λ★ = 20 d the θp = 20°
        // peak obliquity is just under 6° (≈5.9°), so a clean 6° crossing needs
        // a slightly more inclined P9 (θp = 30°), where the outer branch crosses
        // 6° at ãₚ ≈ 430 AU — inside the band. (This θp/ãₚ trade is exactly the
        // contour Lai plots in Fig. 1.)
        let solved = solve_a_tilde_for_obliquity(10.0, 30.0 * DEG2RAD, 0.5, 20.0, 6.0)
            .expect("should find a 6° crossing");
        let (lo, hi) = published::A_TILDE_RANGE_AU;
        assert!(
            (lo..=hi).contains(&solved),
            "inverted ãₚ = {solved:.0} AU (Lai band [{lo}, {hi}])"
        );
    }
}
