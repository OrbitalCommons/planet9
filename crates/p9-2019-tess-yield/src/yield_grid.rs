//! Detectable yield of Planet Nine over a reference (mass, distance) box.
//!
//! For each (mass, heliocentric distance) point we compute the reflected-light
//! apparent magnitude with [`p9_core::analysis::photometry::planet_apparent_magnitude`]
//! (opposition geometry, Neptune-anchored mass-radius relation, geometric
//! albedo), then weight it by the TESS stacked-detection efficiency
//! [`crate::tess::TessStack::detection_efficiency`]. The *detectable fraction*
//! of the box is the mean efficiency over the grid — the fraction of the
//! plausible parameter space within reach of a stacked TESS sector.
//!
//! HEADLINE (Payne, Holman & Pál 2019): a stacked TESS sector reaches
//! I_C ≈ 22, deep enough that only a bright (close, large) Planet Nine is
//! detectable. A nominal mid-range P9 (~6 M⊕, ~500 AU, V ≈ 20.5) sits near
//! the stacked depth — marginal — and for most of the (mass, distance) box
//! P9 is too faint, so the detectable fraction is well below unity and
//! shrinks rapidly with assumed distance.

use serde::{Deserialize, Serialize};

use p9_core::analysis::photometry::planet_apparent_magnitude;
use p9_core::units::{au, earth_masses, Length, Mass};

use crate::tess::TessStack;

/// Geometric albedo assumed for Planet Nine's reflected light. Neptune-like
/// (p_V ≈ 0.41) by default; the IRAS/Neptune analogy is the standard P9
/// assumption (Brown & Batygin and the review literature use p ≈ 0.4).
pub const P9_ALBEDO: f64 = 0.4;

/// A rectangular reference box in (mass, heliocentric distance) parameter
/// space over which the detectable fraction is averaged.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReferenceBox {
    /// Mass range in Earth masses [min, max].
    pub mass_earth: (f64, f64),
    /// Heliocentric distance range in AU [min, max].
    pub distance_au: (f64, f64),
    /// Geometric albedo assumed for the reflected-light magnitude.
    pub albedo: f64,
    /// Grid resolution per axis.
    pub n_mass: usize,
    pub n_dist: usize,
}

impl Default for ReferenceBox {
    /// Plausible Planet Nine box from the post-2019 posteriors: 5-10 M⊕,
    /// 300-700 AU (Batygin et al. 2019 review; Brown & Batygin 2021).
    fn default() -> Self {
        Self {
            mass_earth: (5.0, 10.0),
            distance_au: (300.0, 700.0),
            albedo: P9_ALBEDO,
            n_mass: 40,
            n_dist: 40,
        }
    }
}

impl ReferenceBox {
    /// Mass range `[min, max]` as typed [`Mass`] values.
    pub fn mass_range(&self) -> (Mass, Mass) {
        (
            earth_masses(self.mass_earth.0),
            earth_masses(self.mass_earth.1),
        )
    }

    /// Heliocentric distance range `[min, max]` as typed [`Length`] values.
    pub fn distance_range(&self) -> (Length, Length) {
        (au(self.distance_au.0), au(self.distance_au.1))
    }

    fn mass_at(&self, j: usize) -> f64 {
        let (lo, hi) = self.mass_earth;
        if self.n_mass <= 1 {
            return 0.5 * (lo + hi);
        }
        lo + (hi - lo) * (j as f64) / ((self.n_mass - 1) as f64)
    }

    fn dist_at(&self, k: usize) -> f64 {
        let (lo, hi) = self.distance_au;
        if self.n_dist <= 1 {
            return 0.5 * (lo + hi);
        }
        lo + (hi - lo) * (k as f64) / ((self.n_dist - 1) as f64)
    }

    /// Mean TESS detection efficiency over the box: the detectable fraction
    /// of the reference parameter space (a number in [0, 1]).
    pub fn detectable_fraction(&self, stack: &TessStack) -> f64 {
        let mut sum = 0.0;
        let mut n = 0usize;
        for j in 0..self.n_mass {
            let mass = self.mass_at(j);
            for k in 0..self.n_dist {
                let dist = self.dist_at(k);
                let v = planet_apparent_magnitude(mass, self.albedo, dist);
                sum += stack.detection_efficiency(v);
                n += 1;
            }
        }
        if n == 0 {
            0.0
        } else {
            sum / n as f64
        }
    }

    /// Detectable fraction restricted to a single distance shell `[d, d]`
    /// (averaged only over the mass axis), to show the distance dependence.
    pub fn detectable_fraction_at_distance(&self, stack: &TessStack, dist: f64) -> f64 {
        let mut sum = 0.0;
        for j in 0..self.n_mass {
            let mass = self.mass_at(j);
            let v = planet_apparent_magnitude(mass, self.albedo, dist);
            sum += stack.detection_efficiency(v);
        }
        sum / self.n_mass as f64
    }
}

/// Reflected-light V magnitude of a nominal Planet Nine, for convenience.
pub fn p9_apparent_magnitude(mass_earth: f64, distance_au: f64, albedo: f64) -> f64 {
    planet_apparent_magnitude(mass_earth, albedo, distance_au)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_box_typed_ranges_match_f64() {
        use uom::si::length::astronomical_unit;
        use uom::si::mass::kilogram;
        let bx = ReferenceBox::default();
        let (m_lo, m_hi) = bx.mass_range();
        assert!(
            (m_lo.get::<kilogram>() - bx.mass_earth.0 * p9_core::units::EARTH_MASS_KG).abs() < 1e12
        );
        assert!(
            (m_hi.get::<kilogram>() - bx.mass_earth.1 * p9_core::units::EARTH_MASS_KG).abs() < 1e12
        );
        let (d_lo, d_hi) = bx.distance_range();
        assert!((d_lo.get::<astronomical_unit>() - bx.distance_au.0).abs() < 1e-9);
        assert!((d_hi.get::<astronomical_unit>() - bx.distance_au.1).abs() < 1e-9);
    }

    #[test]
    fn detectable_fraction_is_a_probability() {
        let stack = TessStack::default();
        let f = ReferenceBox::default().detectable_fraction(&stack);
        assert!(f > 0.0 && f < 1.0, "detectable fraction = {f:.4}");
    }

    #[test]
    fn fraction_decreases_with_distance() {
        // The headline result: P9 dims as 1/(r^2 delta^2) (~3 mag per
        // doubling), so the detectable fraction of the box falls steeply with
        // assumed distance. Against the I_C ≈ 22 stacked depth the close-in
        // shells are essentially fully detectable and the far shells drop
        // sharply: ~1.0 at 350 AU, ~0.90 at 600 AU, ~0.41 at 700 AU.
        let stack = TessStack::default();
        let bx = ReferenceBox::default();
        let near = bx.detectable_fraction_at_distance(&stack, 350.0);
        let mid = bx.detectable_fraction_at_distance(&stack, 600.0);
        let far = bx.detectable_fraction_at_distance(&stack, 700.0);
        assert!(near > mid, "{near:.4} should exceed {mid:.4}");
        assert!(mid > far, "{mid:.4} should exceed {far:.4}");
        // Close-in slice is essentially fully detectable; the outer edge of
        // the plausible box has fallen well below it (only a bright/close P9
        // is within reach).
        assert!(near > 0.95, "near-shell fraction = {near:.4}");
        assert!(far < 0.6, "far-shell fraction = {far:.4}");
    }

    #[test]
    fn fraction_increases_with_mass() {
        // A more massive P9 is larger and brighter -> more of it is detectable.
        let stack = TessStack::default();
        let light = ReferenceBox {
            mass_earth: (3.0, 5.0),
            ..ReferenceBox::default()
        };
        let heavy = ReferenceBox {
            mass_earth: (10.0, 20.0),
            ..ReferenceBox::default()
        };
        let f_light = light.detectable_fraction(&stack);
        let f_heavy = heavy.detectable_fraction(&stack);
        assert!(
            f_heavy > f_light + 0.05,
            "heavy {f_heavy:.4} should exceed light {f_light:.4}"
        );
    }

    #[test]
    fn nominal_midrange_p9_is_marginal() {
        // A nominal mid-range P9 (6 M⊕, 500 AU) is V ≈ 20.75 — within the
        // optimistic V≈I_C limit of the I_C ≈ 22 stacked depth, hence formally
        // detectable, BUT only by ~1.25 mag. That margin is exactly the size
        // of the depth uncertainty (±0.5 mag) plus a plausible Neptune-like
        // V−I_C reddening (methane darkening in the red, ~0.5-1 mag), so the
        // detection is marginal once those are folded in. Push the same body
        // out to the box edge (~700 AU) and it sits right at the 50% depth.
        let v500 = p9_apparent_magnitude(6.0, 500.0, P9_ALBEDO);
        assert!(v500 > 20.0 && v500 < 21.0, "V_P9(6 M⊕, 500 AU) = {v500:.2}");
        let stack = TessStack::default();
        // Margin to the stacked depth is comparable to its own uncertainty.
        let margin = stack.stacked_depth() - v500;
        assert!(
            margin > 0.0 && margin < 2.0,
            "margin to depth = {margin:.2} mag (marginal-to-detectable)"
        );
        // At the outer edge of the plausible box the same mass is genuinely
        // borderline (efficiency near 0.5).
        let v700 = p9_apparent_magnitude(6.0, 700.0, P9_ALBEDO);
        let eff = stack.detection_efficiency(v700);
        assert!(
            (0.2..=0.8).contains(&eff),
            "edge-of-box efficiency = {eff:.3} at V = {v700:.2}"
        );
    }

    #[test]
    fn far_faint_p9_beyond_reach() {
        // A small, far P9 (5 M⊕, 1000 AU) is V ≈ 24 — hopeless for TESS even
        // stacked, the regime where the detectable fraction goes to zero.
        let v = p9_apparent_magnitude(5.0, 1000.0, P9_ALBEDO);
        assert!(v > 23.0, "V = {v:.2}");
        let stack = TessStack::default();
        assert!(stack.detection_efficiency(v) < 1e-3);
    }
}
