//! Reproduction of Payne, Holman & Pál (2019),
//! "A TESS Search for Distant Solar System Objects: Yield Estimates"
//! (arXiv:1911.03676).
//!
//! TESS full-frame images, digitally shifted-and-stacked along trial orbits,
//! reach a 50%-completeness depth of I_C ≈ 22.0 ± 0.5 per sector — deep
//! enough that "TESS could discover the hypothesized Planet Nine." But P9's
//! reflected-light brightness (V ≈ 19-25 over the plausible mass/distance
//! box) means only a *bright* P9 — close and/or massive — is within reach.
//!
//! This crate computes three real quantities (no hard-coded answers):
//!
//! 1. The reflected-light apparent magnitude of a P9 of given mass/distance,
//!    via [`p9_core::analysis::photometry`] (Neptune-anchored mass-radius,
//!    opposition geometry). See [`yield_grid::p9_apparent_magnitude`].
//! 2. The achievable stacked TESS limiting magnitude, from a sky-noise
//!    co-addition scaling law reproducing the published I_C ≈ 22.0.
//!    See [`tess::TessStack`].
//! 3. The *detectable fraction* of a reference (mass, distance) box: the mean
//!    TESS detection efficiency over the grid. See [`yield_grid::ReferenceBox`].
//!
//! The detectable fraction lies in (0, 1), increases with assumed mass, and
//! decreases with distance — the central message of the paper.

pub mod tess;
pub mod yield_grid;

use rand::Rng;
use rand_distr::{Distribution, Uniform};

use crate::tess::TessStack;
use crate::yield_grid::{p9_apparent_magnitude, ReferenceBox};

/// Monte-Carlo estimate of the detectable fraction of a reference box, drawing
/// `n` (mass, distance) samples uniformly from the box with a seeded RNG. A
/// cross-check on the deterministic grid average in
/// [`ReferenceBox::detectable_fraction`]; both estimate the same integral.
pub fn detectable_fraction_mc<R: Rng>(
    bx: &ReferenceBox,
    stack: &TessStack,
    n: usize,
    rng: &mut R,
) -> f64 {
    if n == 0 {
        return 0.0;
    }
    let mass_dist = Uniform::new_inclusive(bx.mass_earth.0, bx.mass_earth.1);
    let dist_dist = Uniform::new_inclusive(bx.distance_au.0, bx.distance_au.1);
    let mut sum = 0.0;
    for _ in 0..n {
        let mass = mass_dist.sample(rng);
        let dist = dist_dist.sample(rng);
        let v = p9_apparent_magnitude(mass, dist, bx.albedo);
        sum += stack.detection_efficiency(v);
    }
    sum / n as f64
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;

    #[test]
    fn monte_carlo_matches_grid_average() {
        let bx = ReferenceBox::default();
        let stack = TessStack::default();
        let grid = bx.detectable_fraction(&stack);
        let mut rng = rand::rngs::StdRng::seed_from_u64(2019);
        let mc = detectable_fraction_mc(&bx, &stack, 200_000, &mut rng);
        assert!((mc - grid).abs() < 0.02, "MC {mc:.4} vs grid {grid:.4}");
        assert!(mc > 0.0 && mc < 1.0);
    }
}
