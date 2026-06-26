//! Reproduction of Hu, Huang, Gladman & Zhu (2025), "Early Stellar Flybys are
//! Unlikely: Improved Constraints from Sednoids and Large-q TNOs"
//! (arXiv:2505.16317).
//!
//! # The hypothesis under test
//!
//! A close stellar flyby is one proposed mechanism for *detaching* sednoids —
//! trans-Neptunian objects with large semi-major axes and perihelia lifted well
//! beyond Neptune's gravitational reach. The paper confronts this hypothesis
//! with two constraints that earlier flyby studies did not fold in together:
//!
//! 1. **Low inclinations.** The detached extreme TNOs have an inclination
//!    profile concentrated at `i < 30°`. A flyby that detaches them must not
//!    over-excite their inclinations past this envelope.
//! 2. **Primordial apsidal alignment.** The same population shows clustered
//!    longitudes of perihelion, which the flyby must not destroy.
//!
//! The authors' N-body experiments show that a flyby satisfies *both*
//! constraints only for special encounter geometries — a near-coplanar orbit
//! (`i_* ≈ 0°`) or an ecliptic-symmetric one (`ω_* ≈ 0°`, `i_* ≈ 90°`) — which
//! occupy a small fraction of the isotropic 4π of random encounter
//! orientations. Folding that small geometric acceptance together with the
//! occurrence rate of a sufficiently close encounter
//! (`q_* ≤ 1000 au`) in the Sun's birth cluster, the probability that such a
//! flyby actually happened **and** matched the data is `≲ 5%`.
//!
//! # What this crate reproduces
//!
//! An *analytic occurrence-rate* reconstruction of that `≲ 5%` headline — not
//! an N-body integration. The bound is built as a transparent product of
//! independent factors, each a documented estimate following the paper's logic:
//!
//! ```text
//!   P_total = P_encounter × f_geometry × f_success
//! ```
//!
//! - [`encounter`] — `P_encounter`: probability of at least one close
//!   (`q_* ≤ 1000 au`) flyby during cluster residence, from a focused
//!   encounter rate `Γ = n_* σ v` over a `~10 Myr` open-cluster lifetime.
//! - [`geometry`] — `f_geometry`: fraction of isotropic encounter orientations
//!   inside the coplanar-or-ecliptic-symmetric acceptance bands.
//! - `f_success` ([`F_SUCCESS`]) — conditional probability that an encounter of
//!   the right geometry actually produces a *sufficient* detached population
//!   without over-exciting it; a generous estimate from the paper's discussion.
//!
//! Each factor is recomputed from documented constants; the headline is
//! *computed*, never hard-coded. The goal is a reproducible bound, not a claim
//! of precision.

pub mod encounter;
pub mod geometry;

pub use encounter::{
    encounter_rate_per_day, expected_close_encounters, p_close_encounter, Q_STAR_MAX_AU,
};
pub use geometry::{coplanar_fraction, geometry_fraction, symmetric_fraction};

/// Published upper bound on the probability that a close-in flyby satisfying all
/// constraints actually happened (Hu et al. 2025). Reference constant only.
pub const PUBLISHED_MAX_PROBABILITY: f64 = 0.05;

/// Inclination envelope of the detached extreme TNOs (degrees). The flyby must
/// not pump inclinations beyond this. Reference constant only.
pub const INCLINATION_CONSTRAINT_DEG: f64 = 30.0;

/// Conditional success factor `f_success`: given a close encounter of the right
/// geometry, the probability it actually yields a sufficient, suitably aligned
/// sednoid population (rather than too few objects, or wrecking the alignment).
/// The paper finds this regime narrow but not vanishing; we take a generous
/// 0.5 so the final bound is not driven artificially low by this term.
pub const F_SUCCESS: f64 = 0.5;

/// Combined probability that a close-in (`q_* ≤ 1000 au`) flyby satisfying all
/// constraints actually occurred in the early Solar System:
///
/// `P_total = P_encounter × f_geometry × f_success`.
pub fn total_probability() -> f64 {
    p_close_encounter() * geometry_fraction() * F_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn total_probability_is_a_probability() {
        let p = total_probability();
        assert!((0.0..=1.0).contains(&p), "P_total out of range: {p}");
    }

    #[test]
    fn factors_multiply_to_total() {
        let expect = p_close_encounter() * geometry_fraction() * F_SUCCESS;
        assert!((total_probability() - expect).abs() < 1e-15);
    }
}
