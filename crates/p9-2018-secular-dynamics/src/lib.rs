//! Reproduction of Li, Hadden, Payne & Holman (2018), "The secular dynamics of
//! TNOs and Planet Nine" (arXiv:1806.06867).
//!
//! Headline: a semianalytic *secular* (orbit-averaged) theory of the ETNO–P9
//! interaction. A coplanar test ETNO in the combined field of the giant
//! planets (an axisymmetric J2 quadrupole that drives prograde apsidal
//! precession) and an eccentric Planet Nine (a Δϖ-dependent torque) acquires a
//! *forced eccentricity* and its longitude of perihelion ϖ either **librates**
//! (apsidal confinement / anti-alignment) or **circulates**.
//!
//! Which one happens depends on the test particle's semi-major axis: the
//! giants' free apsidal precession rate falls off steeply with a, while P9's
//! secular coupling rises with a (α = a/a₉ → 1). Above a **critical
//! semi-major axis** the P9 torque wins and Δϖ librates about the
//! anti-aligned equilibrium; below it the giants' precession dominates and the
//! apse circulates. Li et al. find pericenter libration is maintained for
//! a ≳ 250 AU (their abstract).
//!
//! What this crate computes, all from real (seeded) calculations reusing
//! `p9_core::analysis::secular` and `p9_core::forces::j2_secular` — never
//! re-implementing the secular machinery and with no hard-coded answers:
//!
//! - [`secular_model`]: the conserved coplanar secular Hamiltonian
//!   H(e, Δϖ; a) = P9's numerical Gauss-ring average + the giants' J2 term;
//!   the forced eccentricity e_forced(a) (the libration-centre eccentricity);
//!   the free apsidal precession rate; and the libration-vs-circulation
//!   classification by tracing the level curve through the (near-)circular
//!   initial condition. The critical semi-major axis a_crit is the bisection
//!   root of that transition.
//!
//! Documented model reductions (see `REPRODUCTION_NOTES.md`): this is the
//! coplanar (i = 0) reduction of Li et al.'s theory — it captures the apsidal
//! confinement / forced-eccentricity physics and the circulation→libration
//! transition, but not the inclination (Kozai/flip) sector that the full paper
//! also treats. The giants are the orbit-averaged J2 field of J+S+U
//! (`combined_j2_jsu`), not four bodies integrated directly.

pub mod secular_model;

/// Published reference values from Li, Hadden, Payne & Holman (2018), kept as
/// labelled constants for cross-checks — never asserted as if they were the
/// computed output.
pub mod published {
    /// Pericenter (Δϖ) libration is maintained above this approximate
    /// semi-major axis for the relevant perihelion range (paper abstract:
    /// "30 ≲ rₚ ≲ 80 AU & a ≳ 250 AU"). Used only as a scale check on the
    /// computed `a_crit`.
    pub const LIBRATION_A_THRESHOLD_AU: f64 = 250.0;

    /// Nominal Planet Nine mass considered by the paper (~10 M⊕).
    pub const P9_MASS_EARTH: f64 = 10.0;

    /// The paper favours an eccentric (e₉ ≳ 0.4) outer planet on a wide
    /// (a₉ ≈ 400–800 AU) orbit. We adopt e₉ = 0.6, a₉ = 500 AU as a nominal
    /// point inside that box.
    pub const P9_ECC_MIN: f64 = 0.4;
    pub const P9_A_MIN_AU: f64 = 400.0;
    pub const P9_A_MAX_AU: f64 = 800.0;
}
