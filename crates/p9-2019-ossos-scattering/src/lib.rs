//! Reproduction of Kaib et al. (2019), "OSSOS XV: Probing the Distant Solar
//! System with Observed Scattering TNOs" (arXiv:1905.09286).
//!
//! Headline: the observed actively-scattering TNO population — large
//! semi-major axis objects whose perihelia lie near Neptune (q ≈ 30–38 AU) and
//! which undergo a semimajor-axis random walk — constrains Planet Nine. A
//! massive distant P9 lifts perihelia at large a (a secular effect) and so
//! sculpts the scattering → detached transition; the properties of the
//! observed scattering population (its a-distribution and the perihelion at
//! which objects detach from Neptune) DISFAVOR some P9 configurations.
//!
//! What this crate computes, all from real (seeded) calculations reusing
//! `p9_core`, with no hard-coded answers:
//!
//! - `boundary`: the scattering/detached divide in the (a, q) plane from the
//!   Chirikov resonance-overlap criterion in `p9_core::analysis::resonance`
//!   (overlap parameter K and the critical perihelion q_crit(a), the K = 1
//!   root). Maps the boundary for a ≈ 50–1000 AU.
//! - `planet_nine`: P9's secular perihelion lifting, from the numerical
//!   Gauss-ring secular Hamiltonian in `p9_core::analysis::secular`.
//! - `population`: a seeded synthetic scattering population and the detached
//!   fraction with vs without P9 — the perihelion lift pushes a computed,
//!   non-zero fraction across the divide.
//!
//! Documented model reductions (see `REPRODUCTION_NOTES.md`): the boundary is
//! the analytic 2:j overlap criterion, not a fitted N-body diffusion map; the
//! P9 lift is a coplanar test-particle secular cycle, capturing apsidal-driven
//! perihelion raising but not Neptune's simultaneous a random walk or Kozai
//! inclination coupling.

pub mod boundary;
pub mod planet_nine;
pub mod population;

/// Published reference values from Kaib et al. (2019) and the P9 literature,
/// kept as labelled constants for cross-checks (never asserted as if computed).
pub mod published {
    /// Scattering/detached boundary scale: the observed scattering population
    /// detaches with perihelia near the Neptune region, q ≈ 37 AU at the
    /// relevant large semi-major axes (Kaib et al. 2019; the scattering disk
    /// is conventionally defined out to q ≈ 37–38 AU). Used to sanity-check the
    /// computed q_crit(a) scale.
    pub const SCATTERING_Q_BOUNDARY_AU: f64 = 37.0;

    /// Qualitative constraint: a sufficiently massive, distant P9 over-detaches
    /// the large-a scattering population relative to what OSSOS observes, which
    /// is the basis for disfavoring some P9 configurations. There is no single
    /// published scalar for this; recorded as a flag for documentation.
    pub const P9_DISFAVORS_OVER_DETACHMENT: bool = true;
}
