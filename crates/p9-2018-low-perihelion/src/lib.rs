//! Reproduction of Cáceres & Gomes (2018)
//! "The evolution of TNOs under a low-perihelion Planet Nine"
//! (arXiv:1808.01248).
//!
//! HEADLINE. The canonical Brown & Batygin (2016) Planet Nine sits on a
//! detached orbit with a large perihelion, q9 = a9(1−e9) ≈ 280 AU (a9 = 700,
//! e9 = 0.6). Cáceres & Gomes argue that a Planet Nine with a *lower*
//! perihelion — a larger eccentricity at the same semi-major axis — can STILL
//! produce and maintain the observed apsidal confinement of the extreme TNOs
//! (a ≳ 250 AU, q ≳ 40 AU). Lowering q9 therefore *widens* the dynamically
//! allowed (a9, e9) region rather than ruling P9 out.
//!
//! WHAT THIS CRATE COMPUTES (no hard-coded answers). Using the doubly-averaged
//! Gauss-ring secular Hamiltonian of `p9_core::analysis::secular`
//! (`numerical_secular_hamiltonian`), we treat an ETNO as a coplanar inner test
//! particle perturbed by Planet Nine and measure the strength of apsidal
//! confinement as a function of the P9 perihelion q9:
//!
//! - [`confinement`] builds the secular phase portrait H(e, Δϖ) for a test
//!   ETNO at fixed a, locates the favoured apsidal equilibrium (the energy
//!   minimum over Δϖ ∈ {0, π}), and measures the half-width of its Δϖ-libration
//!   island plus the dimensionless apsidal well depth — the confinement metric.
//!   It also reports the forced eccentricity and the torque amplitude. (The
//!   single-ring test particle favours the *aligned* apse; the OBSERVED
//!   anti-alignment needs P9's precession and the giants' field — see the
//!   `confinement` module's sign caveat.)
//! - [`sweep`] sweeps e9 at fixed a9, mapping each (a9, e9) to its perihelion
//!   q9 = a9(1−e9), and tabulates how the confinement metric, forced
//!   eccentricity and torque vary as q9 decreases. The headline result is that
//!   the metric stays non-zero (libration persists) down to a documented
//!   low-q9 minimum.
//! - [`reference`] holds the canonical B&B q9 and the explored low-q9 range as
//!   labelled constants, plus the vetted ETNO sample reused from p9-core.
//!
//! The secular forcing rate is cross-checked against the p9-core numerical ring
//! average directly (same kernel) and against the analytic coplanar quadrupole
//! in the small-α limit; residuals are documented in the tests.

pub mod confinement;
pub mod reference;
pub mod sweep;
