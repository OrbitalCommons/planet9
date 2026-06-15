//! Reproduction of Beust (2016),
//! "Orbital clustering of distant Kuiper belt objects by Planet 9 —
//! secular resonance" (arXiv:1605.02473).
//!
//! # Headline mechanism
//!
//! The apsidal clustering of the extreme trans-Neptunian objects (ETNOs)
//! arises because a distant test particle can be **captured in a secular
//! resonance** with Planet Nine: the particle's free apsidal precession rate
//! matches Planet Nine's, so the relative longitude of perihelion
//!
//! ```text
//!     ψ = Δϖ = ϖ_particle − ϖ_9
//! ```
//!
//! does not circulate through all values but instead **librates** — it is
//! confined to oscillate about a fixed center. Beust's phase-space portraits
//! (orbit-averaged secular dynamics, computed with no expansion in the
//! semi-major-axis ratio because α = a/a₉ is not small) show *apsidally
//! anti-aligned, high-eccentricity libration islands*: closed level curves of
//! the secular Hamiltonian that encircle the anti-aligned equilibrium at
//! Δϖ = π. A particle started inside such an island stays anti-aligned with
//! Planet Nine; one started outside circulates and shows no preferred ϖ.
//!
//! # What this crate computes (real, from p9-core)
//!
//! The secular dynamics is one degree of freedom at fixed semi-major axis
//! (`a` is a secular invariant). Reusing
//! [`p9_core::analysis::secular::numerical_secular_hamiltonian`] — the exact
//! doubly-averaged 1/Δ interaction by Gauss-ring quadrature, the same engine
//! Batygin & Brown (2016) used because the multipole series diverges for the
//! relevant α — we build the conserved Hamiltonian `H(e, Δϖ)` and:
//!
//! * locate the anti-aligned secular equilibrium (the libration center) and
//!   the eccentricity band where the libration island exists ([`portrait`]);
//! * integrate the 1-DOF secular equations of motion in the canonical pair
//!   `(Δϖ, action)` and classify a trajectory as libration or circulation
//!   ([`libration`]);
//! * show the resonance/island strengthens with Planet Nine's mass;
//! * cross-check the Hamiltonian against the p9-core ring average.
//!
//! No published number is hard-coded as an answer — the libration center and
//! amplitude are read off the computed portrait. The paper's qualitative
//! claims are kept as labelled reference constants in [`published`].

pub mod hamiltonian;
pub mod libration;
pub mod portrait;
pub mod published;

pub use hamiltonian::SecularModel;
pub use libration::{classify, Trajectory, TrajectoryKind};
pub use portrait::{libration_island, Island, ResonantHamiltonian};
