//! Reproduction of Lai (2016)
//! "Solar Obliquity Induced by Planet Nine: Analytic Theory" (arXiv:1608.01421).
//!
//! Where the sibling crate [`p9_2016_obliquity`] integrates the full vector
//! secular Hamiltonian numerically, Lai (2016) writes down a *closed-form*
//! analytic theory for the same physics and this crate reproduces it.
//!
//! The picture is three coupled torques:
//!   1. Planet Nine forces the canonical-planet plane `l̂` to precess about the
//!      total angular momentum `ĵ` at the rate `Ω_L cos θp` (Eqs. 2–5).
//!   2. The Sun's spin axis `Ŝ★` precesses about `l̂` at the rate `Ω★ps`,
//!      driven by the planets' torque on the solar rotational quadrupole
//!      (J₂) (Eqs. 7–8).
//!   3. The induced obliquity `θsl` follows from solving the spin equation in
//!      the frame corotating with `l̂` (Eqs. 9–15). It is largest near the
//!      secular spin–orbit resonance `Ω★ps ≈ Ω_L cos θp (L/Lp + cos θp)`,
//!      where the precession-frequency mismatch `Ωz → 0`.
//!
//! All headline quantities are computed from the giant + terrestrial planet
//! masses and the solar J₂; nothing is hard-coded. The observed 6° obliquity
//! and Lai's published coefficients (87.5 Gyr, 55.8 Gyr, 462 AU prefactor)
//! are kept as labelled constants and pinned against the computation in tests.
//!
//! The planet angular momenta reuse the sibling crate's
//! [`p9_2016_obliquity::solar_model`] giant-planet table and angular-momentum
//! helper; the solar J₂ and masses come from `p9_core::constants`.

pub mod planets;
pub mod precession;
pub mod survey;
