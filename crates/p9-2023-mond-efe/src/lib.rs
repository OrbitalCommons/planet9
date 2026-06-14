//! Reproduction of Brown & Mathur (2023),
//! "Modified Newtonian Dynamics as an Alternative to the Planet Nine
//! Hypothesis" (arXiv:2304.00576; companion Migaszewski 2023, arXiv:2303.13339).
//!
//! In MOND, the solar system is embedded in the external galactic field, and
//! the resulting External Field Effect (EFE) imposes an anomalous *quadrupole*
//! tidal field on the outer solar system whose symmetry axis points toward the
//! GALACTIC CENTER. This field secularly torques ETNO orbits and drives a
//! *forced apsidal alignment*: the longitudes of perihelion are pushed toward
//! the galactic-center direction — an alternative explanation for the observed
//! ETNO clustering with no Planet Nine.
//!
//! This crate computes the EFE quadrupole field, its orbit-averaged secular
//! disturbing function, the forced apsidal equilibrium ϖ_forced, the libration
//! about it, and the linear scaling of the precession rate with the EFE
//! amplitude — then pins the headline claim (ϖ_forced ≈ galactic-center
//! ecliptic longitude) in seeded tests.
//!
//! Reuses `p9_core::coords::sky` for the galactic-center direction (SPICE
//! `GALACTIC` frame) and `p9_core::constants` for `GM_sun`.

pub mod efe;
pub mod secular;
