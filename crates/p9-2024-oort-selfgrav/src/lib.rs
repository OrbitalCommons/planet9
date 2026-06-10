//! Reproduction of Batygin & Nesvorny (2024)
//! "Self-Gravitational Dynamics Within the Inner Oort Cloud"
//!
//! Examines how the cumulative self-gravity of the Inner Oort Cloud (IOC),
//! which may approach several Earth masses, drives secular oscillations akin
//! to the von Zeipel-Lidov-Kozai (vZLK) resonance. The IOC potential is
//! modeled via the Miyamoto-Nagai axisymmetric disk profile.
//!
//! Key result: vZLK cycles can drive perihelion down to ~100 AU, but the
//! evolutionary timescales (tens of Gyr) vastly exceed the solar system age
//! of 4.5 Gyr, making IOC self-gravity negligible for current architecture.
//!
//! IOC mass: up to 3 Earth masses (from Nesvorny et al. cluster2 model).

pub mod conservation;
pub mod hamiltonian;
pub mod miyamoto_nagai;
pub mod plots;
pub mod vzlk;
