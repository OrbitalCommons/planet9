//! Synthetic scattered disk simulations.
//!
//! Reproduces the planar (Section 5.1) and 3D (Section 5.2) scattered disk
//! experiments from Batygin & Brown (2016).
//!
//! The simulation orchestrator now lives in
//! [`p9_core::initial_conditions::scattered_disk_sim`] (it wraps the
//! particle generators in `p9_core::initial_conditions::scattered_disk`).
//! It is re-exported here for this crate's existing call sites.

pub use p9_core::initial_conditions::scattered_disk_sim::{
    clustering_statistics, j2_jsun_force, run_scattered_disk, DiskSimConfig, DiskSnapshot,
};
