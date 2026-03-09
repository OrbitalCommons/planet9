//! Reproduction of Batygin & Brown (2016)
//! "Evidence for a Distant Giant Planet in the Solar System"
//!
//! This crate reproduces the five main computational results from the paper:
//! 1. KBO stability analysis (Section 2)
//! 2. Analytical secular phase-space portraits (Section 3)
//! 3. N-body phase-space portraits (Section 4)
//! 4. Synthetic scattered disk — planar (Section 5.1)
//! 5. Synthetic scattered disk — 3D (Section 5.2)

pub mod kbo_elements;
pub mod octupole;
pub mod phase_portrait;
pub mod plots;
pub mod resonance;
pub mod scattered_disk_sim;
