//! Reproduction of Siraj, Loeb & Siegel (2025), "Planet Y" (arXiv:2508.14156).
//!
//! Distinct from Planet Nine: a Mercury-to-Earth-mass object (≈ 0.05–1 M_⊕) at
//! a ≈ 100–200 AU, modestly inclined, imposes a **warp** of the Kuiper belt's
//! mean plane — the local Laplace plane / forced inclination tilts up over a
//! localized band of semi-major axis, a ≈ 50–80 AU, where the distant
//! perturber's nodal torque becomes comparable to the giant planets'.
//!
//! - [`laplace_plane`] computes the forced inclination i_forced(a) from the
//!   competition between the giant-planet quadrupole (inner torque, reusing
//!   `p9_core::forces::j2_secular`) and Planet Y's exterior ring torque
//!   (∝ m_Y), and measures the warp amplitude.
//! - [`mean_plane`] expresses the forced plane as an observable mean-plane pole
//!   (reusing p9-core pole geometry). (A previous "debiasing" here was an
//!   x·s/s identity around a free efficiency-slope knob — deleted rather
//!   than dressed up; a real ε(i) would have to be derived over the
//!   librating population.)

pub mod laplace_plane;
pub mod mean_plane;
