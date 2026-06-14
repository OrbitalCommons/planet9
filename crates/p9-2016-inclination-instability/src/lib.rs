//! Reproduction of Madigan & McCourt (2016),
//! "The inclination instability in the Outer Solar System"
//! (arXiv:1509.08920), with the Das & Batygin (2023) suppression study
//! (arXiv:2307.00378).
//!
//! # The instability (Madigan & McCourt 2016)
//!
//! A massive disc of eccentric planetesimals (several to tens of Earth masses)
//! orbiting the Sun is subject to a collective *secular* instability: the
//! orbital inclinations rise exponentially while the arguments of perihelion
//! and the longitudes of the ascending node cluster. The instability is driven
//! entirely by the disc's own self-gravity — **no planet is required**. Because
//! it is secular (orbit-averaged), the relevant clock is the disc's secular
//! precession time, which is the Keplerian orbital time scaled by the
//! disc-to-Sun mass ratio:
//!
//!   t_sec ~ (M_sun / M_d) * P_orbit,    omega_sec ~ (M_d / M_sun) * n,
//!
//! with n = sqrt(GM_sun / a^3) the mean motion at the disc's characteristic
//! semi-major axis a (Madigan & McCourt 2016, Eq. 1 and Section 2; the
//! self-gravity secular frequency of a near-Keplerian disc). The instability
//! grows on this secular time, so the e-folding growth rate scales as
//!
//!   gamma = C_growth * (M_d / M_sun) * n,
//!
//! a *linear-theory* growth rate proportional to the disc mass (collective
//! self-gravity). C_growth is an order-unity dimensionless eigenvalue of the
//! linearised secular problem; Madigan & McCourt's N-body experiments give
//! e-folding times of order a few disc secular times, and we adopt the labelled
//! constant [`growth::C_GROWTH`] documented there.
//!
//! See [`growth`] for the growth rate / e-folding time, and [`suppression`]
//! for the Das & Batygin (2023) critical-mass criterion.

pub mod growth;
pub mod suppression;
