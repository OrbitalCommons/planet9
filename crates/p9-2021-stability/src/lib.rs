//! Reproduction of Batygin, Mardling & Nesvorny (2021)
//! "The Stability Boundary of the Distant Scattered Disk"
//!
//! Derives an analytical model for the dynamics of distant TNOs showing the
//! scattered disk architecture is shaped by an infinite chain of 2:j resonances
//! with Neptune. Resonance overlap drives chaotic motion when perihelion drops
//! below a critical value.
//!
//! Key result: q_crit = a_N * sqrt(ln((24^2/5)(m_N/M_sun)(a/a_N)^(5/2)))
//!
//! Diffusion coefficient D_a ~ (8/5pi)(m_N/M_sun)*sqrt(GM_sun*a_N)*exp(-(q/a_N)^2/2)
//! is independent of semi-major axis.

pub mod chirikov;
pub mod diffusion;
pub mod hansen;
pub mod plots;
pub mod resonance_chain;
pub mod stability;
