//! Reproduction of Becker et al. (2018),
//! "Discovery and Dynamical Analysis of an Extreme Trans-Neptunian Object with
//! a High Orbital Inclination" — 2015 BP519 (arXiv:1805.05355).
//!
//! Headline: 2015 BP519 is a highly inclined (i ≈ 54°), extremely eccentric
//! (e ≈ 0.92), large-a (a ≈ 449 AU) *detached* ETNO. Such high-inclination
//! detached objects are naturally produced by a Planet Nine, which secularly
//! pumps the inclination of initially low-inclination scattered objects; and
//! BP519's longitude of perihelion is consistent with the P9-clustered
//! population.
//!
//! Modules:
//! * [`bp519`] — BP519's orbital elements with documented provenance
//!   (Becker 2018 discovery solution and the current JPL solution).
//! * [`pumping`] — secular inclination pumping of a scattered test particle by
//!   an inclined P9 (reusing `p9_core::analysis::secular`); an inclined P9
//!   raises a low-i particle to BP519-like inclinations while a coplanar
//!   control stays low.
//! * [`clustering`] — BP519's ϖ in the context of the Brown (2017) clustered
//!   sample (reusing `p9_core::analysis::circular`).

pub mod bp519;
pub mod clustering;
pub mod pumping;
