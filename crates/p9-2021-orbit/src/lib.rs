//! Reproduction of Brown & Batygin (2021)
//! "The Orbit of Planet Nine"
//!
//! Uses an MCMC analysis of 11 distant KBOs (a > 250 AU) to derive
//! posterior distributions for Planet Nine's orbital parameters.
//! Key results: m = 6.2 (+2.2/-1.3) Earth masses, a = 380 (+140/-80) AU,
//! i = 16 +/- 5 degrees, q = 300 (+85/-60) AU.

pub mod plots;
pub mod posterior;
pub mod reference_population;
pub mod statistical_measures;
