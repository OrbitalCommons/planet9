//! Reproduction of Phan, Goto, Yamamura et al. (2025)
//! "A Search for Planet Nine with IRAS and AKARI Data"
//!
//! Searches for Planet Nine candidates by cross-matching far-infrared
//! sources from the IRAS Faint Source Catalogue (~1983) and the AKARI
//! Monthly Unconfirmed Source List (~2006), separated by ~23 years.
//!
//! At 500–700 AU, a 7–17 Earth-mass Planet Nine would drift ~3'/yr,
//! producing angular separations of 42'–69.6' between the two survey
//! epochs. The paper finds 13 candidate pairs after selection cuts,
//! with one good candidate surviving image inspection.
//!
//! Key parameters reproduced here:
//! - Mass range: 7–17 M_Earth
//! - Distance range: 500–700 AU
//! - Effective temperature: ~30–50 K (internal heat dominated)
//! - IRAS 60 µm sensitivity: ~0.2 Jy
//! - AKARI 90 µm sensitivity: ~0.55 Jy
//! - Epoch separation: ~23 years
//! - Circular-orbit transverse rate: ~1.2'–1.9'/yr (h/r²; eccentric higher)
//! - Angular separation window: 42'–69.6' (heliocentric motion at e ≤ 0.7
//!   plus geocentric parallax, ±6.9' per epoch at 500 AU)

pub mod candidate_search;
pub mod orbital_constraints;
pub mod plots;
pub mod proper_motion;
pub mod survey_model;
pub mod thermal_model;
