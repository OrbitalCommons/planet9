//! Reproduction of Sheppard & Trujillo (2016)
//! "New Extreme Trans-Neptunian Objects: Toward a Super-Earth in the Outer
//! Solar System" (arXiv:1608.08772).
//!
//! Headline claim reproduced here: the newly discovered extreme TNOs
//! (2014 SR349, 2013 FT28, 2014 SS349, 2013 UH15, 2013 FS28, 2014 FE72)
//! corroborate the orbital clustering attributed to a massive distant
//! perturber. The high-perihelion (q ≳ 40 AU) new objects cluster in argument
//! of perihelion ω near ~340°, and adding the new objects to the vetted
//! `p9_core::data::etno::BROWN_2017_SAMPLE` preserves significant
//! longitude-of-perihelion ϖ clustering. 2013 FT28 is the notable
//! anti-aligned object, sitting ~180° from the cluster in ϖ.
//!
//! All statistics are computed from the encoded elements with
//! `p9_core::analysis::circular`; the core ETNO table is reused verbatim, not
//! copied. Published reference values are kept as labelled constants in
//! [`discoveries::published`]. Residuals are documented in the module tests
//! and the workspace REPRODUCTION_NOTES.md.

pub mod clustering;
pub mod discoveries;

pub use clustering::{headline, ClusterStats, Headline};
pub use discoveries::SHEPPARD_2016_DISCOVERIES;
