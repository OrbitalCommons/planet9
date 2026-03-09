//! Reproduction of Brown, Holman & Batygin (2024)
//! "A Pan-STARRS1 Search for Planet Nine"
//!
//! Extends the Planet Nine survey exclusion using the Pan-STARRS1 (PS1)
//! data set. PS1 adds depth (V~21.5) and sky coverage complementary to
//! ZTF and DES. The combined three-survey analysis excludes approximately
//! 78% of the prior parameter space for Planet Nine.

pub mod combined_exclusion;
pub mod detection_pipeline;
pub mod plots;
pub mod survey_model;
