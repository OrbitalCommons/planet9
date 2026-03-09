//! Reproduction of Brown & Batygin (2021)
//! "A Search for Planet Nine Using the Zwicky Transient Facility Public Archive"
//!
//! Searches the ZTF public data release for moving objects consistent with
//! Planet Nine. Synthetic P9 orbits are injected into the ZTF cadence to
//! calibrate detection efficiency. No Planet Nine candidate is found,
//! excluding approximately 56% of the prior parameter space.

pub mod detection_efficiency;
pub mod exclusion;
pub mod plots;
pub mod survey_model;
