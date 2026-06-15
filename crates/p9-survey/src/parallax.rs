//! Observer-baseline parallax helpers.
//!
//! The geometry lives in [`p9_core::coords::parallax`]; this module re-exports
//! it so survey code keeps a single import path.

pub use p9_core::constants::{AU_KM, EARTH_RADIUS_KM};
pub use p9_core::coords::parallax::{
    annual_parallax_arcsec, leo_baseline_km, leo_parallax_arcsec, leo_parallax_crossover_au,
    parallax_arcsec, ARCSEC_PER_RAD, LEO_ALTITUDE_KM,
};
