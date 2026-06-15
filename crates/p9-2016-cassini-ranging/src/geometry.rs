//! Planet Nine heliocentric geometry along the Brown & Batygin orbit.
//!
//! The Cassini constraint is a function of WHERE P9 is on its orbit, i.e. of
//! its true anomaly `ν`. The generic element→position machinery now lives in
//! [`p9_core::types`]; this module just re-exports it under the paper-local
//! names the perturbation code uses.

pub use p9_core::types::{
    helio_distance_at_true_anomaly, position_at_true_anomaly as p9_position_at_true_anomaly,
    true_to_mean_anomaly, OrbitGeometry as P9Geometry,
};
