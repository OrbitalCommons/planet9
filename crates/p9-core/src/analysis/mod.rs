pub mod elements;
pub mod secular;

/// Snapshot of orbital elements for all active particles at a given time.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct Snapshot {
    /// Time in days from epoch
    pub t: f64,
    /// Semi-major axes (AU)
    pub a: Vec<f64>,
    /// Eccentricities
    pub e: Vec<f64>,
    /// Inclinations (radians)
    pub i: Vec<f64>,
    /// Arguments of perihelion (radians)
    pub omega: Vec<f64>,
    /// Longitudes of ascending node (radians)
    pub omega_big: Vec<f64>,
    /// Mean anomalies (radians)
    pub mean_anomaly: Vec<f64>,
    /// Particle IDs
    pub ids: Vec<u64>,
}
