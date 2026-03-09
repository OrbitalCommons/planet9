//! Orbital element computation and snapshot recording.

use super::Snapshot;
use crate::constants::GM_SUN;
use crate::types::{cartesian_to_elements, StateVector};

/// Compute orbital elements for all active particles and record a snapshot.
pub fn record_snapshot(
    particles: &[StateVector],
    active: &[bool],
    ids: &[u64],
    t: f64,
) -> Snapshot {
    let mut snap = Snapshot {
        t,
        a: Vec::new(),
        e: Vec::new(),
        i: Vec::new(),
        omega: Vec::new(),
        omega_big: Vec::new(),
        mean_anomaly: Vec::new(),
        ids: Vec::new(),
    };

    for (idx, particle) in particles.iter().enumerate() {
        if !active[idx] {
            continue;
        }

        let elem = cartesian_to_elements(particle, GM_SUN);
        snap.a.push(elem.a);
        snap.e.push(elem.e);
        snap.i.push(elem.i);
        snap.omega.push(elem.omega);
        snap.omega_big.push(elem.omega_big);
        snap.mean_anomaly.push(elem.mean_anomaly);
        snap.ids.push(ids[idx]);
    }

    snap
}

/// Compute longitude of perihelion (varpi = Omega + omega) for all elements
/// in a snapshot.
pub fn longitudes_of_perihelion(snap: &Snapshot) -> Vec<f64> {
    snap.omega_big
        .iter()
        .zip(snap.omega.iter())
        .map(|(&big_o, &o)| (big_o + o) % crate::constants::TWO_PI)
        .collect()
}

/// Compute the difference in longitude of perihelion relative to Planet Nine.
/// Delta_varpi = varpi_particle - varpi_p9 (mod 2*pi, centered on [-pi, pi]).
pub fn delta_varpi(snap: &Snapshot, varpi_p9: f64) -> Vec<f64> {
    let varpis = longitudes_of_perihelion(snap);
    varpis
        .iter()
        .map(|&v| {
            let mut dv = v - varpi_p9;
            while dv > std::f64::consts::PI {
                dv -= crate::constants::TWO_PI;
            }
            while dv < -std::f64::consts::PI {
                dv += crate::constants::TWO_PI;
            }
            dv
        })
        .collect()
}
