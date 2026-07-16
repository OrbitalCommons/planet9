//! Orbital-pole geometry: unit pole vectors from (i, Ω) and directional
//! statistics over sets of poles.
//!
//! The pole of an orbit with inclination `i` and longitude of ascending node
//! `Ω` (ecliptic frame) is
//!
//! ```text
//!   p̂ = (sin i sin Ω, −sin i cos Ω, cos i)
//! ```
//!
//! Clustering of poles (e.g. the ETNO mean-plane arguments) is measured with
//! the vector resultant, the 3D analogue of the circular mean-resultant
//! length in [`super::circular`].

use crate::constants::RAD2DEG;

/// Orbital pole unit vector (ecliptic frame) from inclination `i` and
/// longitude of ascending node `omega_big`, both in radians.
pub fn pole_vector(i: f64, omega_big: f64) -> [f64; 3] {
    let (si, ci) = (i.sin(), i.cos());
    let (so, co) = (omega_big.sin(), omega_big.cos());
    [si * so, -si * co, ci]
}

/// Mean resultant length of a set of pole unit vectors: |Σ p̂| / n, in [0, 1].
/// 1 = perfectly aligned poles, → 0 for isotropically scattered poles.
pub fn pole_resultant_length(poles: &[[f64; 3]]) -> f64 {
    if poles.is_empty() {
        return 0.0;
    }
    let s = vector_sum(poles);
    (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt() / poles.len() as f64
}

/// Normalized mean pole direction of a set of pole unit vectors. Returns the
/// +z axis if the resultant is degenerate (near-zero).
pub fn mean_pole_direction(poles: &[[f64; 3]]) -> [f64; 3] {
    let s = vector_sum(poles);
    let mag = (s[0] * s[0] + s[1] * s[1] + s[2] * s[2]).sqrt();
    if mag < 1e-12 {
        return [0.0, 0.0, 1.0];
    }
    [s[0] / mag, s[1] / mag, s[2] / mag]
}

/// Angular separation between two pole directions, in degrees.
pub fn pole_separation_deg(p1: [f64; 3], p2: [f64; 3]) -> f64 {
    let dot = p1[0] * p2[0] + p1[1] * p2[1] + p1[2] * p2[2];
    dot.clamp(-1.0, 1.0).acos() * RAD2DEG
}

fn vector_sum(poles: &[[f64; 3]]) -> [f64; 3] {
    let mut s = [0.0_f64; 3];
    for p in poles {
        s[0] += p[0];
        s[1] += p[1];
        s[2] += p[2];
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{FRAC_PI_2, PI};

    #[test]
    fn pole_of_zero_inclination_is_z() {
        let p = pole_vector(0.0, 1.2345);
        assert!((p[0]).abs() < 1e-15 && (p[1]).abs() < 1e-15);
        assert!((p[2] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn pole_is_unit_and_node_sets_azimuth() {
        for &(i, om) in &[(0.3, 0.0), (1.0, 2.0), (FRAC_PI_2, PI)] {
            let p = pole_vector(i, om);
            let n2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
            assert!((n2 - 1.0).abs() < 1e-14);
            assert!((p[2] - i.cos()).abs() < 1e-14);
        }
    }

    #[test]
    fn resultant_and_mean_direction_behave() {
        let aligned = vec![pole_vector(0.1, 0.0); 5];
        assert!((pole_resultant_length(&aligned) - 1.0).abs() < 1e-14);
        // Two opposite poles cancel.
        let opposite = vec![[0.0, 0.0, 1.0], [0.0, 0.0, -1.0]];
        assert!(pole_resultant_length(&opposite) < 1e-14);
        // Mean of two poles symmetric about z is z.
        let sym = vec![pole_vector(0.4, 0.0), pole_vector(0.4, PI)];
        let m = mean_pole_direction(&sym);
        assert!(pole_separation_deg(m, [0.0, 0.0, 1.0]) < 1e-6);
    }

    #[test]
    fn separation_is_a_metric_on_the_sphere() {
        let a = pole_vector(0.2, 0.3);
        assert!(pole_separation_deg(a, a) < 1e-9);
        let b = pole_vector(0.2 + 0.5, 0.3);
        assert!((pole_separation_deg(a, b) - 0.5 * RAD2DEG).abs() < 1e-9);
        assert!((pole_separation_deg([0.0, 0.0, 1.0], [0.0, 0.0, -1.0]) - 180.0).abs() < 1e-12);
    }
}
