//! Modified Chirikov overlap criterion with asymmetric resonance widths.
//!
//! The weaker resonance interacts with the stronger only up to the full-width
//! extent of the weaker resonance. Uses Farey sequence for spatial density
//! between resonances.

use serde::{Deserialize, Serialize};

use crate::resonance::Resonance;

/// Check if two adjacent resonances overlap.
///
/// Modified criterion: weaker resonance only interacts up to its own full width.
pub fn resonances_overlap(r1: &Resonance, r2: &Resonance) -> bool {
    let separation = (r1.a_nominal - r2.a_nominal).abs();
    let min_width = r1.delta_a.min(r2.delta_a);
    // Overlap when separation < sum of half-widths, but limited by weaker
    separation < r1.delta_a + min_width
}

/// Overlap parameter for adjacent resonances.
///
/// K = (delta_a1 + delta_a2) / |a1 - a2|
pub fn overlap_parameter(r1: &Resonance, r2: &Resonance) -> f64 {
    let separation = (r1.a_nominal - r2.a_nominal).abs();
    if separation < 1e-10 {
        return f64::INFINITY;
    }
    (r1.delta_a + r2.delta_a) / separation
}

/// Count the number of overlapping pairs in a resonance chain.
pub fn count_overlaps(chain: &[Resonance]) -> usize {
    let mut count = 0;
    for i in 1..chain.len() {
        if resonances_overlap(&chain[i - 1], &chain[i]) {
            count += 1;
        }
    }
    count
}

/// MEGNO chaos indicator classification.
///
/// log2(Y) < 1.01: regular
/// 1.01 - 3.0: intermediate
/// > 4.0: chaotic
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum MegnoClass {
    Regular,
    Intermediate,
    Chaotic,
}

pub fn classify_megno(log2_y: f64) -> MegnoClass {
    if log2_y < 1.01 {
        MegnoClass::Regular
    } else if log2_y < 3.0 {
        MegnoClass::Intermediate
    } else {
        MegnoClass::Chaotic
    }
}

/// Farey sequence mediant for finding resonance density.
///
/// Given fractions p1/q1 and p2/q2, the mediant is (p1+p2)/(q1+q2).
pub fn farey_mediant(p1: i64, q1: i64, p2: i64, q2: i64) -> (i64, i64) {
    (p1 + p2, q1 + q2)
}

/// Compute the number of resonances between two given m:j resonances
/// using the Farey sequence density.
pub fn farey_density(j1: i64, j2: i64, max_order: i64) -> usize {
    if (j2 - j1).abs() <= 1 {
        return 0;
    }

    let mut count = 0;
    for j in (j1 + 1)..j2 {
        if j <= max_order {
            count += 1;
        }
    }
    count
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_res(a: f64, da: f64) -> Resonance {
        Resonance {
            m: 2,
            j: 1,
            l: 2,
            a_nominal: a,
            delta_a: da,
        }
    }

    #[test]
    fn test_overlap_touching() {
        let r1 = make_res(100.0, 5.0);
        let r2 = make_res(108.0, 5.0);
        assert!(resonances_overlap(&r1, &r2));
    }

    #[test]
    fn test_no_overlap_far() {
        let r1 = make_res(100.0, 1.0);
        let r2 = make_res(200.0, 1.0);
        assert!(!resonances_overlap(&r1, &r2));
    }

    #[test]
    fn test_overlap_parameter() {
        let r1 = make_res(100.0, 5.0);
        let r2 = make_res(110.0, 5.0);
        let k = overlap_parameter(&r1, &r2);
        assert!((k - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_classify_megno_regular() {
        assert_eq!(classify_megno(0.5), MegnoClass::Regular);
    }

    #[test]
    fn test_classify_megno_chaotic() {
        assert_eq!(classify_megno(5.0), MegnoClass::Chaotic);
    }

    #[test]
    fn test_farey_mediant() {
        let (p, q) = farey_mediant(1, 2, 1, 3);
        assert_eq!((p, q), (2, 5));
    }
}
