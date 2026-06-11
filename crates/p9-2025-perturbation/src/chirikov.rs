//! Modified Chirikov overlap criterion with asymmetric resonance widths,
//! Farey enumeration of the intermediate-resonance "fine comb", and a
//! MEGNO chaos indicator computed from trajectory-divergence series.

use serde::{Deserialize, Serialize};

use crate::resonance::Resonance;

/// Check if two adjacent resonances overlap.
///
/// Modified Chirikov criterion (Belyakov & Batygin 2025): the *narrower*
/// resonance interacts with the wider one only up to its own full-width
/// extent (2× its half-width), capped at the wider one's half-width. The
/// implementation is symmetric in its arguments (the previous form used
/// `r1`'s half-width unconditionally, so `overlap(r1, r2) ≠
/// overlap(r2, r1)` and results depended on chain iteration direction).
pub fn resonances_overlap(r1: &Resonance, r2: &Resonance) -> bool {
    let separation = (r1.a_nominal - r2.a_nominal).abs();
    let wide = r1.delta_a.max(r2.delta_a);
    let narrow = r1.delta_a.min(r2.delta_a);
    separation < wide + (2.0 * narrow).min(wide)
}

/// Overlap parameter for adjacent resonances (symmetric):
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
/// log2(Y) < 1.01: regular (⟨Y⟩ → 2 for quasi-periodic orbits)
/// 1.01 - 3.0: intermediate
/// > 3.0: chaotic (⟨Y⟩ grows linearly with the Lyapunov exponent)
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

/// Mean MEGNO ⟨Y⟩ from a tangent/shadow-trajectory separation series
/// δ(t_k) sampled at uniform spacing `dt` (two-particle approximation,
/// Cincotta & Simó 2000):
///
///   Y(t) = (2/t) ∫₀ᵗ (δ̇/δ) s ds,   ⟨Y⟩ = (1/t) ∫₀ᵗ Y ds
///
/// evaluated with midpoint discretization of d(ln δ). For regular orbits
/// (power-law divergence) ⟨Y⟩ → 2; for chaotic orbits ⟨Y⟩ ≈ λt/2 grows
/// without bound. Feed `classify_megno(megno.log2())`.
pub fn megno_mean_from_separation(deltas: &[f64], dt: f64) -> f64 {
    assert!(deltas.len() >= 3, "need at least 3 separation samples");
    assert!(dt > 0.0);

    let n = deltas.len();
    let mut y_sum = 0.0; // running ∫ (δ̇/δ) s ds
    let mut mean_acc = 0.0;
    let mut count = 0usize;

    for k in 1..n {
        let dln = (deltas[k] / deltas[k - 1]).ln();
        let t_mid = (k as f64 - 0.5) * dt;
        y_sum += dln * t_mid; // (δ̇/δ)·s·ds with midpoint s
        let t_k = k as f64 * dt;
        let y_k = 2.0 * y_sum / t_k;
        mean_acc += y_k;
        count += 1;
    }
    mean_acc / count as f64
}

/// Farey mediant: given fractions p1/q1 and p2/q2, the mediant is
/// (p1+p2)/(q1+q2) — the unique lowest-denominator fraction strictly
/// between two Farey neighbours.
pub fn farey_mediant(p1: i64, q1: i64, p2: i64, q2: i64) -> (i64, i64) {
    (p1 + p2, q1 + q2)
}

fn gcd(mut a: i64, mut b: i64) -> i64 {
    while b != 0 {
        (a, b) = (b, a % b);
    }
    a.abs()
}

/// Enumerate every reduced fraction m/j with m ≤ `max_m` strictly between
/// the frequency ratios m1/j1 and m2/j2 — the intermediate m:j resonances
/// that form the paper's "fine comb" between two low-order resonances.
/// (n/n_N = m/j, so each fraction is an m:j mean-motion resonance with
/// Neptune.) Returned sorted by frequency ratio (ascending), i.e. by
/// descending semi-major axis.
///
/// This is the Farey sequence of order-(max_m numerator) restricted to
/// the open interval — a complete enumeration, replacing the previous
/// stub that just counted integers.
pub fn enumerate_intermediate_resonances(
    m1: i64,
    j1: i64,
    m2: i64,
    j2: i64,
    max_m: i64,
) -> Vec<(i64, i64)> {
    let f_lo = (m1 as f64 / j1 as f64).min(m2 as f64 / j2 as f64);
    let f_hi = (m1 as f64 / j1 as f64).max(m2 as f64 / j2 as f64);
    assert!(f_lo > 0.0, "frequency ratios must be positive");

    let mut out = Vec::new();
    for m in 1..=max_m {
        // m/j ∈ (f_lo, f_hi)  ⇔  m/f_hi < j < m/f_lo
        let j_min = (m as f64 / f_hi).floor() as i64 + 1;
        let j_max = ((m as f64 / f_lo).ceil() as i64 - 1).max(0);
        for j in j_min..=j_max {
            if j <= 0 {
                continue;
            }
            let f = m as f64 / j as f64;
            if f <= f_lo || f >= f_hi {
                continue; // endpoint or boundary fraction
            }
            if gcd(m, j) == 1 {
                out.push((m, j));
            }
        }
    }
    out.sort_by(|a, b| {
        let fa = a.0 as f64 / a.1 as f64;
        let fb = b.0 as f64 / b.1 as f64;
        fa.partial_cmp(&fb).unwrap()
    });
    out
}

/// Number of intermediate m:j resonances (m ≤ `max_m`) between two
/// resonances — the fine-comb density.
pub fn farey_density(m1: i64, j1: i64, m2: i64, j2: i64, max_m: i64) -> usize {
    enumerate_intermediate_resonances(m1, j1, m2, j2, max_m).len()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resonance::build_2j_chain;
    use p9_core::analysis::resonance::chirikov_overlap_parameter;

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
    fn test_overlap_commutative() {
        // M4: the criterion must not depend on argument order. Probe
        // asymmetric width pairs around the overlap threshold.
        for &(d1, d2, sep) in &[
            (5.0, 1.0, 6.5_f64),
            (5.0, 1.0, 7.5),
            (1.0, 5.0, 6.5),
            (0.5, 8.0, 8.7),
            (3.0, 3.0, 5.9),
        ] {
            let r1 = make_res(100.0, d1);
            let r2 = make_res(100.0 + sep, d2);
            assert_eq!(
                resonances_overlap(&r1, &r2),
                resonances_overlap(&r2, &r1),
                "overlap not commutative for widths ({d1}, {d2}) at sep {sep}"
            );
            assert!((overlap_parameter(&r1, &r2) - overlap_parameter(&r2, &r1)).abs() < 1e-12);
        }
    }

    #[test]
    fn test_weaker_resonance_reach_limited() {
        // Wide (5.0) + narrow (1.0): the narrow side contributes its full
        // width (2.0), so the threshold is 7.0 — not the naive 6.0 sum
        // and not the order-dependent 10.0.
        let wide = make_res(100.0, 5.0);
        let narrow_in = make_res(106.9, 1.0);
        let narrow_out = make_res(107.1, 1.0);
        assert!(resonances_overlap(&wide, &narrow_in));
        assert!(!resonances_overlap(&wide, &narrow_out));
    }

    #[test]
    fn test_overlap_parameter() {
        let r1 = make_res(100.0, 5.0);
        let r2 = make_res(110.0, 5.0);
        let k = overlap_parameter(&r1, &r2);
        assert!((k - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_chain_overlap_consistent_with_p9_core_chirikov() {
        // H3 cross-crate consistency: for the 2:j chain the pairwise
        // width-based overlap measure (Δa_j + Δa_{j+1})/(2·spacing)
        // reproduces p9-core's closed-form Chirikov parameter
        // K = (24/√5)(a/a_N)^{5/4}√(m_N) e^{−q²/2a_N²} at the midpoint.
        for &q in &[30.0, 35.0, 40.0] {
            let chain = build_2j_chain(15, 45, q);
            for pair in chain.windows(2) {
                let spacing = pair[1].a_nominal - pair[0].a_nominal;
                let k_chain = (pair[0].delta_a + pair[1].delta_a) / (2.0 * spacing);
                let a_mid = 0.5 * (pair[0].a_nominal + pair[1].a_nominal);
                let k_core = chirikov_overlap_parameter(a_mid, q);
                assert!(
                    (k_chain - k_core).abs() < 5e-3 * k_core,
                    "q={q}, j={}: chain K = {k_chain:.5} vs core K = {k_core:.5}",
                    pair[0].j
                );
            }
        }
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
    fn test_megno_regular_orbit_converges_to_two() {
        // Linear (power-law) divergence δ = 1 + αt → Y(t) → 2: regular.
        let dt = 0.5;
        let deltas: Vec<f64> = (0..20_000).map(|k| 1.0 + 0.1 * (k as f64 * dt)).collect();
        let y = megno_mean_from_separation(&deltas, dt);
        assert!(y > 1.5 && y < 2.05, "⟨Y⟩ = {y:.3}, expected → 2");
        assert_eq!(classify_megno(y.log2()), MegnoClass::Regular);
    }

    #[test]
    fn test_megno_chaotic_orbit_grows() {
        // Exponential divergence δ = e^{λt} → Y(t) = λt: chaotic.
        let (dt, lambda) = (0.5, 0.05);
        let deltas: Vec<f64> = (0..2_000).map(|k| (lambda * k as f64 * dt).exp()).collect();
        let y = megno_mean_from_separation(&deltas, dt);
        // ⟨Y⟩ ≈ λT/2 = 0.05·1000/2 = 25
        assert!(y > 20.0 && y < 30.0, "⟨Y⟩ = {y:.1}, expected ~25");
        assert_eq!(classify_megno(y.log2()), MegnoClass::Chaotic);
    }

    #[test]
    fn test_farey_mediant() {
        let (p, q) = farey_mediant(1, 2, 1, 3);
        assert_eq!((p, q), (2, 5));
    }

    #[test]
    fn test_fine_comb_between_adjacent_2j() {
        // Between 2:26 and 2:25 (ratios 0.0769–0.0800) with m ≤ 6 the
        // complete Farey comb is {3/38, 4/51, 5/64, 5/63, 6/77}.
        let comb = enumerate_intermediate_resonances(2, 26, 2, 25, 6);
        assert_eq!(
            comb,
            vec![(6, 77), (5, 64), (4, 51), (3, 38), (5, 63)],
            "comb = {comb:?}"
        );
        assert_eq!(farey_density(2, 26, 2, 25, 6), 5);

        // All reduced, strictly inside the interval.
        for &(m, j) in &comb {
            assert_eq!(gcd(m, j), 1);
            let f = m as f64 / j as f64;
            assert!(f > 2.0 / 26.0 && f < 2.0 / 25.0);
        }

        // Mediant property: the mediant of the endpoints (reduced) is the
        // lowest-numerator member of the comb.
        let (pm, qm) = farey_mediant(1, 13, 2, 25); // 2/26 reduces to 1/13
        assert_eq!((pm, qm), (3, 38));
        assert!(comb.contains(&(3, 38)));
        let min_m = comb.iter().map(|&(m, _)| m).min().unwrap();
        assert_eq!(min_m, 3);
    }

    #[test]
    fn test_fine_comb_density_grows_with_order() {
        // Higher allowed numerators → denser comb (the paper's point).
        let d4 = farey_density(2, 26, 2, 25, 4);
        let d8 = farey_density(2, 26, 2, 25, 8);
        let d16 = farey_density(2, 26, 2, 25, 16);
        assert!(d4 < d8 && d8 < d16, "densities {d4}, {d8}, {d16}");
    }

    #[test]
    fn test_no_resonances_between_adjacent_fractions() {
        // 1/2 and 1/3 are Farey neighbours: nothing with m = 1 between.
        assert_eq!(farey_density(1, 2, 1, 3, 1), 0);
    }
}
