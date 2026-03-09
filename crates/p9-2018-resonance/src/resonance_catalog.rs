//! Extended mean-motion resonance catalog and Farey sequence generation.
//!
//! Builds on the 11 resonances from Batygin & Brown (2016) with the full
//! spectrum of high-order resonances identified in Bailey+ (2018).
//!
//! Key concept: the Farey sequence F_n contains all fractions p/q with
//! q ≤ n in lowest terms. Millholland & Laughlin (2017) used F_5 (only
//! denominators ≤ 5), while Bailey+ (2018) shows the actual resonance
//! occupation includes many higher-order terms.

use std::f64::consts::PI;

use p9_core::constants::TWO_PI;
use p9_core::types::OrbitalElements;

/// A mean-motion resonance specification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub struct Resonance {
    /// Particle's coefficient (e.g., p in a p:q resonance)
    pub p: u32,
    /// Planet Nine's coefficient
    pub q: u32,
}

impl Resonance {
    pub fn new(p: u32, q: u32) -> Self {
        Self { p, q }
    }

    /// Period ratio T_particle / T_p9 = q/p.
    pub fn period_ratio(&self) -> f64 {
        self.q as f64 / self.p as f64
    }

    /// Semi-major axis of the resonance given perturber semi-major axis.
    /// a_res = a_p9 * (q/p)^(2/3)
    pub fn semimajor_axis(&self, a_p9: f64) -> f64 {
        a_p9 * self.period_ratio().powf(2.0 / 3.0)
    }

    /// Whether this is an N/1 period ratio (integer period ratio, p=1).
    pub fn is_n_over_1(&self) -> bool {
        self.p == 1
    }

    /// Whether this is an N/2 period ratio (half-integer period ratio, p=2).
    pub fn is_n_over_2(&self) -> bool {
        self.p == 2
    }

    /// Whether this is N/1 or N/2 (the "simple" resonances).
    pub fn is_simple(&self) -> bool {
        self.is_n_over_1() || self.is_n_over_2()
    }

    /// Order of the resonance: |p - q|.
    pub fn order(&self) -> u32 {
        if self.p > self.q {
            self.p - self.q
        } else {
            self.q - self.p
        }
    }
}

impl std::fmt::Display for Resonance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}:{}", self.p, self.q)
    }
}

/// Generate the Farey sequence F_n: all fractions p/q in (0, 1] with q ≤ n,
/// in lowest terms, plus interior/exterior extensions.
///
/// For resonances, we want period ratios both < 1 (interior) and > 1 (exterior).
/// Returns resonances with period ratios from 1/max_p to max_p/1.
pub fn farey_resonances(max_denominator: u32, max_p: u32) -> Vec<Resonance> {
    let mut resonances = Vec::new();

    for q in 1..=max_denominator {
        for p in 1..=max_p {
            if gcd(p, q) == 1 && p != q {
                resonances.push(Resonance::new(p, q));
            }
        }
    }

    resonances.sort_by(|a, b| {
        a.period_ratio()
            .partial_cmp(&b.period_ratio())
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    resonances.dedup();
    resonances
}

/// Farey F_5 sequence used by Millholland & Laughlin (2017).
pub fn farey_f5() -> Vec<Resonance> {
    farey_resonances(5, 30)
}

/// Extended resonance catalog including high-order terms (denominators up to 20).
pub fn extended_catalog() -> Vec<Resonance> {
    farey_resonances(20, 35)
}

/// Most populated resonances from Bailey+ (2018) Table/Figure 4.
pub fn most_populated_resonances() -> Vec<Resonance> {
    vec![
        Resonance::new(1, 2), // Exterior 1:2
        Resonance::new(1, 1), // 1:1 (co-orbital)
        Resonance::new(3, 2), // 3:2
        Resonance::new(2, 1), // 2:1
        Resonance::new(3, 1), // 3:1
        Resonance::new(5, 3), // 5:3
        Resonance::new(7, 4), // 7:4
        Resonance::new(10, 1),
        Resonance::new(11, 1),
        Resonance::new(13, 3),
        Resonance::new(13, 4),
        Resonance::new(20, 1),
        Resonance::new(22, 7),
    ]
}

fn gcd(a: u32, b: u32) -> u32 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// Identify which resonance (if any) a particle is near.
///
/// Searches through the given catalog and returns the closest resonance
/// within the fractional tolerance.
pub fn identify_resonance(
    a: f64,
    a_p9: f64,
    catalog: &[Resonance],
    tolerance_frac: f64,
) -> Option<(Resonance, f64)> {
    let mut best: Option<(Resonance, f64)> = None;

    for &res in catalog {
        let a_res = res.semimajor_axis(a_p9);
        let delta = (a - a_res).abs() / a_res;

        if delta < tolerance_frac {
            if best.is_none() || delta < best.unwrap().1 {
                best = Some((res, delta));
            }
        }
    }

    best
}

/// Compute the resonant angle for a p:q resonance.
///
/// φ = p*λ_particle - q*λ_p9 - (p-q)*ϖ_particle
pub fn resonant_angle(
    elem_particle: &OrbitalElements,
    elem_p9: &OrbitalElements,
    res: &Resonance,
) -> f64 {
    let lambda_part = elem_particle.mean_anomaly + elem_particle.omega + elem_particle.omega_big;
    let lambda_p9 = elem_p9.mean_anomaly + elem_p9.omega + elem_p9.omega_big;
    let varpi_part = elem_particle.omega + elem_particle.omega_big;

    let diff = res.p as i64 - res.q as i64;
    let mut phi = res.p as f64 * lambda_part - res.q as f64 * lambda_p9 - diff as f64 * varpi_part;

    phi = phi % TWO_PI;
    if phi > PI {
        phi -= TWO_PI;
    }
    if phi < -PI {
        phi += TWO_PI;
    }
    phi
}

/// Check if a sequence of resonant angles indicates libration.
///
/// Uses circular statistics: r_bar close to 1 means tightly clustered (librating).
pub fn is_librating(angles: &[f64]) -> (bool, f64) {
    if angles.len() < 10 {
        return (false, TWO_PI);
    }

    let sin_sum: f64 = angles.iter().map(|&a| a.sin()).sum();
    let cos_sum: f64 = angles.iter().map(|&a| a.cos()).sum();
    let n = angles.len() as f64;
    let r_bar = (sin_sum * sin_sum + cos_sum * cos_sum).sqrt() / n;

    let is_lib = r_bar > 0.5;
    let amplitude = if r_bar > 0.0 {
        2.0 * (-2.0 * r_bar.ln()).sqrt()
    } else {
        TWO_PI
    };

    (is_lib, amplitude.min(TWO_PI))
}

/// Count how many elements in a snapshot fall near each resonance.
pub fn resonance_census(
    elements: &[OrbitalElements],
    a_p9: f64,
    catalog: &[Resonance],
    tolerance: f64,
) -> Vec<(Resonance, usize)> {
    catalog
        .iter()
        .filter_map(|&res| {
            let a_res = res.semimajor_axis(a_p9);
            let count = elements
                .iter()
                .filter(|e| ((e.a - a_res) / a_res).abs() < tolerance)
                .count();
            if count > 0 {
                Some((res, count))
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_resonance_semimajor_axis() {
        let res = Resonance::new(2, 1);
        let a = res.semimajor_axis(700.0);
        // 2:1 → a = 700 * (1/2)^(2/3) ≈ 441
        assert!((a - 441.0).abs() < 1.0, "2:1 at {:.1} AU", a);
    }

    #[test]
    fn test_resonance_classification() {
        // p=1 → N/1 (integer period ratio)
        assert!(Resonance::new(1, 3).is_n_over_1());
        // p=2 → N/2 (half-integer period ratio)
        assert!(Resonance::new(2, 5).is_n_over_2());
        // p=5 → high-order
        assert!(!Resonance::new(5, 3).is_simple());
        assert_eq!(Resonance::new(5, 3).order(), 2);
    }

    #[test]
    fn test_farey_f5() {
        let f5 = farey_f5();
        assert!(!f5.is_empty());

        // All denominators should be ≤ 5
        for res in &f5 {
            assert!(res.q <= 5, "F5 should have q ≤ 5, got {}", res);
        }

        // Should contain standard resonances
        assert!(f5.contains(&Resonance::new(2, 1)));
        assert!(f5.contains(&Resonance::new(3, 1)));
        assert!(f5.contains(&Resonance::new(3, 2)));
    }

    #[test]
    fn test_extended_catalog_is_larger() {
        let f5 = farey_f5();
        let extended = extended_catalog();
        assert!(
            extended.len() > f5.len(),
            "Extended ({}) should be larger than F5 ({})",
            extended.len(),
            f5.len()
        );
    }

    #[test]
    fn test_gcd() {
        assert_eq!(gcd(12, 8), 4);
        assert_eq!(gcd(7, 3), 1);
        assert_eq!(gcd(6, 6), 6);
    }

    #[test]
    fn test_identify_resonance() {
        let catalog = most_populated_resonances();
        let a_p9 = 600.0;

        // Near 2:1 resonance: a = 600 * (1/2)^(2/3) ≈ 378
        let a_21 = Resonance::new(2, 1).semimajor_axis(a_p9);
        let result = identify_resonance(a_21 + 1.0, a_p9, &catalog, 0.02);
        assert!(result.is_some());
        assert_eq!(result.unwrap().0, Resonance::new(2, 1));

        // Far from any resonance
        let result = identify_resonance(500.0, a_p9, &catalog, 0.005);
        assert!(result.is_none());
    }

    #[test]
    fn test_libration_detection() {
        // Librating angles near 0
        let librating: Vec<f64> = (0..100).map(|i| 0.3 * (i as f64 * 0.1).sin()).collect();
        let (is_lib, amp) = is_librating(&librating);
        assert!(is_lib);
        assert!(amp < PI);

        // Circulating angles
        let circulating: Vec<f64> = (0..100).map(|i| (i as f64 / 100.0 * TWO_PI) - PI).collect();
        let (is_lib, _) = is_librating(&circulating);
        assert!(!is_lib);
    }

    #[test]
    fn test_resonance_census() {
        let catalog = vec![Resonance::new(2, 1), Resonance::new(3, 1)];
        let a_p9 = 600.0;
        let a_21 = Resonance::new(2, 1).semimajor_axis(a_p9);

        let elements = vec![
            OrbitalElements {
                a: a_21 + 0.5,
                e: 0.5,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
            OrbitalElements {
                a: a_21 - 0.5,
                e: 0.6,
                i: 0.0,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            },
        ];

        let census = resonance_census(&elements, a_p9, &catalog, 0.02);
        assert_eq!(census.len(), 1);
        assert_eq!(census[0].0, Resonance::new(2, 1));
        assert_eq!(census[0].1, 2);
    }
}
