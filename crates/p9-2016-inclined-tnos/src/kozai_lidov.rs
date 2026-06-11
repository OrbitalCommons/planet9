//! Kozai-Lidov oscillation detection and pathway analysis.
//!
//! The paper identifies a two-step pathway for generating high-i TNOs:
//! 1. Planet Nine drives Kozai-Lidov oscillations (coupled e-i oscillations)
//! 2. Neptune scattering reduces semi-major axis at high inclination
//!
//! Classification works on per-particle *time series* (see
//! `simulation::particle_histories`), not single-epoch (a, i) cuts: the KL
//! signature is the e–i anti-correlation over time, and the high-i pathway
//! is keyed on a semi-major-axis decrease accompanied by Neptune-encounter
//! proximity. (A previous static classifier also rejected exactly the
//! non-conserving trajectories that produce orbit flips by demanding the
//! Kozai integral stay constant.)

use p9_core::constants::*;
use p9_core::types::OrbitalElements;

/// Time series record for a single particle.
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct ParticleHistory {
    pub id: usize,
    pub times: Vec<f64>,
    pub elements: Vec<OrbitalElements>,
}

/// Classification of a particle's dynamical pathway.
#[derive(Debug, Clone, Copy, PartialEq, serde::Serialize, serde::Deserialize)]
pub enum Pathway {
    /// Still in the scattered disk (no KL signature, no decoupling)
    ScatteredDisk,
    /// Undergoing Kozai-Lidov oscillations (e–i anti-correlated in time)
    KozaiLidov,
    /// Decoupled high-i object: ends at small a and high i after an
    /// a-decrease with Neptune-encounter proximity (Drac/Niku analogs)
    HighInclination,
    /// Ejected or accreted (history is empty or ends before the run)
    Removed,
}

/// Minimum inclination swing (degrees) for a KL detection.
pub const KL_MIN_I_AMPLITUDE_DEG: f64 = 20.0;

/// Minimum |Pearson r| between e and cos i for a KL detection. In a KL
/// cycle e is maximal where i is minimal, so e and cos i are strongly
/// correlated for prograde orbits (anti-correlated e–i); the correlation
/// flips sign across i = 90°, hence the absolute value.
pub const KL_MIN_E_COSI_CORRELATION: f64 = 0.5;

/// Perihelion distance (AU) inside which a particle is in Neptune's
/// dynamical reach (a_N ≈ 30.1 AU plus the strong-scattering zone;
/// the classic scattered-disk criterion q ≲ 38 AU).
pub const NEPTUNE_SCATTERING_Q_AU: f64 = 38.0;

/// Pearson correlation coefficient of two equal-length series.
fn pearson_r(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let mx = x.iter().sum::<f64>() / n;
    let my = y.iter().sum::<f64>() / n;
    let mut sxy = 0.0;
    let mut sxx = 0.0;
    let mut syy = 0.0;
    for (&xi, &yi) in x.iter().zip(y) {
        sxy += (xi - mx) * (yi - my);
        sxx += (xi - mx) * (xi - mx);
        syy += (yi - my) * (yi - my);
    }
    if sxx <= 0.0 || syy <= 0.0 {
        return 0.0;
    }
    sxy / (sxx * syy).sqrt()
}

/// Detect Kozai-Lidov oscillations in a particle's history.
///
/// KL signature: a large inclination swing (> `KL_MIN_I_AMPLITUDE_DEG`)
/// whose e and cos i track each other (|Pearson r| >
/// `KL_MIN_E_COSI_CORRELATION`), i.e. e–i anti-correlation in time.
pub fn detect_kozai_lidov(history: &ParticleHistory) -> bool {
    if history.elements.len() < 10 {
        return false;
    }

    let i_values: Vec<f64> = history.elements.iter().map(|e| e.i * RAD2DEG).collect();
    let i_min = i_values.iter().cloned().fold(f64::INFINITY, f64::min);
    let i_max = i_values.iter().cloned().fold(f64::NEG_INFINITY, f64::max);

    if i_max - i_min < KL_MIN_I_AMPLITUDE_DEG {
        return false;
    }

    let e_values: Vec<f64> = history.elements.iter().map(|e| e.e).collect();
    let cos_i_values: Vec<f64> = history.elements.iter().map(|e| e.i.cos()).collect();

    pearson_r(&e_values, &cos_i_values).abs() > KL_MIN_E_COSI_CORRELATION
}

/// Classify a particle's pathway from its full time series.
///
/// High-inclination decoupling requires all three signatures of the
/// two-step mechanism: a high-i, small-a endpoint, a substantial semi-major
/// axis *decrease* along the way, and perihelion dipping into Neptune's
/// scattering zone (the a-decrease must be Neptune-mediated, not numerical).
pub fn classify_pathway(history: &ParticleHistory) -> Pathway {
    let Some(last) = history.elements.last() else {
        return Pathway::Removed;
    };

    let a_max = history
        .elements
        .iter()
        .map(|e| e.a)
        .fold(f64::NEG_INFINITY, f64::max);
    let q_min = history
        .elements
        .iter()
        .map(|e| e.a * (1.0 - e.e))
        .fold(f64::INFINITY, f64::min);

    let high_i_final = last.i > 50.0 * DEG2RAD;
    let a_decreased = last.a < 0.5 * a_max;
    let neptune_proximity = q_min < NEPTUNE_SCATTERING_Q_AU;

    if last.a < 100.0 && high_i_final && a_decreased && neptune_proximity {
        return Pathway::HighInclination;
    }

    if detect_kozai_lidov(history) {
        return Pathway::KozaiLidov;
    }

    Pathway::ScatteredDisk
}

/// Count particles in each pathway classification over their histories.
///
/// Returns (scattered, kozai_lidov, high_inclination, removed).
pub fn pathway_census(histories: &[ParticleHistory]) -> (usize, usize, usize, usize) {
    let mut scattered = 0;
    let mut kozai = 0;
    let mut high_i = 0;
    let mut removed = 0;

    for history in histories {
        match classify_pathway(history) {
            Pathway::ScatteredDisk => scattered += 1,
            Pathway::KozaiLidov => kozai += 1,
            Pathway::HighInclination => high_i += 1,
            Pathway::Removed => removed += 1,
        }
    }

    (scattered, kozai, high_i, removed)
}

/// Check if a particle matches the orbital region of known high-i TNOs.
///
/// Returns true if a < 100 AU, i > 60°, and q < 40 AU.
pub fn matches_high_i_region(elem: &OrbitalElements) -> bool {
    let q = elem.a * (1.0 - elem.e);
    elem.a < 100.0 && elem.i > 60.0 * DEG2RAD && q < 40.0
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Synthetic KL cycle: i oscillates, e follows from a conserved Kozai
    /// integral (e maximal at i minimal).
    fn kl_history(a: f64, i_mean_deg: f64, i_amp_deg: f64, n: usize) -> ParticleHistory {
        let mut elements = Vec::new();
        let mut times = Vec::new();
        for j in 0..n {
            let phase = TWO_PI * j as f64 / n as f64;
            let i = (i_mean_deg + i_amp_deg * phase.sin()) * DEG2RAD;
            let kozai_const = 0.5;
            let e = (1.0 - (kozai_const / i.cos()).powi(2))
                .sqrt()
                .clamp(0.01, 0.99);
            elements.push(OrbitalElements {
                a,
                e,
                i,
                omega: phase,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            });
            times.push(j as f64 * 1e6 * YEAR_DAYS);
        }
        ParticleHistory {
            id: 0,
            times,
            elements,
        }
    }

    #[test]
    fn test_kozai_lidov_detection() {
        // e–i anti-correlated with a 40° inclination swing: detected.
        let history = kl_history(300.0, 30.0, 20.0, 50);
        assert!(detect_kozai_lidov(&history));
        assert_eq!(classify_pathway(&history), Pathway::KozaiLidov);
    }

    #[test]
    fn test_uncorrelated_oscillation_rejected() {
        // Same i swing, but e oscillates at twice the frequency (uncorrelated
        // with cos i): no KL signature even though the old amplitude-only
        // cut would have fired.
        let n = 50;
        let mut elements = Vec::new();
        let mut times = Vec::new();
        for j in 0..n {
            let phase = TWO_PI * j as f64 / n as f64;
            elements.push(OrbitalElements {
                a: 300.0,
                e: 0.5 + 0.3 * (2.0 * phase).sin(),
                i: (30.0 + 20.0 * phase.sin()) * DEG2RAD,
                omega: phase,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            });
            times.push(j as f64 * 1e6 * YEAR_DAYS);
        }
        let history = ParticleHistory {
            id: 0,
            times,
            elements,
        };
        assert!(!detect_kozai_lidov(&history));
        assert_eq!(classify_pathway(&history), Pathway::ScatteredDisk);
    }

    #[test]
    fn test_small_amplitude_rejected() {
        // Correlated but tiny i swing: not KL.
        let history = kl_history(300.0, 30.0, 4.0, 50);
        assert!(!detect_kozai_lidov(&history));
    }

    #[test]
    fn test_high_inclination_pathway() {
        // Two-step mechanism: KL raises i, then a Neptune encounter
        // (q dips below 38 AU) cuts a from 300 to 45 AU at high i.
        let mut history = kl_history(300.0, 60.0, 25.0, 40);
        for j in 0..10 {
            let frac = j as f64 / 9.0;
            let a = 300.0 - 255.0 * frac; // a: 300 → 45 AU
            let e = 1.0 - 32.0 / a; // q held at 32 AU < 38 AU
            history.elements.push(OrbitalElements {
                a,
                e,
                i: 100.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            });
            history.times.push((40 + j) as f64 * 1e6 * YEAR_DAYS);
        }
        assert_eq!(classify_pathway(&history), Pathway::HighInclination);
    }

    #[test]
    fn test_high_i_without_neptune_proximity_is_not_decoupled() {
        // Ends at high i and small-ish a, but the perihelion never enters
        // Neptune's scattering zone: cannot be the Neptune-assisted pathway.
        let n = 20;
        let mut elements = Vec::new();
        let mut times = Vec::new();
        for j in 0..n {
            let frac = j as f64 / (n - 1) as f64;
            let a = 300.0 - 220.0 * frac;
            elements.push(OrbitalElements {
                a,
                e: 0.3, // q = 0.7a ≥ 56 AU throughout
                i: 80.0 * DEG2RAD,
                omega: 0.0,
                omega_big: 0.0,
                mean_anomaly: 0.0,
            });
            times.push(j as f64 * 1e6 * YEAR_DAYS);
        }
        let history = ParticleHistory {
            id: 0,
            times,
            elements,
        };
        assert_ne!(classify_pathway(&history), Pathway::HighInclination);
    }

    #[test]
    fn test_empty_history_is_removed() {
        let history = ParticleHistory {
            id: 0,
            times: vec![],
            elements: vec![],
        };
        assert_eq!(classify_pathway(&history), Pathway::Removed);
    }

    #[test]
    fn test_pathway_census() {
        let kl = kl_history(300.0, 30.0, 20.0, 50);
        let quiet = ParticleHistory {
            id: 1,
            times: (0..20).map(|j| j as f64 * 1e6 * YEAR_DAYS).collect(),
            elements: (0..20)
                .map(|_| OrbitalElements {
                    a: 300.0,
                    e: 0.7,
                    i: 20.0 * DEG2RAD,
                    omega: 0.0,
                    omega_big: 0.0,
                    mean_anomaly: 0.0,
                })
                .collect(),
        };
        let removed = ParticleHistory {
            id: 2,
            times: vec![],
            elements: vec![],
        };

        let (scattered, kozai, high_i, gone) = pathway_census(&[kl, quiet, removed]);
        assert_eq!(scattered, 1);
        assert_eq!(kozai, 1);
        assert_eq!(high_i, 0);
        assert_eq!(gone, 1);
    }

    #[test]
    fn test_matches_high_i_region() {
        let drac = OrbitalElements {
            a: 41.4,
            e: 0.49,
            i: 103.4 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert!(matches_high_i_region(&drac));

        let distant = OrbitalElements {
            a: 300.0,
            e: 0.7,
            i: 100.0 * DEG2RAD,
            omega: 0.0,
            omega_big: 0.0,
            mean_anomaly: 0.0,
        };
        assert!(!matches_high_i_region(&distant));
    }
}
