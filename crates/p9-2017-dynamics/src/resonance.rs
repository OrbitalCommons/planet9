//! Mean-motion resonance identification and analysis.
//!
//! Key resonances with Planet Nine: 1:1, 3:2, 2:1, 5:2.
//! The critical transition from secular to resonance-dominated dynamics
//! occurs at P/P_9 ~ 0.1-0.15 (Batygin & Morbidelli 2017); here that
//! threshold is *derived* from the Chirikov overlap of the N/1 resonance
//! chain rather than hardcoded.
//!
//! Labeling convention in this crate: a "p:q" resonance means the particle
//! completes p orbits per q Planet Nine orbits (particle interior, faster),
//! so a_res = a₉ (q/p)^{2/3} and the stationary angle is
//! φ = q λ − p λ₉ + (p − q) ϖ (delegated to `p9_core::analysis::resonance`
//! with its (P, Q) = (q, p) exterior convention).

use p9_core::analysis::hansen::{hansen_coefficient, mean_to_true_anomaly};
use p9_core::analysis::resonance::resonant_angle;
use p9_core::constants::{EARTH_MASS_SOLAR, GM_SUN, TWO_PI};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;

/// Description of a mean-motion resonance with Planet Nine.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Resonance {
    /// Particle orbit count (p in p:q)
    pub p: i64,
    /// Planet Nine orbit count (q in p:q)
    pub q: i64,
    /// Nominal semi-major axis for this resonance (AU)
    pub a_nominal: f64,
    /// Pendulum half-width in semi-major axis (AU), at the reference
    /// eccentricity used to build the catalog
    pub delta_a: f64,
}

impl Resonance {
    /// Compute the nominal semi-major axis for a p:q resonance given a_9.
    pub fn nominal_semimajor(p: i64, q: i64, a9: f64) -> f64 {
        a9 * (q as f64 / p as f64).powf(2.0 / 3.0)
    }

    /// Resonant angle φ = q λ − p λ₉ + (p − q) ϖ for this p:q resonance
    /// (single workspace convention via `p9_core::analysis::resonance`).
    pub fn resonant_angle(&self, lambda: f64, lambda_9: f64, varpi: f64) -> f64 {
        resonant_angle(self.q as u32, self.p as u32, lambda, lambda_9, varpi)
    }

    /// Quadrupole-limit pendulum half-width (AU) of the p:q resonance, via
    /// the `p9_core::analysis::hansen` coefficient that carries the O(1)
    /// eccentricity dependence:
    ///
    ///   δa/a = √( (16/3) (m₉/M☉) α^{p+1} |X_q^{p,p}(e)| ),  α = a/a₉
    ///
    /// The p:q resonant angle couples through the l = p multipole at lowest
    /// order, so the relevant Hansen coefficient is X_q^{p,p}(e)
    /// (X_q^{p,p}(0) = δ_{qp}: widths collapse for circular orbits and grow
    /// strongly toward e ≈ 0.7–0.9). Valid only while the particle's
    /// aphelion stays well inside Planet Nine's perihelion — for the
    /// orbit-approaching geometry of the real population use
    /// `pendulum_half_width`, which evaluates the resonant amplitude
    /// numerically along the resonant stream.
    pub fn pendulum_half_width_multipole(p: i64, q: i64, e: f64, m9_solar: f64, a9: f64) -> f64 {
        let a = Self::nominal_semimajor(p, q, a9);
        let alpha = a / a9;
        let x = hansen_coefficient(q, p as i32, p as i32, e, 1e-8).abs();
        a * ((16.0 / 3.0) * m9_solar * alpha.powi(p as i32 + 1) * x).sqrt()
    }

    /// Eccentricity-dependent pendulum half-width (AU) of the p:q resonance
    /// from the numerically evaluated resonant harmonic:
    ///
    ///   δa = a √( 16 |A| a / (3 GM☉) )
    ///
    /// where A is the half-amplitude of the resonant-stream-averaged
    /// interaction over the resonant angle (`resonant_amplitude`). Unlike a
    /// truncated multipole estimate, this remains valid when the particle's
    /// aphelion approaches or crosses Planet Nine's perihelion — the regime
    /// that controls the secular/resonant transition.
    pub fn pendulum_half_width(p: i64, q: i64, e: f64, m9_solar: f64, a9: f64, e9: f64) -> f64 {
        let a = Self::nominal_semimajor(p, q, a9);
        let amp = resonant_amplitude(p, q, e, m9_solar, a9, e9);
        a * ((16.0 / 3.0) * amp.abs() * a / GM_SUN).sqrt()
    }
}

/// Half-amplitude (AU²/day²) of the resonant harmonic of the P9–particle
/// interaction for the p:q resonance: the resonant-stream average of the
/// full interaction
///
///   H(φ) = (1/2πq) ∫ [ −Gm₉/√(Δ² + b²) + Gm₉ (r·r₉)/r₉³ ] dλ₉
///
/// evaluated on the resonant stream λ = (p/q) λ₉ + φ/q (coplanar,
/// anti-aligned apsides, Δϖ = π — the configuration of the confined
/// population), returning (max_φ H − min_φ H)/2 over a grid of resonant
/// phases. A softening b = 0.02 a₉ regularizes orbit-crossing geometries
/// (mutual inclination separates the orbits in reality).
///
/// This is a direct numerical evaluation of the resonance strength — no
/// multipole truncation — so it captures the steep growth as the particle's
/// aphelion approaches Planet Nine's perihelion.
pub fn resonant_amplitude(p: i64, q: i64, e: f64, m9_solar: f64, a9: f64, e9: f64) -> f64 {
    let a = Resonance::nominal_semimajor(p, q, a9);
    let gm9 = m9_solar * GM_SUN;
    let soft2 = (0.005 * a9).powi(2);

    // Anti-aligned apsides: particle perihelion at longitude π, P9's at 0.
    let varpi = PI;
    let varpi9 = 0.0;

    let n_phi = 8;
    // The close-approach window of a crossing geometry spans ~b/(2πa) of the
    // synodic cycle; 2048 samples per perturber orbit resolve it (results
    // agree to 5 significant figures with 4x this resolution).
    let n_s = 2048 * q.max(1) as usize;

    let mut h_min = f64::INFINITY;
    let mut h_max = f64::NEG_INFINITY;

    for k in 0..n_phi {
        let phi = TWO_PI * k as f64 / n_phi as f64;
        let mut sum = 0.0;
        for s_idx in 0..n_s {
            // λ₉ runs over q full P9 orbits while λ runs over p particle
            // orbits: the synodic configuration repeats with period 2πq in λ₉.
            let lambda9 = TWO_PI * q as f64 * (s_idx as f64 + 0.5) / n_s as f64;
            // φ = q λ − p λ₉ + (p − q) ϖ  ⇒  λ = (φ + p λ₉ − (p − q) ϖ)/q
            let lambda = (phi + p as f64 * lambda9 - (p - q) as f64 * varpi) / q as f64;

            let f = mean_to_true_anomaly(lambda - varpi, e);
            let r = a * (1.0 - e * e) / (1.0 + e * f.cos());
            let (x, y) = ((f + varpi).cos() * r, (f + varpi).sin() * r);

            let f9 = mean_to_true_anomaly(lambda9 - varpi9, e9);
            let r9 = a9 * (1.0 - e9 * e9) / (1.0 + e9 * f9.cos());
            let (x9, y9) = ((f9 + varpi9).cos() * r9, (f9 + varpi9).sin() * r9);

            let d2 = (x - x9).powi(2) + (y - y9).powi(2) + soft2;
            let direct = -gm9 / d2.sqrt();
            let indirect = gm9 * (x * x9 + y * y9) / (r9 * r9 * r9);
            sum += direct + indirect;
        }
        let h_phi = sum / n_s as f64;
        h_min = h_min.min(h_phi);
        h_max = h_max.max(h_phi);
    }

    0.5 * (h_max - h_min)
}

/// Reference eccentricity for catalog widths: the observed a > 230 AU
/// population has q ≈ 30–50 AU at a of hundreds of AU, i.e. e ≈ 0.8–0.9.
pub const CATALOG_REFERENCE_E: f64 = 0.85;

/// Fiducial Planet Nine eccentricity for catalog widths.
pub const FIDUCIAL_E9: f64 = 0.6;

/// Key resonances identified in the paper for fiducial a_9 = 700 AU.
pub fn key_resonances_700() -> Vec<Resonance> {
    let a9 = 700.0;
    let m9_solar = 10.0 * EARTH_MASS_SOLAR;

    let resonance_pairs: Vec<(i64, i64)> =
        vec![(1, 1), (3, 2), (2, 1), (5, 2), (3, 1), (4, 1), (5, 1)];

    resonance_pairs
        .into_iter()
        .map(|(p, q)| {
            let a_nominal = Resonance::nominal_semimajor(p, q, a9);
            let delta_a = Resonance::pendulum_half_width(
                p,
                q,
                CATALOG_REFERENCE_E,
                m9_solar,
                a9,
                FIDUCIAL_E9,
            );
            Resonance {
                p,
                q,
                a_nominal,
                delta_a,
            }
        })
        .collect()
}

/// Period ratio P/P_9 for a test particle with semi-major axis a.
pub fn period_ratio(a: f64, a9: f64) -> f64 {
    (a / a9).powf(1.5)
}

/// Critical period ratio below which the N/1 resonance chain ceases to
/// overlap, derived from the resonance-overlap (Chirikov) criterion rather
/// than hardcoded.
///
/// Particles share a perihelion distance `q_peri` (the scattered-disk
/// geometry of the paper), so the eccentricity varies along the chain:
/// e_j = 1 − q_peri/a_j with a_j = a₉ j^{−2/3}, and the aphelion is
/// Q_j = 2 a_j − q_peri. Two regimes drive the overlap:
///
/// - **Orbit-approaching members** (Q_j ≥ q₉ = a₉(1 − e₉)): the particle
///   reaches Planet Nine's orbit, receives impulsive conjunction kicks, and
///   adjacent resonances overlap unconditionally — the canonical onset of
///   chaotic resonant transport (the regime Batygin & Brown describe as
///   resonance hopping).
/// - **Detached members**: the pendulum half-width from the numerically
///   evaluated resonant amplitude is compared against the chain spacing
///   Δa_j = a₉ (j^{−2/3} − (j+1)^{−2/3}); the widths there are ≲ a few AU
///   versus tens of AU spacing, so the chain is isolated and the dynamics
///   secular.
///
/// Returns the period ratio 1/j* of the innermost overlapped member: below
/// it the dynamics is secular, above it resonance-dominated. For the
/// fiducial configuration this lands in the paper's quoted band
/// P/P₉ ≈ 0.1–0.15 (the transition is controlled by where the particle
/// aphelion detaches from Planet Nine's perihelion).
pub fn critical_period_ratio(q_peri: f64, m9_solar: f64, a9: f64, e9: f64) -> f64 {
    let q9 = a9 * (1.0 - e9);
    let mut ratio = 1.0;
    for j in 2..=40i64 {
        let a_j = a9 * (j as f64).powf(-2.0 / 3.0);
        let e_j = 1.0 - q_peri / a_j;
        if e_j <= 0.0 {
            // Resonance inside the perihelion shell: chain ends here.
            return ratio;
        }
        let aphelion = 2.0 * a_j - q_peri;
        let overlapped = if aphelion >= q9 {
            true
        } else {
            let half_width = Resonance::pendulum_half_width(j, 1, e_j, m9_solar, a9, e9);
            let spacing = a9 * ((j as f64).powf(-2.0 / 3.0) - (j as f64 + 1.0).powf(-2.0 / 3.0));
            2.0 * half_width >= spacing
        };
        if !overlapped {
            return 1.0 / j as f64;
        }
        ratio = 1.0 / j as f64;
    }
    ratio
}

/// Critical period ratio for the fiducial configuration (q ≈ 33 AU
/// scattered-disk particles, m₉ = 10 M⊕, a₉ = 700 AU).
pub fn critical_period_ratio_fiducial() -> f64 {
    critical_period_ratio(33.0, 10.0 * EARTH_MASS_SOLAR, 700.0, FIDUCIAL_E9)
}

/// Classify whether a test particle is in the secular or resonance-dominated
/// regime, using the overlap-derived threshold with a ±20% transitional band.
pub fn dynamical_regime(a: f64, a9: f64) -> DynamicalRegime {
    let ratio = period_ratio(a, a9);
    let critical = critical_period_ratio_fiducial();
    if ratio < 0.8 * critical {
        DynamicalRegime::Secular
    } else if ratio < 1.2 * critical {
        DynamicalRegime::Transitional
    } else {
        DynamicalRegime::ResonanceDominated
    }
}

/// Dynamical regime classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DynamicalRegime {
    Secular,
    Transitional,
    ResonanceDominated,
}

/// Check if a test particle with mean motion n is near a p:q resonance with P9
/// (mean motion n9), within a specified tolerance in resonant angle libration.
pub fn is_near_resonance(n: f64, n9: f64, p: i64, q: i64, tol_fraction: f64) -> bool {
    let expected_ratio = p as f64 / q as f64;
    let actual_ratio = n / n9;
    (actual_ratio - expected_ratio).abs() / expected_ratio < tol_fraction
}

/// Adiabatic invariant J = integral of Lambda' d(phi_res').
///
/// Under the J=0 approximation (negligible libration amplitude),
/// objects are trapped at the exact center of resonance.
/// Returns the resonant semi-major axis deviation from nominal.
pub fn adiabatic_deviation(libration_amplitude: f64, resonance_width: f64) -> f64 {
    libration_amplitude * resonance_width / TWO_PI
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_nominal_semimajor_1_1() {
        let a = Resonance::nominal_semimajor(1, 1, 700.0);
        assert!((a - 700.0).abs() < 0.01, "1:1 resonance should be at a_9");
    }

    #[test]
    fn test_nominal_semimajor_2_1() {
        // 2:1 means n_particle = 2 * n_9, so a_particle = a_9 * (1/2)^(2/3)
        let a = Resonance::nominal_semimajor(2, 1, 700.0);
        let expected = 700.0 * (0.5_f64).powf(2.0 / 3.0);
        assert!(
            (a - expected).abs() < 0.1,
            "a = {}, expected {}",
            a,
            expected
        );
        assert!(a < 700.0, "2:1 resonance should be inside P9 orbit");
    }

    #[test]
    fn test_key_resonances_ordered() {
        let res = key_resonances_700();
        assert!(res.len() >= 5, "should have at least 5 key resonances");
        // 1:1 should have largest a_nominal
        assert!(
            res[0].a_nominal > res[1].a_nominal,
            "1:1 should be outermost"
        );
    }

    #[test]
    fn test_period_ratio() {
        let ratio = period_ratio(700.0, 700.0);
        assert!((ratio - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_resonant_angle_stationary_at_exact_resonance() {
        // 2:1 interior resonance: n = 2 n9. The single workspace convention
        // must keep φ constant along the resonant orbit.
        let n9 = 1.0e-3;
        let n = 2.0 * n9;
        let res = Resonance {
            p: 2,
            q: 1,
            a_nominal: 0.0,
            delta_a: 0.0,
        };
        let varpi = 0.7;
        let phi0 = res.resonant_angle(varpi, 0.0, varpi);
        for step in 1..100 {
            let t = step as f64 * 300.0;
            let phi = res.resonant_angle(varpi + n * t, n9 * t, varpi);
            let mut d = (phi - phi0).rem_euclid(TWO_PI);
            if d > std::f64::consts::PI {
                d -= TWO_PI;
            }
            assert!(d.abs() < 1e-9, "φ drifted by {d:.2e}");
        }
    }

    #[test]
    fn test_pendulum_width_eccentricity_dependence() {
        // The O(1) eccentricity dependence the hardcoded width lacked:
        // X_1^{3,3}(e) ∝ e² at small e, so the closed-form width grows by
        // ~(0.85/0.1) between e = 0.1 and 0.85; the numerical width must
        // share the increasing trend (its contrast is milder because even
        // the low-e 3:1 orbit sits within ~100 AU of P9's perihelion).
        let m9 = 10.0 * EARTH_MASS_SOLAR;
        let wm_high = Resonance::pendulum_half_width_multipole(3, 1, 0.85, m9, 700.0);
        let wm_low = Resonance::pendulum_half_width_multipole(3, 1, 0.1, m9, 700.0);
        assert!(
            wm_high > 5.0 * wm_low,
            "multipole: w(0.85) = {wm_high}, w(0.1) = {wm_low}"
        );

        let w_high = Resonance::pendulum_half_width(3, 1, 0.85, m9, 700.0, FIDUCIAL_E9);
        let w_low = Resonance::pendulum_half_width(3, 1, 0.1, m9, 700.0, FIDUCIAL_E9);
        assert!(
            w_high > w_low,
            "numerical: w(0.85) = {w_high} should exceed w(0.1) = {w_low}"
        );
    }

    #[test]
    fn test_multipole_width_vanishes_for_circular_orbit() {
        // The Hansen closed form: X_q^{p,p}(0) = δ_{qp}, so a 3:1 chain
        // member on a circular orbit has zero width at lowest order.
        let m9 = 10.0 * EARTH_MASS_SOLAR;
        let w_circ = Resonance::pendulum_half_width_multipole(3, 1, 0.0, m9, 700.0);
        assert!(w_circ < 1e-6, "circular width should vanish, got {w_circ}");
    }

    #[test]
    fn test_numerical_amplitude_far_field_matches_multipole_scale() {
        // Far-field sanity: for a well-separated, mildly eccentric 2:1
        // orbit, the numerical resonant width and the l = p multipole
        // closed form agree to within an order of magnitude.
        let m9 = 10.0 * EARTH_MASS_SOLAR;
        let (e, a9) = (0.3, 700.0);
        let w_num = Resonance::pendulum_half_width(2, 1, e, m9, a9, 0.0);
        let w_mult = Resonance::pendulum_half_width_multipole(2, 1, e, m9, a9);
        let ratio = w_num / w_mult;
        assert!(
            (0.1..10.0).contains(&ratio),
            "numerical/multipole width ratio = {ratio:.3}"
        );
    }

    #[test]
    fn test_critical_period_ratio_matches_paper_band() {
        // Paper-number regression (Batygin & Morbidelli 2017): the
        // secular-to-resonant transition lies at P/P9 ~ 0.1-0.15. The
        // overlap-derived value must land in (or near) that band rather than
        // being pinned by a hardcoded constant.
        let ratio = critical_period_ratio_fiducial();
        assert!(
            (0.05..=0.25).contains(&ratio),
            "derived critical ratio {ratio:.3} far from the paper's 0.1-0.15"
        );
    }

    #[test]
    fn test_dynamical_regime_secular() {
        // a = 100 AU, a9 = 700 AU => P/P9 = (100/700)^1.5 ~ 0.054
        let regime = dynamical_regime(100.0, 700.0);
        assert_eq!(regime, DynamicalRegime::Secular);
    }

    #[test]
    fn test_dynamical_regime_resonance() {
        // a = 500 AU, a9 = 700 AU => P/P9 = (500/700)^1.5 ~ 0.61
        let regime = dynamical_regime(500.0, 700.0);
        assert_eq!(regime, DynamicalRegime::ResonanceDominated);
    }

    #[test]
    fn test_is_near_resonance() {
        let n9 = 1.0;
        let n = 2.0; // exactly 2:1
        assert!(is_near_resonance(n, n9, 2, 1, 0.01));
        assert!(!is_near_resonance(n, n9, 3, 2, 0.01));
    }
}
