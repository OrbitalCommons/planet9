//! Posterior over (m, a_p) and its maximum-a-posteriori solution.
//!
//! The forcing likelihood alone is degenerate: any (m, a_p) on the ridge
//! S(m,a_p) = m/a_p³ = S_obs fits the apsidal confinement equally well. Siraj,
//! Chyba & Tremaine (2024) break that degeneracy with an *independent*
//! constraint from the same dynamics — the perturber inclination forces a
//! node confinement and a perihelion-distance scale that, together, localize
//! a_p. We encode that independent constraint as a Gaussian prior on a_p
//! centered on the inclination-derived semi-major axis (the published 290 AU),
//! with a width broad enough that the forcing likelihood still shapes the MAP.
//!
//! Posterior (up to a constant):
//!
//!   −2 ln P(m, a_p) = [ln S(m,a_p) − ln S_obs]² / σ_lnS²
//!                   + [(a_p − a_ref) / σ_a]²
//!
//! The first term is the apsidal-confinement likelihood (a Gaussian in ln S,
//! producing the m ∝ a_p³ ridge); the second is the independent a_p prior.
//! The MAP is found by a numerical grid search — nothing about the answer is
//! hard-coded; `k_cal` is the *only* calibration, fixed so that the observed
//! ETNO confinement maps the ridge through the published reference orbit.

use crate::confinement::observed_forcing;
use crate::forcing::forcing_strength;
use p9_core::data::ephemeris_constraint::{SIRAJ_2024_A_AU, SIRAJ_2024_MASS_EARTH};
use p9_core::units::{au, earth_masses, Length, Mass};

/// A point in the (mass, semi-major axis) inference space.
#[derive(Debug, Clone, Copy)]
pub struct OrbitPoint {
    /// Perturber mass in Earth masses.
    pub mass_earth: f64,
    /// Perturber semi-major axis in AU.
    pub a_p: f64,
}

impl OrbitPoint {
    /// Perturber mass as a dimension-checked [`Mass`].
    pub fn mass(&self) -> Mass {
        earth_masses(self.mass_earth)
    }

    /// Perturber semi-major axis as a dimension-checked [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_p)
    }
}

/// The (m, a_p) posterior model.
#[derive(Debug, Clone)]
pub struct ForcingPosterior {
    /// Required observed forcing strength S_obs (M⊕ / AU³).
    pub s_obs: f64,
    /// 1σ width of the ln-S likelihood (fractional κ̂ uncertainty).
    pub sigma_ln_s: f64,
    /// Center of the independent a_p prior (AU).
    pub a_ref: f64,
    /// 1σ width of the independent a_p prior (AU).
    pub sigma_a: f64,
}

impl ForcingPosterior {
    /// Build the posterior from the ETNO ϖ sample.
    ///
    /// `k_cal` calibrates κ̂ → S_obs (see [`crate::confinement`]); `a_ref` and
    /// `sigma_a` describe the independent inclination-derived a_p constraint.
    pub fn from_sample(angles: &[f64], k_cal: f64, a_ref: f64, sigma_a: f64) -> Self {
        let (s_obs, frac) = observed_forcing(angles, k_cal);
        Self {
            s_obs,
            sigma_ln_s: frac,
            a_ref,
            sigma_a,
        }
    }

    /// −2 ln P at a point (smaller is better).
    pub fn neg2_log_post(&self, p: OrbitPoint) -> f64 {
        let s = forcing_strength(p.mass_earth, p.a_p);
        let d_ln_s = (s / self.s_obs).ln() / self.sigma_ln_s;
        let d_a = (p.a_p - self.a_ref) / self.sigma_a;
        d_ln_s * d_ln_s + d_a * d_a
    }

    /// Log-posterior (up to an additive constant).
    pub fn log_post(&self, p: OrbitPoint) -> f64 {
        -0.5 * self.neg2_log_post(p)
    }

    /// Numerically locate the MAP over a (mass, a_p) grid spanning the search
    /// box, then refine with a local golden-section-style contraction. No
    /// closed-form answer is used.
    pub fn map_estimate(&self) -> OrbitPoint {
        let mut best = OrbitPoint {
            mass_earth: 1.0,
            a_p: self.a_ref,
        };
        let mut best_val = f64::INFINITY;

        // Coarse grid over a physically generous box.
        let mut m_lo = 0.5_f64;
        let mut m_hi = 30.0_f64;
        let mut a_lo = 100.0_f64;
        let mut a_hi = 700.0_f64;

        for _refine in 0..6 {
            let n = 81;
            for i in 0..n {
                let a_p = a_lo + (a_hi - a_lo) * i as f64 / (n - 1) as f64;
                for j in 0..n {
                    let mass_earth = m_lo + (m_hi - m_lo) * j as f64 / (n - 1) as f64;
                    let p = OrbitPoint { mass_earth, a_p };
                    let v = self.neg2_log_post(p);
                    if v < best_val {
                        best_val = v;
                        best = p;
                    }
                }
            }
            // Contract the box around the current best.
            let da = (a_hi - a_lo) * 0.15;
            let dm = (m_hi - m_lo) * 0.15;
            a_lo = (best.a_p - da).max(50.0);
            a_hi = best.a_p + da;
            m_lo = (best.mass_earth - dm).max(0.1);
            m_hi = best.mass_earth + dm;
        }
        best
    }
}

/// Fixed calibration anchor: the von Mises concentration of the clustered
/// ETNO ϖ sample AT THE PAPER'S EPOCH, frozen as a labelled constant
/// (derived from the paper's own headline: 2.7σ Rayleigh significance at
/// n = 21 means p = exp(−nR̄²) ≈ 0.0069, so R̄ ≈ 0.487 and the Mardia-Jupp
/// inverse gives κ ≈ 1.09 — a number taken from the publication, not from
/// this repo's sample).
///
/// This one-point calibration ties (κ_REF → the published reference orbit's
/// forcing m_ref/a_ref³) once. Crucially it is NOT recomputed from the input
/// angles: a previous version solved k_cal from the live sample's κ̂, which
/// made κ̂ cancel identically out of S_obs — the "inference" returned the
/// published orbit for ANY input angles. With the frozen anchor, S_obs = κ̂ ·
/// (m_ref/a_ref³)/κ_REF moves with the data (see the falsifiability tests).
pub const KAPPA_REFERENCE: f64 = 1.09;

/// Build the Siraj et al. (2024) posterior from the clustered-ETNO ϖ sample.
/// The ridge scale is anchored by the frozen one-point calibration
/// [`KAPPA_REFERENCE`] → (4.4 M⊕/290 AU³); the input sample's κ̂ then sets
/// S_obs relative to that anchor, so a weaker-clustered sample infers a
/// proportionally weaker perturber.
pub fn calibrated_posterior(angles: &[f64]) -> ForcingPosterior {
    let k_cal = KAPPA_REFERENCE * SIRAJ_2024_A_AU.powi(3) / SIRAJ_2024_MASS_EARTH;
    // Independent a_p prior: centered on the paper's inclination/node-derived
    // semi-major-axis scale (a labelled paper input, distinct from the ridge
    // normalization), with a 20% width — broad enough that the forcing
    // likelihood meaningfully shapes the MAP, tight enough to break the
    // m ∝ a_p³ degeneracy.
    let sigma_a = 0.20 * SIRAJ_2024_A_AU;
    ForcingPosterior::from_sample(angles, k_cal, SIRAJ_2024_A_AU, sigma_a)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::data::etno::longitudes_of_perihelion;
    use uom::si::length::astronomical_unit;
    use uom::si::mass::kilogram;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let p = OrbitPoint {
            mass_earth: SIRAJ_2024_MASS_EARTH,
            a_p: SIRAJ_2024_A_AU,
        };
        assert_relative_eq!(
            p.semi_major_axis().get::<astronomical_unit>(),
            p.a_p,
            epsilon = 1e-12
        );
        // Mass reads back in Earth masses via the kg ratio.
        let earth_kg = earth_masses(1.0).get::<kilogram>();
        assert_relative_eq!(
            p.mass().get::<kilogram>() / earth_kg,
            p.mass_earth,
            max_relative = 1e-12
        );
    }

    #[test]
    fn test_neg2_log_post_minimized_on_ridge_at_aref() {
        let post = calibrated_posterior(&longitudes_of_perihelion());
        // At a_ref, the best mass is the one on the DATA's ridge
        // (m = S_obs·a_ref³); perturbing mass away raises the objective.
        let m_ridge = crate::forcing::mass_for_strength(post.s_obs, SIRAJ_2024_A_AU);
        let on = OrbitPoint {
            mass_earth: m_ridge,
            a_p: SIRAJ_2024_A_AU,
        };
        let off = OrbitPoint {
            mass_earth: m_ridge * 2.0,
            a_p: SIRAJ_2024_A_AU,
        };
        assert!(post.neg2_log_post(on) < post.neg2_log_post(off));
    }

    #[test]
    fn test_inference_is_falsifiable() {
        use p9_core::constants::TWO_PI;
        // A uniform (unclustered) ϖ sample must NOT return the published
        // orbit: with κ̂ ≈ 0 the inferred mass collapses far below 4.4 M⊕.
        // (The pre-fix calibration returned the published orbit for ANY
        // input — the data cancelled identically.)
        let uniform: Vec<f64> = (0..21).map(|k| k as f64 * TWO_PI / 21.0).collect();
        let map_u = calibrated_posterior(&uniform).map_estimate();
        assert!(
            map_u.mass_earth < 0.5 * SIRAJ_2024_MASS_EARTH,
            "uniform sample inferred {:.2} M⊕",
            map_u.mass_earth
        );

        // A more concentrated sample infers a stronger perturber.
        let tight: Vec<f64> = (0..21).map(|k| 0.9 + 0.02 * k as f64).collect();
        let map_t = calibrated_posterior(&tight).map_estimate();
        let map_obs = calibrated_posterior(&longitudes_of_perihelion()).map_estimate();
        assert!(
            map_t.mass_earth > map_obs.mass_earth,
            "tight {:.2} vs observed {:.2} M⊕",
            map_t.mass_earth,
            map_obs.mass_earth
        );
    }
}
