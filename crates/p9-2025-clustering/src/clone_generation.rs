//! TNO clone generation from observational uncertainties.
//!
//! 25 clones per TNO (paper scale), orbital elements sampled around the
//! measured values. Dataset: the vetted `p9_core::data::etno` table
//! (Brown 2017 / JPL SBDB cross-checked) filtered to the paper's
//! selection (a > 150 AU or a + σ_a > 150 AU, σ_a ≤ 100 AU) — replacing
//! the previous 7 hard-coded TNOs.
//!
//! Uncertainty model (documented assumption — the vetted table carries no
//! covariances): σ_a = 2% of a, σ_e from the (a, e) correlation below,
//! σ_i = σ_ω = σ_Ω = 0.1°. ETNO angles are tightly constrained by
//! astrometry relative to a and e.
//!
//! Correlated (a, e) draws: for a multi-opposition ETNO the perihelion
//! distance q is far better constrained than a, so the dominant error
//! mode moves along the q ≈ const curve: we draw δa, set
//! e = 1 − q₀/a + N(0, σ_e) with small independent jitter σ_e = 0.005.
//!
//! Clones that would cross Neptune (q < a_N) when the parent does not are
//! *resampled* rather than silently clamped; if resampling fails after
//! 1000 draws the clone keeps the parent's (a, e) and is flagged via
//! `neptune_crossing` (no silent class changes).

use p9_core::constants::{A_NEPTUNE_AU, DEG2RAD, TWO_PI};
use p9_core::data::etno::BROWN_2017_SAMPLE;
use p9_core::types::OrbitalElements;
use p9_core::units::{au, radians, Angle, Length};
use rand::Rng;
use rand_distr::{Distribution, Normal, Uniform};
use serde::{Deserialize, Serialize};

/// Fractional semi-major-axis uncertainty assumed for the vetted sample.
const SIGMA_A_FRAC: f64 = 0.02;
/// Independent eccentricity jitter on top of the q-preserving correlation.
const SIGMA_E: f64 = 0.005;
/// Angle uncertainties (radians): 0.1° for i, ω, Ω.
const SIGMA_ANGLE: f64 = 0.1 * DEG2RAD;

/// Observed TNO with uncertainties.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ObservedTno {
    pub name: &'static str,
    /// Semi-major axis (AU)
    pub a: f64,
    /// Uncertainty in a (AU)
    pub sigma_a: f64,
    /// Eccentricity
    pub e: f64,
    /// Independent eccentricity jitter (the dominant e error is the
    /// q-preserving correlation with a)
    pub sigma_e: f64,
    /// Inclination (rad)
    pub i: f64,
    /// Uncertainty in i (rad)
    pub sigma_i: f64,
    /// Argument of perihelion (rad)
    pub omega: f64,
    /// Longitude of ascending node (rad)
    pub omega_big: f64,
    /// Uncertainty applied to ω and Ω (rad)
    pub sigma_angle: f64,
}

impl ObservedTno {
    /// Semi-major axis as a dimension-checked [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }

    /// Perihelion distance `q = a(1 − e)` as a dimension-checked [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        au(self.a * (1.0 - self.e))
    }

    /// Inclination as a dimension-checked [`Angle`].
    pub fn inclination_typed(&self) -> Angle {
        radians(self.i)
    }

    /// Longitude of perihelion ϖ = ω + Ω as a dimension-checked [`Angle`],
    /// wrapped to [0, 2π).
    pub fn varpi_typed(&self) -> Angle {
        radians((self.omega + self.omega_big).rem_euclid(TWO_PI))
    }
}

/// A clone of a TNO with sampled orbital elements (all six, including the
/// angles and a uniform mean anomaly — the previous version copied the
/// parent's angles verbatim).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TnoClone {
    pub parent_name: &'static str,
    pub clone_index: usize,
    pub a: f64,
    pub e: f64,
    pub i: f64,
    /// Argument of perihelion (rad)
    pub omega: f64,
    /// Longitude of ascending node (rad)
    pub omega_big: f64,
    /// Mean anomaly (rad), drawn uniformly (unconstrained by the table)
    pub mean_anomaly: f64,
    /// True when resampling could not keep the clone outside Neptune's
    /// orbit even though the parent is (flagged, never silently clamped)
    pub neptune_crossing: bool,
}

impl TnoClone {
    /// Perihelion distance `q = a(1 − e)` as a dimension-checked [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        au(self.a * (1.0 - self.e))
    }

    /// Longitude of perihelion ϖ = ω + Ω as a dimension-checked [`Angle`],
    /// wrapped to [0, 2π).
    pub fn varpi_typed(&self) -> Angle {
        radians((self.omega + self.omega_big).rem_euclid(TWO_PI))
    }

    /// Full orbital element set for integration.
    pub fn elements(&self) -> OrbitalElements {
        OrbitalElements {
            a: self.a,
            e: self.e,
            i: self.i,
            omega: self.omega,
            omega_big: self.omega_big,
            mean_anomaly: self.mean_anomaly,
        }
    }
}

/// Generate clones for a single TNO (see module docs for the sampling
/// model: correlated (a, e), perturbed angles, resampling instead of
/// clamping).
pub fn generate_clones<R: Rng>(tno: &ObservedTno, n_clones: usize, rng: &mut R) -> Vec<TnoClone> {
    let a_dist = Normal::new(0.0, tno.sigma_a.max(1e-6)).unwrap();
    let e_dist = Normal::new(0.0, tno.sigma_e.max(1e-9)).unwrap();
    let i_dist = Normal::new(0.0, tno.sigma_i.max(1e-9)).unwrap();
    let ang_dist = Normal::new(0.0, tno.sigma_angle.max(1e-9)).unwrap();
    let m_dist = Uniform::new(0.0, TWO_PI);

    let q0 = (tno.perihelion_typed() / p9_core::units::au(1.0)).value;
    let parent_crossing = q0 < A_NEPTUNE_AU;

    let mut clones = Vec::with_capacity(n_clones);

    for idx in 0..n_clones {
        let mut accepted: Option<(f64, f64)> = None;
        for _attempt in 0..1000 {
            let a = tno.a + a_dist.sample(rng);
            // (a, e) correlation: preserve q, plus small independent jitter
            let e = 1.0 - q0 / a + e_dist.sample(rng);
            let valid_shape = a > 0.0 && e > 0.0 && e < 0.999;
            let crosses = a * (1.0 - e) < A_NEPTUNE_AU;
            if valid_shape && (parent_crossing || !crosses) {
                accepted = Some((a, e));
                break;
            }
        }
        // Resampling exhausted: keep the parent's (a, e) and flag.
        let (a, e, neptune_crossing) = match accepted {
            Some((a, e)) => (a, e, false),
            None => (tno.a, tno.e, !parent_crossing),
        };

        // Inclination: redraw rather than clamp at 0.
        let mut i = tno.i + i_dist.sample(rng);
        while i < 0.0 {
            i = tno.i + i_dist.sample(rng);
        }

        clones.push(TnoClone {
            parent_name: tno.name,
            clone_index: idx,
            a,
            e,
            i,
            omega: (tno.omega + ang_dist.sample(rng)).rem_euclid(TWO_PI),
            omega_big: (tno.omega_big + ang_dist.sample(rng)).rem_euclid(TWO_PI),
            mean_anomaly: m_dist.sample(rng),
            neptune_crossing,
        });
    }

    clones
}

/// Selection criteria from the paper: a > 150 AU (or a + σ_a > 150 AU),
/// excluding σ_a > 100 AU.
pub fn passes_selection(tno: &ObservedTno) -> bool {
    (tno.a > 150.0 || tno.a + tno.sigma_a > 150.0) && tno.sigma_a <= 100.0
}

/// Sample of distant TNOs used in the analysis: the vetted
/// `p9_core::data::etno` table with the documented uncertainty model,
/// filtered to the paper's selection.
pub fn distant_tno_sample() -> Vec<ObservedTno> {
    BROWN_2017_SAMPLE
        .iter()
        .map(|k| ObservedTno {
            name: k.name,
            a: k.a,
            sigma_a: SIGMA_A_FRAC * k.a,
            e: k.e,
            sigma_e: SIGMA_E,
            i: k.i_deg * DEG2RAD,
            sigma_i: SIGMA_ANGLE,
            omega: k.omega_deg * DEG2RAD,
            omega_big: k.omega_big_deg * DEG2RAD,
            sigma_angle: SIGMA_ANGLE,
        })
        .filter(passes_selection_ref)
        .collect()
}

fn passes_selection_ref(tno: &ObservedTno) -> bool {
    passes_selection(tno)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;
    use uom::si::angle::radian;
    use uom::si::length::astronomical_unit;

    #[test]
    fn typed_accessors_match_f64_sources() {
        let tno = &distant_tno_sample()[0];
        assert_relative_eq!(
            tno.semi_major_axis().get::<astronomical_unit>(),
            tno.a,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            tno.perihelion_typed().get::<astronomical_unit>(),
            tno.a * (1.0 - tno.e),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            tno.inclination_typed().get::<radian>(),
            tno.i,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            tno.varpi_typed().get::<radian>(),
            (tno.omega + tno.omega_big).rem_euclid(TWO_PI),
            epsilon = 1e-12
        );

        let mut rng = rand::rngs::StdRng::seed_from_u64(1);
        let clone = &generate_clones(tno, 1, &mut rng)[0];
        assert_relative_eq!(
            clone.perihelion_typed().get::<astronomical_unit>(),
            clone.a * (1.0 - clone.e),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            clone.varpi_typed().get::<radian>(),
            (clone.omega + clone.omega_big).rem_euclid(TWO_PI),
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_sample_from_vetted_table() {
        let sample = distant_tno_sample();
        // The Brown (2017) table has 10 ETNOs; all pass the a > 150 AU
        // selection (smallest is 2000 CR105 at 228 AU).
        assert_eq!(sample.len(), 10, "expanded sample should have 10 TNOs");
        let names: Vec<&str> = sample.iter().map(|t| t.name).collect();
        assert!(names.contains(&"Sedna"));
        assert!(names.contains(&"2013 RF98"));
    }

    #[test]
    fn test_all_pass_selection() {
        for tno in distant_tno_sample() {
            assert!(
                passes_selection(&tno),
                "{} should pass selection criteria",
                tno.name
            );
        }
    }

    #[test]
    fn test_generate_clones_count() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let tno = &distant_tno_sample()[0];
        let clones = generate_clones(tno, 25, &mut rng);
        assert_eq!(clones.len(), 25);
    }

    #[test]
    fn test_clones_near_parent() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);
        let tno = &distant_tno_sample()[0]; // Sedna
        let clones = generate_clones(tno, 25, &mut rng);

        for c in &clones {
            assert!(
                (c.a - tno.a).abs() < 5.0 * tno.sigma_a,
                "clone a = {} too far from parent a = {}",
                c.a,
                tno.a
            );
        }
    }

    #[test]
    fn test_clone_angles_and_anomaly_perturbed() {
        // M9: all elements vary across clones, including angles and M.
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);
        let tno = &distant_tno_sample()[0];
        let clones = generate_clones(tno, 50, &mut rng);

        let distinct = |get: fn(&TnoClone) -> f64| {
            let v0 = get(&clones[0]);
            clones.iter().any(|c| (get(c) - v0).abs() > 1e-9)
        };
        assert!(distinct(|c| c.omega), "ω never perturbed");
        assert!(distinct(|c| c.omega_big), "Ω never perturbed");
        assert!(distinct(|c| c.i), "i never perturbed");
        assert!(distinct(|c| c.mean_anomaly), "M never drawn");

        // Mean anomalies roughly cover the circle (uniform draw).
        let spread = clones
            .iter()
            .map(|c| c.mean_anomaly)
            .fold((f64::MAX, f64::MIN), |(lo, hi), m| (lo.min(m), hi.max(m)));
        assert!(spread.1 - spread.0 > 3.0, "M spread = {:?}", spread);
    }

    #[test]
    fn test_correlated_a_e_preserves_perihelion() {
        // The (a, e) draws move along q ≈ const: clone perihelia scatter
        // by ~σ_e·a (≈ 2.5 AU for Sedna), much tighter than the naive
        // independent-draw scatter σ_a (≈ 10 AU).
        let mut rng = rand::rngs::StdRng::seed_from_u64(11);
        let tno = &distant_tno_sample()[0]; // Sedna, q = 75.9 AU
        let clones = generate_clones(tno, 200, &mut rng);
        let q0 = tno.perihelion_typed().get::<astronomical_unit>();
        for c in &clones {
            let q = c.perihelion_typed().get::<astronomical_unit>();
            assert!(
                (q - q0).abs() < 5.0 * SIGMA_E * c.a,
                "{}: clone q = {:.1} vs parent q = {:.1}",
                c.parent_name,
                q,
                q0
            );
        }
    }

    #[test]
    fn test_clone_eccentricity_valid_and_no_silent_crossers() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(99);
        for tno in distant_tno_sample() {
            let clones = generate_clones(&tno, 100, &mut rng);
            for c in &clones {
                assert!(c.e > 0.0 && c.e < 1.0, "e = {} invalid", c.e);
                // Parent is detached (q > a_N): clones either stay
                // outside Neptune or carry the explicit flag.
                if !c.neptune_crossing {
                    let q = c.perihelion_typed().get::<astronomical_unit>();
                    assert!(
                        q > A_NEPTUNE_AU,
                        "{}: silent Neptune crosser q = {:.1}",
                        c.parent_name,
                        q
                    );
                }
            }
        }
    }

    #[test]
    fn test_clone_elements_roundtrip() {
        let mut rng = rand::rngs::StdRng::seed_from_u64(3);
        let tno = &distant_tno_sample()[1];
        let clone = &generate_clones(tno, 1, &mut rng)[0];
        let elem = clone.elements();
        assert!((elem.a - clone.a).abs() < 1e-12);
        assert!(
            (elem.perihelion_typed().get::<astronomical_unit>()
                - clone.perihelion_typed().get::<astronomical_unit>())
            .abs()
                < 1e-9
        );
    }
}
