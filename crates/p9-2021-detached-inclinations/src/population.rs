//! Seeded synthetic population of detached Kuiper belt objects and the
//! inclination statistics that distinguish a Planet-Nine Solar System from a
//! Planet-Nine-free one.
//!
//! The headline of Anderson & Kaib (2021): an inclined Planet Nine secularly
//! excites and **broadens** the inclination distribution of the detached
//! Kuiper belt. We realize this as follows.
//!
//! Each particle is born in (or near) the giant-planet invariable plane with a
//! small *free* inclination `i_free` drawn from a Rayleigh distribution (the
//! cold detached disk). Its orbital pole then precesses on a cone of opening
//! angle `i_free` about the **forced pole** — which, with Planet Nine, is
//! tilted away from the invariable plane by the linear forced inclination
//! `i_forced(a, e, i_free, P9)` (see [`crate::forcing`]). Sampling a uniform
//! secular precession phase `φ` for each particle, its instantaneous
//! inclination to the invariable plane is the spherical-triangle combination
//!
//!   cos i = cos i_forced cos i_free + sin i_forced sin i_free cos φ.
//!
//! With **no** Planet Nine, i_forced = 0 and the observed inclination is just
//! the free inclination — a tight, low-i distribution. **With** an inclined
//! Planet Nine, the forced tilt offsets every pole, and because i_forced grows
//! with semi-major axis the population fans out: the dispersion broadens and a
//! tail of moderately/highly inclined detached objects appears that the
//! 8-planet model does not produce. That broadening is the testable signature.
//!
//! The population is built from a fixed `StdRng` seed so the statistics are
//! reproducible and can be pinned in tests.

use p9_core::constants::DEG2RAD;
use p9_core::types::P9Params;
use p9_core::units::{au, radians, Angle, Length};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_distr::{Distribution, Uniform};

use crate::forcing::forced_inclination;

/// A single detached test particle.
#[derive(Debug, Clone, Copy)]
pub struct DetachedParticle {
    /// Semi-major axis (AU).
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
    /// Free (intrinsic) inclination to the invariable plane (rad).
    pub i_free: f64,
    /// Secular precession phase of the node about the forced pole (rad).
    pub phase: f64,
}

impl DetachedParticle {
    /// Perihelion distance q = a(1 − e) (AU).
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }

    /// Observed inclination to the invariable plane (rad) given a forced
    /// inclination `i_forced`: the spherical combination of the free
    /// inclination on a cone about the forced pole at this precession phase.
    pub fn observed_inclination(&self, i_forced: f64) -> f64 {
        let cos_i = i_forced.cos() * self.i_free.cos()
            + i_forced.sin() * self.i_free.sin() * self.phase.cos();
        cos_i.clamp(-1.0, 1.0).acos()
    }

    /// Semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a)
    }

    /// Perihelion distance as a typed [`Length`].
    pub fn perihelion_typed(&self) -> Length {
        au(self.perihelion())
    }

    /// Observed inclination as a typed [`Angle`] (dimension-checked view of
    /// [`observed_inclination`](Self::observed_inclination)).
    pub fn observed_inclination_typed(&self, i_forced: f64) -> Angle {
        radians(self.observed_inclination(i_forced))
    }
}

/// Parameters of the detached-disk realization.
#[derive(Debug, Clone, Copy)]
pub struct PopulationConfig {
    pub n: usize,
    /// Minimum semi-major axis of the detached region (AU).
    pub a_min: f64,
    /// Maximum semi-major axis (AU).
    pub a_max: f64,
    /// Minimum perihelion to count as "detached" (decoupled from Neptune, AU).
    pub q_min: f64,
    /// Rayleigh scale of the free-inclination distribution (rad). The cold
    /// detached disk: σ ≈ 10–15°.
    pub i_free_sigma: f64,
    /// RNG seed.
    pub seed: u64,
}

impl Default for PopulationConfig {
    fn default() -> Self {
        Self {
            n: 4000,
            a_min: 250.0,
            a_max: 500.0,
            q_min: 40.0,
            i_free_sigma: 12.0 * DEG2RAD,
            seed: 0xA21_DE7A,
        }
    }
}

/// Sample a free inclination from a Rayleigh-in-inclination distribution
/// f(i) ∝ i exp(−i²/2σ²) via inverse transform from a uniform deviate.
fn sample_free_inclination(u: f64, sigma: f64) -> f64 {
    // i = σ √(−2 ln(1 − u)); guard u → 1.
    sigma * (-2.0 * (1.0 - u).max(1e-12).ln()).sqrt()
}

/// Build a seeded detached population. Semi-major axes are uniform in `a`,
/// perihelia uniform in q ∈ [q_min, a) (so e follows), free inclinations
/// Rayleigh-distributed, and secular phases uniform on [0, 2π).
pub fn build_population(cfg: &PopulationConfig) -> Vec<DetachedParticle> {
    let mut rng = StdRng::seed_from_u64(cfg.seed);
    let ua = Uniform::new(cfg.a_min, cfg.a_max);
    let uphase = Uniform::new(0.0, std::f64::consts::TAU);

    let mut out = Vec::with_capacity(cfg.n);
    while out.len() < cfg.n {
        let a = ua.sample(&mut rng);
        // Detached: q between q_min and (just below) a; draw q, derive e.
        let q_max = (a - 1.0).max(cfg.q_min + 1.0);
        let q = rng.gen_range(cfg.q_min..q_max);
        let e = 1.0 - q / a;
        if !(0.0..1.0).contains(&e) {
            continue;
        }
        let i_free = sample_free_inclination(rng.gen::<f64>(), cfg.i_free_sigma);
        let phase = uphase.sample(&mut rng);
        out.push(DetachedParticle {
            a,
            e,
            i_free,
            phase,
        });
    }
    out
}

/// Observed inclinations (rad) of the whole population under a given Planet
/// Nine (its inclination forcing applied per particle).
pub fn observed_inclinations(pop: &[DetachedParticle], p9: &P9Params) -> Vec<f64> {
    pop.iter()
        .map(|p| {
            let i_forced = forced_inclination(p.a, p.e, p.i_free, p9);
            p.observed_inclination(i_forced)
        })
        .collect()
}

/// Observed inclinations (rad) with **no** Planet Nine: the forced inclination
/// vanishes, so the observed inclination is the free inclination.
pub fn observed_inclinations_no_p9(pop: &[DetachedParticle]) -> Vec<f64> {
    pop.iter().map(|p| p.i_free).collect()
}

/// Population mean of a slice (rad).
pub fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().sum::<f64>() / values.len() as f64
}

/// Population standard deviation (the inclination *dispersion*, rad).
pub fn dispersion(values: &[f64]) -> f64 {
    if values.len() < 2 {
        return 0.0;
    }
    let m = mean(values);
    let var = values.iter().map(|v| (v - m).powi(2)).sum::<f64>() / values.len() as f64;
    var.sqrt()
}

/// Fraction of the population with inclination above a threshold (rad) — the
/// moderately/highly inclined tail that Planet Nine populates.
pub fn fraction_above(values: &[f64], threshold: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.iter().filter(|&&v| v > threshold).count() as f64 / values.len() as f64
}

/// Mean forced inclination (rad) over the population for a given Planet Nine —
/// a compact summary of how strongly the disk is being tilted.
pub fn mean_forced_inclination(pop: &[DetachedParticle], p9: &P9Params) -> f64 {
    let forced: Vec<f64> = pop
        .iter()
        .map(|p| forced_inclination(p.a, p.e, p.i_free, p9))
        .collect();
    mean(&forced)
}

#[cfg(test)]
mod tests {
    use super::*;
    use p9_core::constants::RAD2DEG;

    fn nominal() -> P9Params {
        P9Params::nominal_2016()
    }

    #[test]
    fn typed_particle_accessors_match_f64() {
        use approx::assert_relative_eq;
        use uom::si::angle::radian;
        use uom::si::length::astronomical_unit;
        let p = DetachedParticle {
            a: 400.0,
            e: 0.7,
            i_free: 0.2,
            phase: 0.5,
        };
        assert_relative_eq!(
            p.semi_major_axis().get::<astronomical_unit>(),
            p.a,
            epsilon = 1e-12
        );
        assert_relative_eq!(
            p.perihelion_typed().get::<astronomical_unit>(),
            p.perihelion(),
            epsilon = 1e-12
        );
        assert_relative_eq!(
            p.observed_inclination_typed(0.3).get::<radian>(),
            p.observed_inclination(0.3),
            epsilon = 1e-12
        );
    }

    #[test]
    fn test_population_is_detached_and_reproducible() {
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        assert_eq!(pop.len(), cfg.n);
        for p in &pop {
            assert!(p.perihelion() >= cfg.q_min - 1e-9, "q = {}", p.perihelion());
            assert!(p.a >= cfg.a_min && p.a <= cfg.a_max);
            assert!((0.0..1.0).contains(&p.e));
        }
        // Same seed → identical population.
        let pop2 = build_population(&cfg);
        assert_eq!(pop[0].a, pop2[0].a);
        assert_eq!(pop[cfg.n - 1].i_free, pop2[cfg.n - 1].i_free);
    }

    #[test]
    fn test_dispersion_larger_with_p9() {
        // THE HEADLINE: an inclined Planet Nine broadens the inclination
        // dispersion of the detached disk relative to a P9-free Solar System.
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let with = observed_inclinations(&pop, &nominal());
        let without = observed_inclinations_no_p9(&pop);
        let sd_with = dispersion(&with);
        let sd_without = dispersion(&without);
        assert!(
            sd_with > sd_without,
            "dispersion with P9 ({:.2}°) must exceed without ({:.2}°)",
            sd_with * RAD2DEG,
            sd_without * RAD2DEG
        );
    }

    #[test]
    fn test_high_inclination_tail_populated_by_p9() {
        // More moderately/highly inclined detached objects exist with P9.
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let with = observed_inclinations(&pop, &nominal());
        let without = observed_inclinations_no_p9(&pop);
        let thresh = 30.0 * DEG2RAD;
        let f_with = fraction_above(&with, thresh);
        let f_without = fraction_above(&without, thresh);
        assert!(
            f_with > f_without,
            "i>30° fraction with P9 ({:.3}) must exceed without ({:.3})",
            f_with,
            f_without
        );
    }

    #[test]
    fn test_mean_inclination_raised_by_p9() {
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let with = mean(&observed_inclinations(&pop, &nominal()));
        let without = mean(&observed_inclinations_no_p9(&pop));
        assert!(
            with > without,
            "mean i with P9 ({:.2}°) must exceed without ({:.2}°)",
            with * RAD2DEG,
            without * RAD2DEG
        );
    }

    #[test]
    fn test_mean_forced_inclination_grows_with_mass_and_i9() {
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let light = P9Params {
            mass_earth: 5.0,
            i: 15.0 * DEG2RAD,
            ..nominal()
        };
        let heavy = P9Params {
            mass_earth: 15.0,
            i: 30.0 * DEG2RAD,
            ..nominal()
        };
        let fl = mean_forced_inclination(&pop, &light);
        let fh = mean_forced_inclination(&pop, &heavy);
        assert!(fh > fl, "mean forced i: light {fl} vs heavy {fh}");
    }

    #[test]
    fn test_broadening_present_at_every_p9_mass() {
        // Any inclined Planet Nine broadens the disk relative to the P9-free
        // baseline. (Strict monotonicity in mass holds only in the weak-forcing
        // regime: once i_forced approaches the free-inclination scale the
        // spherical cone saturates and the dispersion plateaus — see
        // REPRODUCTION_NOTES; we pin the robust, paper-relevant claim instead.)
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let baseline = dispersion(&observed_inclinations_no_p9(&pop));
        for &m in &[2.0_f64, 5.0, 10.0, 20.0] {
            let p9 = P9Params {
                mass_earth: m,
                ..nominal()
            };
            let sd = dispersion(&observed_inclinations(&pop, &p9));
            assert!(
                sd > baseline,
                "dispersion at m = {m} ({sd}) must exceed P9-free baseline ({baseline})"
            );
        }
    }

    #[test]
    fn test_dispersion_grows_with_mass_in_weak_forcing_regime() {
        // In the weak-forcing regime (a more distant Planet Nine, so the
        // forced inclination stays below the free-inclination scale) the
        // broadening is monotone in mass, as linear theory predicts.
        let cfg = PopulationConfig::default();
        let pop = build_population(&cfg);
        let distant = P9Params {
            a: 1500.0,
            ..nominal()
        };
        let masses = [0.0_f64, 3.0, 6.0, 12.0];
        let mut prev = -1.0;
        for &m in &masses {
            let sd = if m == 0.0 {
                dispersion(&observed_inclinations_no_p9(&pop))
            } else {
                let p9 = P9Params {
                    mass_earth: m,
                    ..distant
                };
                dispersion(&observed_inclinations(&pop, &p9))
            };
            assert!(
                sd > prev,
                "dispersion not monotone at m = {m}: {sd} <= {prev}"
            );
            prev = sd;
        }
    }
}
