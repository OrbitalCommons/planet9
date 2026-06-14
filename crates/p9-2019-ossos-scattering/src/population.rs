//! Synthetic actively-scattering TNO population and the P9 detachment test.
//!
//! Kaib et al. (2019) constrain P9 with the *observed* scattering population:
//! large-a objects with perihelia near Neptune (q ≈ 30–38 AU) undergoing a
//! semimajor-axis random walk. The high-q edge of that population — the q at
//! which objects stop scattering and detach — is the diagnostic. A massive
//! distant P9 lifts perihelia secularly (`planet_nine`), pushing a fraction of
//! the scattering objects above the resonance-overlap divide (`boundary`) and
//! into the detached class, raising the detached fraction.
//!
//! Here we draw a seeded synthetic scattering population, evolve each object's
//! perihelion under the secular P9 lift, and measure the detached fraction
//! with and without P9. The increase is computed, non-zero, and reproducible.

use crate::boundary::{classify, DynamicalClass};
use crate::planet_nine::{perihelion_lift, PlanetNine};
use p9_core::constants::A_NEPTUNE_AU;
use rand::SeedableRng;
use rand_distr::{Distribution, Uniform};

/// One synthetic scattering TNO.
#[derive(Debug, Clone, Copy)]
pub struct ScatteringTno {
    pub a_au: f64,
    pub q_au: f64,
}

/// Draw `n` synthetic actively-scattering TNOs (seeded).
///
/// Semi-major axes are log-uniform over [`a_min`, `a_max`] (the survey detects
/// objects across decades in a); perihelia are uniform over the Neptune-coupled
/// band q ∈ [q_lo, q_hi] with q_lo just outside Neptune's orbit. We then keep
/// only objects that are actually scattering at draw time (q < q_crit(a)), so
/// the sample is, by construction, the actively-scattering population the
/// paper studies.
pub fn synthetic_scattering_population(
    seed: u64,
    n: usize,
    a_min: f64,
    a_max: f64,
    q_lo: f64,
    q_hi: f64,
) -> Vec<ScatteringTno> {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let log_a = Uniform::new(a_min.ln(), a_max.ln());
    let q_dist = Uniform::new(q_lo, q_hi);

    let mut out = Vec::with_capacity(n);
    // Oversample then filter to the scattering class to reach n members.
    while out.len() < n {
        let a = log_a.sample(&mut rng).exp();
        let q = q_dist.sample(&mut rng);
        if q <= A_NEPTUNE_AU {
            continue; // inside Neptune: not a scattering-disk object
        }
        if classify(a, q) == DynamicalClass::Scattering {
            out.push(ScatteringTno { a_au: a, q_au: q });
        }
    }
    out
}

/// Result of the detachment comparison.
#[derive(Debug, Clone, Copy)]
pub struct DetachmentResult {
    /// Detached fraction with no Planet Nine (control).
    pub detached_fraction_no_p9: f64,
    /// Detached fraction after applying the P9 secular perihelion lift.
    pub detached_fraction_with_p9: f64,
    /// Number of objects in the population.
    pub n: usize,
}

impl DetachmentResult {
    /// The increase in detached fraction attributable to P9.
    pub fn delta(&self) -> f64 {
        self.detached_fraction_with_p9 - self.detached_fraction_no_p9
    }
}

/// For each object, raise perihelion by the secular P9 lift accumulated over
/// `time_days` and re-classify. Returns the detached fraction with and without
/// P9.
pub fn detached_fraction_comparison(
    p9: &PlanetNine,
    population: &[ScatteringTno],
    time_days: f64,
) -> DetachmentResult {
    let n = population.len();
    assert!(n > 0);

    // Without P9 the population is, by construction, fully scattering: 0
    // detached. We still compute it so the comparison is symmetric and would
    // catch any boundary inconsistency.
    let detached_no = population
        .iter()
        .filter(|t| classify(t.a_au, t.q_au) == DynamicalClass::Detached)
        .count();

    let detached_with = population
        .iter()
        .filter(|t| {
            let dq = perihelion_lift(p9, t.a_au, t.q_au, time_days);
            classify(t.a_au, t.q_au + dq) == DynamicalClass::Detached
        })
        .count();

    DetachmentResult {
        detached_fraction_no_p9: detached_no as f64 / n as f64,
        detached_fraction_with_p9: detached_with as f64 / n as f64,
        n,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::planet_nine::SECULAR_SCULPTING_DAYS;

    fn standard_population() -> Vec<ScatteringTno> {
        // Large-a scattering disk: a 150-900 AU, q in the Neptune-coupled band.
        synthetic_scattering_population(1905, 300, 150.0, 900.0, 31.0, 40.0)
    }

    #[test]
    fn synthetic_population_is_all_scattering() {
        let pop = standard_population();
        assert_eq!(pop.len(), 300);
        for t in &pop {
            assert_eq!(classify(t.a_au, t.q_au), DynamicalClass::Scattering);
            assert!(t.q_au > A_NEPTUNE_AU);
        }
    }

    #[test]
    fn population_is_reproducible() {
        let a = standard_population();
        let b = standard_population();
        assert_eq!(a.len(), b.len());
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.a_au, y.a_au);
            assert_eq!(x.q_au, y.q_au);
        }
    }

    #[test]
    fn p9_raises_detached_fraction_by_nonzero_amount() {
        let pop = standard_population();
        let p9 = PlanetNine::default();
        let res = detached_fraction_comparison(&p9, &pop, SECULAR_SCULPTING_DAYS);
        // Control population is fully scattering.
        assert_eq!(res.detached_fraction_no_p9, 0.0);
        // P9 lifts perihelia -> a non-zero fraction detaches.
        assert!(
            res.delta() > 0.0,
            "P9 should raise detached fraction; delta = {}",
            res.delta()
        );
        assert!(res.detached_fraction_with_p9 <= 1.0);
        // OSSOS XV: a nominal P9 over-detaches the large-a scattering
        // population, the basis for disfavoring such configurations.
        assert!(
            res.detached_fraction_with_p9 > 0.5,
            "nominal P9 should strongly deplete the scattering population: {}",
            res.detached_fraction_with_p9
        );
    }

    #[test]
    fn more_massive_p9_detaches_at_least_as_many() {
        let pop = standard_population();
        let light = PlanetNine {
            mass_earth: 2.0,
            ..Default::default()
        };
        let heavy = PlanetNine {
            mass_earth: 6.0,
            ..Default::default()
        };
        let dl = detached_fraction_comparison(&light, &pop, SECULAR_SCULPTING_DAYS)
            .detached_fraction_with_p9;
        let dh = detached_fraction_comparison(&heavy, &pop, SECULAR_SCULPTING_DAYS)
            .detached_fraction_with_p9;
        assert!(
            dh >= dl,
            "heavier P9 should detach >= as many: light {dl} heavy {dh}"
        );
        assert!(dl > 0.0, "even a light P9 should detach some: {dl}");
    }
}
