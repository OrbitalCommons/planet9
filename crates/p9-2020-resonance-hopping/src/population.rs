//! Seeded synthetic distant-TNO population and per-class fraction analysis.
//!
//! We draw a population of distant TNOs in semimajor axis (the resonance
//! landscape lives in `a`), assign each a perihelion / eccentricity in the
//! detached-to-scattering range, classify each against the P9 resonance
//! landscape, and report the fraction in each dynamical class — globally and
//! per semimajor-axis bin, to expose the rise of the hopping fraction with a.

use crate::classification::{classify, Class, ResonanceLandscape};
use crate::published;
use p9_core::constants::EARTH_MASS_SOLAR;
use rand::rngs::StdRng;
use rand::SeedableRng;
use rand_distr::{Distribution, Uniform};
use serde::{Deserialize, Serialize};

/// A synthetic distant TNO: just the orbital quantities the resonance
/// classification needs.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Tno {
    /// Semimajor axis (AU).
    pub a: f64,
    /// Eccentricity.
    pub e: f64,
}

impl Tno {
    /// Perihelion distance q = a(1 − e) (AU).
    pub fn perihelion(&self) -> f64 {
        self.a * (1.0 - self.e)
    }
}

/// Draw a seeded synthetic distant-TNO population.
///
/// `a` is drawn uniformly in `[a_min, a_max]` (AU), spanning the band interior
/// to Planet Nine where its n:1 / n:2 resonances live. Eccentricity is drawn
/// so the perihelion lands in the distant detached/scattering range (q ≳ 30
/// AU): we draw q uniformly in `[q_min, q_max]` and set e = 1 − q/a. This
/// keeps the population physically distant (decoupled from Neptune) and gives
/// the eccentricities that set the resonance widths.
pub fn synthetic_population(
    seed: u64,
    n: usize,
    a_min: f64,
    a_max: f64,
    q_min: f64,
    q_max: f64,
) -> Vec<Tno> {
    let mut rng = StdRng::seed_from_u64(seed);
    let a_dist = Uniform::new(a_min, a_max);
    let q_dist = Uniform::new(q_min, q_max);
    (0..n)
        .map(|_| {
            let a = a_dist.sample(&mut rng);
            let q = q_dist.sample(&mut rng).min(a * 0.999);
            let e = (1.0 - q / a).clamp(0.0, 0.99);
            Tno { a, e }
        })
        .collect()
}

/// Fractions of a population in each dynamical class. The three fields are
/// valid probabilities and sum to 1.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ClassFractions {
    pub resonant: f64,
    pub hopping: f64,
    pub non_resonant: f64,
    /// Number of objects the fractions were computed over.
    pub n: usize,
}

impl ClassFractions {
    /// Sum of the three fractions (must be 1 for a non-empty population).
    pub fn total(&self) -> f64 {
        self.resonant + self.hopping + self.non_resonant
    }
}

/// Classify a population against `landscape` and tabulate class fractions.
///
/// `m_p9_solar` is Planet Nine's mass (solar masses); `width_coeff` the
/// labelled pendulum-width constant passed through to the classifier.
pub fn class_fractions(
    population: &[Tno],
    landscape: &ResonanceLandscape,
    m_p9_solar: f64,
    width_coeff: f64,
) -> ClassFractions {
    let n = population.len();
    if n == 0 {
        return ClassFractions {
            resonant: 0.0,
            hopping: 0.0,
            non_resonant: 0.0,
            n: 0,
        };
    }
    let (mut res, mut hop, mut non) = (0usize, 0usize, 0usize);
    for tno in population {
        match classify(landscape, tno.a, tno.e, m_p9_solar, width_coeff).class {
            Class::Resonant => res += 1,
            Class::Hopping => hop += 1,
            Class::NonResonant => non += 1,
        }
    }
    let nf = n as f64;
    ClassFractions {
        resonant: res as f64 / nf,
        hopping: hop as f64 / nf,
        non_resonant: non as f64 / nf,
        n,
    }
}

/// Hopping fraction per semimajor-axis bin: split `[a_min, a_max]` into
/// `n_bins` equal-width bins and return the (bin-centre, hopping fraction)
/// pairs for the bins that contain objects. Used to show the hopping fraction
/// rising with a.
pub fn hopping_fraction_by_bin(
    population: &[Tno],
    landscape: &ResonanceLandscape,
    m_p9_solar: f64,
    width_coeff: f64,
    a_min: f64,
    a_max: f64,
    n_bins: usize,
) -> Vec<(f64, f64)> {
    let width = (a_max - a_min) / n_bins as f64;
    let mut counts = vec![(0usize, 0usize); n_bins]; // (hopping, total)
    for tno in population {
        if tno.a < a_min || tno.a >= a_max {
            continue;
        }
        let bin = (((tno.a - a_min) / width) as usize).min(n_bins - 1);
        let is_hop = matches!(
            classify(landscape, tno.a, tno.e, m_p9_solar, width_coeff).class,
            Class::Hopping
        );
        counts[bin].1 += 1;
        if is_hop {
            counts[bin].0 += 1;
        }
    }
    counts
        .into_iter()
        .enumerate()
        .filter(|(_, (_, total))| *total > 0)
        .map(|(i, (hop, total))| {
            let centre = a_min + (i as f64 + 0.5) * width;
            (centre, hop as f64 / total as f64)
        })
        .collect()
}

/// Default Planet Nine mass for the analysis (5 M⊕, in the Khain-era nominal
/// range). Exposed so callers/tests can vary it.
pub fn default_m_p9_solar() -> f64 {
    5.0 * EARTH_MASS_SOLAR
}

/// Resonance-spectrum depth: the maximum P9-orbit coefficient q resolved. The
/// dense accumulation of high-order resonances near a₉ needs a deep spectrum
/// for the overlap K to cross 1; q ≤ 20 captures the period-ratio-→-1 crowding
/// that drives the hopping zone in the outer band.
pub const SPECTRUM_Q_MAX: u32 = 20;

/// Convenience: build the canonical landscape for the assumed P9 over the
/// distant band `[a_min, a_max]`.
pub fn default_landscape(a_min: f64, a_max: f64) -> ResonanceLandscape {
    ResonanceLandscape::new(published::A_P9_AU, a_min, a_max, SPECTRUM_Q_MAX)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::classification::resonance_overlap_k;

    const A_MIN: f64 = 150.0;
    const A_MAX: f64 = 499.0;
    // Pendulum-width constant: a labelled O(1) model constant (not tuned to a
    // published fraction). Sets the absolute resonance widths; the class
    // ORDERING and the rise of hopping with a are independent of it.
    const WIDTH: f64 = 1.2;

    fn setup(seed: u64) -> (Vec<Tno>, ResonanceLandscape, f64) {
        let pop = synthetic_population(seed, 4000, A_MIN, A_MAX, 30.0, 120.0);
        let landscape = default_landscape(A_MIN, A_MAX);
        (pop, landscape, default_m_p9_solar())
    }

    #[test]
    fn fractions_are_valid_probabilities_summing_to_one() {
        let (pop, l, m) = setup(2020);
        let f = class_fractions(&pop, &l, m, WIDTH);
        for x in [f.resonant, f.hopping, f.non_resonant] {
            assert!((0.0..=1.0).contains(&x), "fraction {x} out of [0,1]");
        }
        assert!(
            (f.total() - 1.0).abs() < 1e-12,
            "fractions must sum to 1, got {}",
            f.total()
        );
        assert_eq!(f.n, pop.len());
    }

    #[test]
    fn all_three_classes_are_populated() {
        // The reduced model must actually produce a mix, not collapse to one
        // class — resonant sticking, hopping, and non-resonant all present.
        let (pop, l, m) = setup(2020);
        let f = class_fractions(&pop, &l, m, WIDTH);
        assert!(f.resonant > 0.0, "expected some resonant objects");
        assert!(f.hopping > 0.0, "expected some hopping objects");
        assert!(f.non_resonant > 0.0, "expected some non-resonant objects");
    }

    #[test]
    fn hopping_fraction_increases_with_semimajor_axis() {
        // Headline: the hopping fraction rises with the semimajor-axis bin.
        // We assert a monotone NON-DECREASING trend across coarse bins (the
        // mechanism — collapsing n:1 spacing toward a₉ — is monotone; finite
        // sampling can tie adjacent bins but must not reverse materially).
        let (pop, l, m) = setup(2020);
        let bins = hopping_fraction_by_bin(&pop, &l, m, WIDTH, A_MIN, A_MAX, 5);
        assert!(bins.len() >= 3, "need several populated bins");
        // First populated bin must be strictly below the last.
        let first = bins.first().unwrap().1;
        let last = bins.last().unwrap().1;
        assert!(
            last > first,
            "hopping fraction should rise with a: {first:.3} (a≈{:.0}) -> {last:.3} (a≈{:.0})",
            bins.first().unwrap().0,
            bins.last().unwrap().0,
        );
        // And it should be monotone non-decreasing bin-to-bin within noise.
        for w in bins.windows(2) {
            assert!(
                w[1].1 + 0.06 >= w[0].1,
                "hopping fraction reversed: {:.3} (a≈{:.0}) -> {:.3} (a≈{:.0})",
                w[0].1,
                w[0].0,
                w[1].1,
                w[1].0,
            );
        }
    }

    #[test]
    fn resonant_sticking_dominates_at_small_a() {
        // Clean resonant sticking should dominate over hopping in the inner
        // third (well-separated low-order n:1 resonances), and hopping should
        // dominate resonant in the outer third (overlapping resonances). The
        // middle third is the transition where K crosses 1.
        let (pop, l, m) = setup(2020);
        let third = (A_MAX - A_MIN) / 3.0;
        let lo_hi = A_MIN + third;
        let hi_lo = A_MAX - third;
        let inner: Vec<Tno> = pop.iter().copied().filter(|t| t.a < lo_hi).collect();
        let outer: Vec<Tno> = pop.iter().copied().filter(|t| t.a >= hi_lo).collect();
        let fi = class_fractions(&inner, &l, m, WIDTH);
        let fo = class_fractions(&outer, &l, m, WIDTH);
        assert!(
            fi.resonant > fi.hopping,
            "inner third (a<{lo_hi:.0}): resonant {:.3} should beat hopping {:.3}",
            fi.resonant,
            fi.hopping
        );
        assert!(
            fo.hopping > fo.resonant,
            "outer third (a>={hi_lo:.0}): hopping {:.3} should beat resonant {:.3}",
            fo.hopping,
            fo.resonant
        );
    }

    #[test]
    fn results_are_seed_robust() {
        // The class fractions must be stable across seeds (population-level
        // statistics, not a single draw).
        let l = default_landscape(A_MIN, A_MAX);
        let m = default_m_p9_solar();
        let mut fr = Vec::new();
        for seed in [1u64, 7, 42, 2020, 777] {
            let pop = synthetic_population(seed, 4000, A_MIN, A_MAX, 30.0, 120.0);
            fr.push(class_fractions(&pop, &l, m, WIDTH));
        }
        let mean_hop = fr.iter().map(|f| f.hopping).sum::<f64>() / fr.len() as f64;
        for f in &fr {
            assert!(
                (f.hopping - mean_hop).abs() < 0.03,
                "hopping fraction {:.3} deviates from mean {:.3} by > 0.03",
                f.hopping,
                mean_hop
            );
            assert!((f.total() - 1.0).abs() < 1e-12);
        }
    }

    #[test]
    fn overlap_k_crosses_unity_within_the_band() {
        // Sanity: K must actually cross 1 somewhere in the band, otherwise the
        // resonant/hopping split is trivial. Low a ⇒ K < 1, high a ⇒ K > 1.
        let l = default_landscape(A_MIN, A_MAX);
        let m = default_m_p9_solar();
        let k_low = resonance_overlap_k(&l, 180.0, 0.6, m, WIDTH).unwrap();
        let k_high = resonance_overlap_k(&l, 480.0, 0.6, m, WIDTH).unwrap();
        assert!(
            k_low < 1.0,
            "expected isolated resonances at small a, K={k_low:.3}"
        );
        assert!(k_high > 1.0, "expected overlap at large a, K={k_high:.3}");
    }
}
