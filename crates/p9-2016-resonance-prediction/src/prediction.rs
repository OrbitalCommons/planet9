//! Period-ratio resonance localization of Planet Nine.
//!
//! Following Malhotra, Volk & Wang (2016): Planet Nine is the *outer*,
//! longest-period body, and the distant ETNOs are interior objects each
//! trapped in a small-integer mean-motion resonance with it. An ETNO that
//! completes q orbits while Planet Nine completes p (p > q) sits in the
//! interior p:q resonance at
//!
//!   a_ETNO = a9 · (q/p)^(2/3)   ⇔   P9 / P_ETNO = p / q,
//!
//! so the observed ETNO periods are commensurate with — and predict — P9's
//! period. We restrict to N:1 and N:2 resonances (Malhotra et al.'s case),
//! scan a9, and pick the value that places the distant set closest to those
//! low-order ratios.

use p9_core::data::etno::{Etno, BROWN_2017_SAMPLE};

/// Malhotra, Volk & Wang (2016) best-fit Planet Nine semi-major axis (AU).
pub const MALHOTRA_2016_A9_AU: f64 = 665.0;

/// Malhotra, Volk & Wang (2016) corresponding Planet Nine period (years):
/// P = a9^(3/2) ≈ 665^1.5 ≈ 17,117 yr.
pub const MALHOTRA_2016_PERIOD_YR: f64 = 17_117.0;

/// Millholland & Laughlin (2016) best-fit Planet Nine semi-major axis (AU).
pub const MILLHOLLAND_2016_A9_AU: f64 = 654.0;

/// Default highest order on the *p* side of the scanned p:q resonances.
pub const DEFAULT_P_MAX: u32 = 9;

/// Default highest order on the *q* side. Malhotra et al. consider N:1 and
/// N:2 resonances, i.e. q ∈ {1, 2}.
pub const DEFAULT_Q_MAX: u32 = 2;

/// Heliocentric orbital period in years from the semi-major axis (AU), via
/// Kepler's third law for a one-solar-mass primary: P[yr] = a[AU]^(3/2).
///
/// The neglected ~1e-4 mass of Planet Nine and the ETNOs is far below the
/// precision of these predictions, so the bare Kepler relation is used (it is
/// exactly the relation the resonance arithmetic inverts).
pub fn orbital_period_years(a_au: f64) -> f64 {
    a_au.powf(1.5)
}

/// Period ratio of the outer body to the inner body, P_outer / P_inner =
/// (a_outer / a_inner)^(3/2). With Planet Nine outermost this is P9 / P_ETNO.
pub fn period_ratio(a_outer: f64, a_inner: f64) -> f64 {
    orbital_period_years(a_outer) / orbital_period_years(a_inner)
}

/// A candidate small-integer interior resonance p:q (the ETNO completes q
/// orbits per p orbits of Planet Nine) and how well a given ETNO matches it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ResonanceMatch {
    /// Number of Planet Nine orbits per resonance cycle (the larger integer).
    pub p: u32,
    /// Number of ETNO orbits per resonance cycle (the smaller integer).
    pub q: u32,
    /// The observed period ratio P9 / P_ETNO for the test a9.
    pub period_ratio: f64,
    /// Residual |observed ratio − p/q|.
    pub residual: f64,
}

/// For an ETNO at `a_etno` interior to a candidate Planet Nine at `a9`, find
/// the small-integer resonance p:q (p ≤ `p_max`, q ≤ `q_max`, p > q,
/// gcd(p,q)=1) whose ratio p/q is closest to the observed period ratio
/// P9 / P_ETNO. Returns `None` if `a9 <= a_etno` (Planet Nine must be the
/// outer body for the ETNO to be in an interior resonance with it).
pub fn nearest_resonance(a_etno: f64, a9: f64, p_max: u32, q_max: u32) -> Option<ResonanceMatch> {
    if a9 <= a_etno {
        return None;
    }
    let ratio = period_ratio(a9, a_etno);
    let mut best: Option<ResonanceMatch> = None;
    for q in 1..=q_max {
        for p in (q + 1)..=p_max {
            if gcd(p, q) != 1 {
                continue;
            }
            let residual = (ratio - p as f64 / q as f64).abs();
            if best.is_none_or(|b| residual < b.residual) {
                best = Some(ResonanceMatch {
                    p,
                    q,
                    period_ratio: ratio,
                    residual,
                });
            }
        }
    }
    best
}

/// Residual of a single ETNO against its nearest low-order resonance with a
/// candidate Planet Nine at `a9` (0 if perfectly commensurate). ETNOs exterior
/// to a9 contribute no constraint and return `None`.
pub fn resonance_residual(a_etno: f64, a9: f64, p_max: u32, q_max: u32) -> Option<f64> {
    nearest_resonance(a_etno, a9, p_max, q_max).map(|m| m.residual)
}

/// Aggregate goodness-of-fit of a candidate Planet Nine semi-major axis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct A9Score {
    /// Candidate Planet Nine semi-major axis (AU).
    pub a9_au: f64,
    /// Candidate Planet Nine period (years).
    pub period_yr: f64,
    /// Mean over the constraining ETNOs of the distance of each period ratio to
    /// its nearest small-integer ratio. Lower is a better resonant fit.
    pub residual: f64,
    /// Number of ETNOs that were interior to a9 (and so contributed).
    pub n_constraints: usize,
}

/// Score a single candidate a9 against a set of ETNO semi-major axes: the mean
/// per-object residual to the nearest interior p:q (p ≤ `p_max`, q ≤ `q_max`)
/// resonance. Only ETNOs that are interior to a9 contribute; a candidate
/// constrained by *every* member of the set is scored over all of them.
///
/// `requires_all` enforces that every ETNO in the set is interior to a9 (the
/// Malhotra setup, where Planet Nine sits beyond all of them). When set and
/// the condition fails, the residual is `+∞` so such candidates lose the scan.
pub fn score_a9(a9: f64, a_etnos: &[f64], p_max: u32, q_max: u32, requires_all: bool) -> A9Score {
    let mut residual = 0.0;
    let mut n = 0usize;
    for &a in a_etnos {
        if let Some(r) = resonance_residual(a, a9, p_max, q_max) {
            residual += r;
            n += 1;
        }
    }
    let usable = n > 0 && (!requires_all || n == a_etnos.len());
    A9Score {
        a9_au: a9,
        period_yr: orbital_period_years(a9),
        residual: if usable {
            residual / n as f64
        } else {
            f64::INFINITY
        },
        n_constraints: n,
    }
}

/// Scan a9 over `[a_min, a_max]` AU in `steps` uniform steps, scoring each
/// candidate against the ETNO set with the Malhotra constraint that Planet
/// Nine lie beyond every ETNO. Returns the full curve (a9 ascending).
pub fn scan_a9(
    a_etnos: &[f64],
    a_min: f64,
    a_max: f64,
    steps: usize,
    p_max: u32,
    q_max: u32,
) -> Vec<A9Score> {
    assert!(steps >= 2 && a_max > a_min);
    (0..steps)
        .map(|k| {
            let a9 = a_min + (a_max - a_min) * k as f64 / (steps - 1) as f64;
            score_a9(a9, a_etnos, p_max, q_max, true)
        })
        .collect()
}

/// Best-fit Planet Nine semi-major axis: the minimum-residual a9 from a scan
/// over `[a_min, a_max]`, requiring all ETNOs to be interior. The constraining
/// set should be the distant clustered ETNOs that drive the localization.
pub fn best_fit_a9(
    a_etnos: &[f64],
    a_min: f64,
    a_max: f64,
    steps: usize,
    p_max: u32,
    q_max: u32,
) -> A9Score {
    scan_a9(a_etnos, a_min, a_max, steps, p_max, q_max)
        .into_iter()
        .filter(|s| s.residual.is_finite())
        .min_by(|a, b| a.residual.partial_cmp(&b.residual).unwrap())
        .expect("scan must contain a candidate beyond all ETNOs")
}

/// The longest-period clustered ETNOs from `BROWN_2017_SAMPLE`, sorted by
/// descending semi-major axis and truncated to the `n` most distant.
pub fn longest_period_etnos(n: usize) -> Vec<Etno> {
    let mut v: Vec<Etno> = BROWN_2017_SAMPLE.to_vec();
    v.sort_by(|x, y| y.a.partial_cmp(&x.a).unwrap());
    v.truncate(n);
    v
}

/// The four ETNOs Malhotra, Volk & Wang (2016) use to localize Planet Nine:
/// Sedna, 2012 VP113, 2004 VN112 and 2010 GB174 — the longest-period
/// clustered objects known at the time. Looked up by name from the vetted
/// `BROWN_2017_SAMPLE`. Panics if a name is missing (a data-table change that
/// should fail loudly rather than silently dropping a constraint).
pub fn malhotra_constraining_set() -> Vec<Etno> {
    ["Sedna", "2012 VP113", "2004 VN112", "2010 GB174"]
        .iter()
        .map(|name| {
            *BROWN_2017_SAMPLE
                .iter()
                .find(|e| e.name == *name)
                .unwrap_or_else(|| panic!("ETNO {name} missing from BROWN_2017_SAMPLE"))
        })
        .collect()
}

fn gcd(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = b;
        b = a % b;
        a = t;
    }
    a
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use p9_core::analysis::resonance::resonance_semi_major_axis;
    use p9_core::data::etno::semi_major_axes;

    // Malhotra et al.'s four longest-period clustered ETNOs.
    fn distant_set() -> Vec<f64> {
        malhotra_constraining_set().iter().map(|e| e.a).collect()
    }

    #[test]
    fn test_kepler_period_round_trips_published_values() {
        // 665 AU -> ~17,117 yr (Malhotra et al.).
        let p = orbital_period_years(MALHOTRA_2016_A9_AU);
        assert!(
            (p - MALHOTRA_2016_PERIOD_YR).abs() < 100.0,
            "P(665 AU) = {p:.0} yr, expected ~{MALHOTRA_2016_PERIOD_YR}"
        );
        // And the inverse: the published period implies the published a.
        let a = MALHOTRA_2016_PERIOD_YR.powf(2.0 / 3.0);
        assert_relative_eq!(a, MALHOTRA_2016_A9_AU, epsilon = 1.0);
    }

    #[test]
    fn test_resonance_semi_major_axis_round_trips_period_arithmetic() {
        // Interior p:q resonance: a_ETNO = a9 (q/p)^(2/3). Equivalently the
        // exterior helper resonance_semi_major_axis(p,q,a) places the OUTER
        // body p:q beyond an inner one; we use it to cross-check that
        // a9 = a_ETNO (p/q)^(2/3) reproduces the period-ratio inversion to
        // machine precision.
        let a9 = MALHOTRA_2016_A9_AU;
        for (p, q) in [(3u32, 2u32), (2, 1), (3, 1), (5, 2), (4, 1), (7, 2)] {
            // ETNO location for P9 at a9 in the interior p:q resonance.
            let a_etno = a9 * (q as f64 / p as f64).powf(2.0 / 3.0);
            // p9-core's exterior helper run the other direction recovers a9.
            assert_relative_eq!(resonance_semi_major_axis(p, q, a_etno), a9, epsilon = 1e-9);
            // The period ratio P9/P_ETNO at that location is exactly p/q.
            assert_relative_eq!(
                period_ratio(a9, a_etno),
                p as f64 / q as f64,
                epsilon = 1e-9
            );
            // The detector recovers (p, q) for a perfectly resonant object.
            let m = nearest_resonance(a_etno, a9, 9, 2).unwrap();
            assert_eq!((m.p, m.q), (p, q), "recovered {:?}", (m.p, m.q));
            assert!(m.residual < 1e-9);
        }
    }

    #[test]
    fn test_period_ratio_and_exterior_guard() {
        let a9 = 660.0;
        // ETNO interior to P9 -> ratio > 1.
        assert!(period_ratio(a9, 350.0) > 1.0);
        // Object beyond P9 -> no interior-resonance match.
        assert!(nearest_resonance(900.0, a9, 9, 2).is_none());
    }

    #[test]
    fn test_best_fit_a9_lands_near_published() {
        // Scan 400-900 AU restricted to N:1 / N:2 resonances (Malhotra's case);
        // the minimum-residual a9 should land within 10% of the published
        // ~665 AU and have a period near ~17,000 yr.
        let best = best_fit_a9(
            &distant_set(),
            400.0,
            900.0,
            5001,
            DEFAULT_P_MAX,
            DEFAULT_Q_MAX,
        );
        let published = MALHOTRA_2016_A9_AU;
        let rel = (best.a9_au - published).abs() / published;
        assert!(
            rel < 0.10,
            "best-fit a9 = {:.1} AU (residual {:.4}), published {published} AU, rel = {rel:.3}",
            best.a9_au,
            best.residual
        );
        // Period of the best-fit a9 near ~17,000 yr (within 15%).
        let rel_p = (best.period_yr - MALHOTRA_2016_PERIOD_YR).abs() / MALHOTRA_2016_PERIOD_YR;
        assert!(
            rel_p < 0.15,
            "best-fit period = {:.0} yr, published ~{MALHOTRA_2016_PERIOD_YR} yr, rel = {rel_p:.3}",
            best.period_yr
        );
        // Every ETNO in the distant set constrained the fit.
        assert_eq!(best.n_constraints, distant_set().len());
    }

    #[test]
    fn test_best_fit_consistent_with_millholland() {
        // The same scan also sits within 10% of Millholland & Laughlin's 654 AU
        // (the two papers' estimates differ by less than the scan tolerance).
        let best = best_fit_a9(
            &distant_set(),
            400.0,
            900.0,
            5001,
            DEFAULT_P_MAX,
            DEFAULT_Q_MAX,
        );
        let rel = (best.a9_au - MILLHOLLAND_2016_A9_AU).abs() / MILLHOLLAND_2016_A9_AU;
        assert!(
            rel < 0.10,
            "best-fit {:.1} vs Millholland 654 AU",
            best.a9_au
        );
    }

    #[test]
    fn test_distant_set_lands_on_low_order_resonances() {
        // At the best-fit a9 the distant ETNOs sit at simple p:q ratios — the
        // qualitative result of Malhotra et al. (3:2, 2:1, 5:2, 3:1, 4:1, ...).
        let best = best_fit_a9(
            &distant_set(),
            400.0,
            900.0,
            5001,
            DEFAULT_P_MAX,
            DEFAULT_Q_MAX,
        );
        for e in malhotra_constraining_set() {
            let m = nearest_resonance(e.a, best.a9_au, DEFAULT_P_MAX, DEFAULT_Q_MAX).unwrap();
            // q is 1 or 2 (N:1 / N:2), and the match is reasonably close.
            assert!(m.q <= 2, "{}: q = {}", e.name, m.q);
            assert!(
                m.residual < 0.2,
                "{} at {}:{} residual {:.3}",
                e.name,
                m.p,
                m.q,
                m.residual
            );
        }
    }

    #[test]
    fn test_residual_curve_shows_bailey_degeneracy() {
        // HONEST counterpoint to Bailey, Brown & Batygin (2018): the resonance
        // localization does have a minimum near 660 AU, but the residual curve
        // is riddled with near-equally-good a9 values, so the "prediction" is
        // far weaker than a single sharp peak would suggest.
        let set = distant_set();
        let curve = scan_a9(&set, 400.0, 900.0, 5001, DEFAULT_P_MAX, DEFAULT_Q_MAX);
        let finite: Vec<f64> = curve
            .iter()
            .filter(|s| s.residual.is_finite())
            .map(|s| s.residual)
            .collect();
        let min_r = finite.iter().cloned().fold(f64::INFINITY, f64::min);
        let n_good = finite
            .iter()
            .filter(|&&r| r < 1.5 * min_r.max(1e-6))
            .count();
        assert!(
            n_good > 20,
            "expected many near-degenerate minima (Bailey+ 2018), got {n_good} of {}",
            finite.len()
        );
    }

    #[test]
    fn test_constraining_set_membership() {
        // Malhotra et al.'s four objects, drawn from the vetted sample.
        let names: Vec<&str> = malhotra_constraining_set().iter().map(|e| e.name).collect();
        assert_eq!(names.len(), 4);
        for want in ["Sedna", "2012 VP113", "2004 VN112", "2010 GB174"] {
            assert!(names.contains(&want), "missing {want}");
        }
        let all = semi_major_axes();
        for e in malhotra_constraining_set() {
            assert!(e.a > 250.0, "{}: a = {}", e.name, e.a);
            assert!(all.contains(&e.a));
        }
        // longest_period_etnos still returns the strict high-a tail.
        assert_eq!(longest_period_etnos(1)[0].name, "Sedna");
    }
}
