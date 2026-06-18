//! Locations of Neptune's distant exterior mean-motion resonances and the
//! clustered ETNOs that sit near them.
//!
//! The companion 2021 work (arXiv:2105.01065) maps Neptune's distant n:1 and
//! n:2 exterior resonances across the ~100–800 AU range and asks which
//! clustered ETNOs lie near them. An exterior n:m resonance with Neptune (the
//! particle completes m orbits per n Neptune orbits, n > m) sits at
//!
//!   a_res = a_N (n/m)^{2/3}
//!
//! which is exactly `p9_core::analysis::resonance::resonance_semi_major_axis`.
//! This module only catalogs the relevant (n, m) and finds the nearest
//! resonance for each object — it adds no new resonance physics.

use p9_core::analysis::resonance::resonance_semi_major_axis;
use p9_core::constants::A_NEPTUNE_AU;
use p9_core::data::etno::{Etno, BROWN_2017_SAMPLE};
use p9_core::units::{au, Length};

/// A located Neptune exterior resonance.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Resonance {
    /// Period ratio numerator (Neptune orbits).
    pub n: u32,
    /// Period ratio denominator (particle orbits).
    pub m: u32,
    /// Semi-major axis of the resonance (AU).
    pub a_au: f64,
}

impl Resonance {
    fn new(n: u32, m: u32) -> Self {
        Self {
            n,
            m,
            a_au: resonance_semi_major_axis(n, m, A_NEPTUNE_AU),
        }
    }

    /// "n:m" label, e.g. "5:1".
    pub fn label(&self) -> String {
        format!("{}:{}", self.n, self.m)
    }

    /// Resonance semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }
}

/// Neptune's exterior n:1 resonances for n = 2..=12 — the headline chain
/// (a ≈ 48–158 AU) anchoring the lower half of the distant range.
pub fn neptune_n1_resonances() -> Vec<Resonance> {
    (2..=12).map(|n| Resonance::new(n, 1)).collect()
}

/// A selection of Neptune's exterior n:2 resonances (odd n so the ratio is in
/// lowest terms) spanning the distant range.
pub fn neptune_n2_resonances() -> Vec<Resonance> {
    [5u32, 7, 9, 11, 13, 15, 17, 19, 21, 23]
        .iter()
        .map(|&n| Resonance::new(n, 2))
        .collect()
}

/// Neptune's exterior n:1 resonances out to a_au ≈ 800 AU. The n:1 chain runs
/// to n ≈ 130 at 800 AU and the spacing shrinks to a few AU past ~400 AU, so
/// every clustered ETNO (a ≈ 228–506 AU) sits within a few AU of one.
pub fn distant_n1_resonances() -> Vec<Resonance> {
    (2u32..=130)
        .map(|n| Resonance::new(n, 1))
        .take_while(|r| r.a_au <= 800.0)
        .collect()
}

/// All cataloged distant resonances (n:1 to ~800 AU plus the n:2 selection)
/// sorted by semi-major axis. Used for nearest-resonance matching.
pub fn distant_resonances() -> Vec<Resonance> {
    let mut all = distant_n1_resonances();
    all.extend(neptune_n2_resonances());
    all.sort_by(|x, y| x.a_au.partial_cmp(&y.a_au).unwrap());
    all
}

/// The resonance from `catalog` whose semi-major axis is closest to `a_au`,
/// with the absolute offset in AU.
pub fn nearest_resonance(a_au: f64, catalog: &[Resonance]) -> (Resonance, f64) {
    catalog
        .iter()
        .map(|&r| (r, (r.a_au - a_au).abs()))
        .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap())
        .expect("non-empty catalog")
}

/// An ETNO paired with its nearest distant Neptune resonance.
#[derive(Debug, Clone)]
pub struct EtnoResonance {
    pub name: &'static str,
    pub a_au: f64,
    pub resonance: Resonance,
    /// Signed offset a_etno − a_res (AU); negative means interior to the
    /// resonance.
    pub offset_au: f64,
}

impl EtnoResonance {
    /// ETNO semi-major axis as a typed [`Length`].
    pub fn semi_major_axis(&self) -> Length {
        au(self.a_au)
    }

    /// Signed offset from the nearest resonance as a typed [`Length`].
    pub fn offset(&self) -> Length {
        au(self.offset_au)
    }
}

/// Match every member of the Brown (2017) clustered sample to its nearest
/// distant Neptune resonance (n:1 ∪ n:2).
pub fn match_sample_to_resonances() -> Vec<EtnoResonance> {
    let catalog = distant_resonances();
    BROWN_2017_SAMPLE
        .iter()
        .map(|etno: &Etno| {
            let (res, _) = nearest_resonance(etno.a, &catalog);
            EtnoResonance {
                name: etno.name,
                a_au: etno.a,
                resonance: res,
                offset_au: etno.a - res.a_au,
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_n1_resonances_match_kepler_relation() {
        // The headline pin: each n:1 a-value matches the analytic Kepler
        // relation a = a_N · n^(2/3) to within 1e-6 AU, cross-checking the
        // shared `resonance_semi_major_axis`.
        for r in neptune_n1_resonances() {
            assert_eq!(r.m, 1);
            let kepler = A_NEPTUNE_AU * (r.n as f64).powf(2.0 / 3.0);
            assert_relative_eq!(r.a_au, kepler, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_n2_resonances_match_kepler_relation() {
        for r in neptune_n2_resonances() {
            assert_eq!(r.m, 2);
            let kepler = A_NEPTUNE_AU * (r.n as f64 / 2.0).powf(2.0 / 3.0);
            assert_relative_eq!(r.a_au, kepler, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_n1_range_spans_distant_region() {
        let res = neptune_n1_resonances();
        // 2:1 just beyond Neptune, 12:1 out near ~480 AU — the n=2..12 chain
        // covers the lower half of the 100–800 AU band.
        let a2 = res.first().unwrap().a_au;
        let a12 = res.last().unwrap().a_au;
        assert!((a2 - 47.8).abs() < 0.5, "2:1 at {a2:.1} AU");
        assert!(a12 > 150.0 && a12 < 200.0, "12:1 at {a12:.1} AU");
    }

    #[test]
    fn test_resonances_sorted_and_nonempty() {
        let all = distant_resonances();
        assert!(!all.is_empty());
        for w in all.windows(2) {
            assert!(w[0].a_au <= w[1].a_au, "catalog must be a-sorted");
        }
    }

    #[test]
    fn test_typed_accessors_match_f64() {
        use uom::si::length::astronomical_unit;
        let matches = match_sample_to_resonances();
        let m = &matches[0];
        assert_relative_eq!(
            m.semi_major_axis().get::<astronomical_unit>(),
            m.a_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            m.offset().get::<astronomical_unit>(),
            m.offset_au,
            epsilon = 1e-9
        );
        assert_relative_eq!(
            m.resonance.semi_major_axis().get::<astronomical_unit>(),
            m.resonance.a_au,
            epsilon = 1e-9
        );
    }

    #[test]
    fn test_nearest_resonance_picks_closest() {
        let catalog = distant_resonances();
        // Build a point exactly on the 9:1 resonance; nearest must be it.
        let nine = Resonance::new(9, 1);
        let (got, off) = nearest_resonance(nine.a_au, &catalog);
        assert_eq!(got, nine);
        assert!(off < 1e-9, "offset {off}");
    }

    #[test]
    fn test_sample_matches_have_modest_offsets() {
        // Every clustered ETNO should land within a fraction of the local
        // resonance spacing of *some* distant resonance — the offsets are
        // tens of AU, not hundreds, across the 228–506 AU sample.
        let matches = match_sample_to_resonances();
        assert_eq!(matches.len(), BROWN_2017_SAMPLE.len());
        for m in &matches {
            assert!(
                m.offset_au.abs() < 60.0,
                "{} (a={:.0}) nearest {} at offset {:.1} AU",
                m.name,
                m.a_au,
                m.resonance.label(),
                m.offset_au
            );
        }
        // Sedna (a = 506 AU) is the most distant; its nearest distant
        // resonance must itself be far out (a_res > 300 AU).
        let sedna = matches.iter().find(|m| m.name == "Sedna").unwrap();
        assert!(
            sedna.resonance.a_au > 300.0,
            "Sedna nearest resonance {} at {:.0} AU",
            sedna.resonance.label(),
            sedna.resonance.a_au
        );
    }
}
