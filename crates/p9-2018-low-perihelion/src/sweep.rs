//! Sweep e9 at fixed a9, mapping each P9 configuration to its perihelion
//! q9 = a9(1 − e9), and tabulate the confinement metric as q9 decreases.
//!
//! This is the crate headline: confinement of a test ETNO persists — the
//! apsidal libration island stays finite (and in fact tightens) and the
//! apsidal well deepens — as the Planet Nine perihelion is lowered well below
//! the canonical Brown & Batygin q9 ≈ 280 AU, down to the documented low-q9
//! floor. So a low-perihelion (high-e9) Planet Nine still clusters the ETNOs,
//! widening the dynamically-allowed (a9, e9) region.

use crate::confinement::{analyze, Confinement};
use crate::reference::{E9_SWEEP_HI, E9_SWEEP_LO};

/// One row of the q9 sweep.
#[derive(Debug, Clone, Copy)]
pub struct SweepRow {
    pub e9: f64,
    pub q9_au: f64,
    pub confinement: Confinement,
}

/// Sweep e9 ∈ [E9_SWEEP_LO, E9_SWEEP_HI] at fixed `a9` (so q9 = a9(1−e9)
/// DECREASES across the sweep), at `n` points, for a test ETNO of semi-major
/// axis `a` and a P9 mass of `mass_earth`. Rows are ordered by increasing e9,
/// i.e. DECREASING q9.
pub fn sweep_q9(a: f64, a9: f64, mass_earth: f64, n: usize) -> Vec<SweepRow> {
    assert!(n >= 2, "need at least two sweep points");
    (0..n)
        .map(|k| {
            let e9 = E9_SWEEP_LO + (E9_SWEEP_HI - E9_SWEEP_LO) * (k as f64) / (n as f64 - 1.0);
            let confinement = analyze(a, a9, e9, mass_earth);
            SweepRow {
                e9,
                q9_au: a9 * (1.0 - e9),
                confinement,
            }
        })
        .collect()
}

/// The lowest perihelion (AU) in the sweep at which the configuration still
/// confines (finite sub-π libration island AND a real apsidal well).
pub fn lowest_confining_q9(rows: &[SweepRow]) -> Option<f64> {
    rows.iter()
        .filter(|r| r.confinement.is_confined())
        .map(|r| r.q9_au)
        .fold(None, |acc, q| Some(acc.map_or(q, |a: f64| a.min(q))))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::reference::{BB_A9_AU, BB_Q9_AU, LOW_Q9_MIN_AU, P9_MASS_EARTH, TEST_ETNO_A_AU};

    fn run() -> Vec<SweepRow> {
        sweep_q9(TEST_ETNO_A_AU, BB_A9_AU, P9_MASS_EARTH, 7)
    }

    #[test]
    fn sweep_spans_low_perihelion_range() {
        let rows = run();
        // q9 decreases monotonically from the canonical value toward the floor.
        assert!((rows.first().unwrap().q9_au - BB_Q9_AU).abs() < 1e-9);
        assert!((rows.last().unwrap().q9_au - LOW_Q9_MIN_AU).abs() < 1e-9);
        for w in rows.windows(2) {
            assert!(w[1].q9_au < w[0].q9_au, "q9 must decrease across the sweep");
        }
    }

    #[test]
    fn confinement_persists_at_low_perihelion() {
        // HEADLINE. The confinement metric is non-zero and apsidal libration
        // persists across the whole low-q9 range — including q9 far below the
        // canonical 280 AU. Every swept configuration confines.
        let rows = run();
        for r in &rows {
            assert!(
                r.confinement.is_confined(),
                "lost confinement at q9 = {:.1} AU (e9 = {:.2}): {:?}",
                r.q9_au,
                r.e9,
                r.confinement
            );
            assert!(r.confinement.libration_half_width > 0.0);
            assert!(r.confinement.torque_amplitude > 0.0);
        }
        // The lowest confining perihelion reaches the documented floor.
        let qmin = lowest_confining_q9(&rows).expect("at least one confining row");
        assert!(
            qmin <= LOW_Q9_MIN_AU + 1e-9,
            "expected confinement down to {LOW_Q9_MIN_AU} AU, only reached {qmin}"
        );
        assert!(
            qmin < BB_Q9_AU,
            "low-perihelion confinement must extend below the canonical 280 AU"
        );
    }

    #[test]
    fn well_depth_increases_monotonically_as_perihelion_drops() {
        // The cleanest confinement-strength metric. Lowering q9 (raising e9) at
        // fixed a9 brings P9 closer to the ETNOs at perihelion and deepens its
        // eccentric (octupole+) apsidal well: the dimensionless well depth must
        // rise smoothly and monotonically as q9 falls. This is the
        // "varies smoothly/monotonically with q9/e9" gate.
        let rows = run();
        for w in rows.windows(2) {
            // rows ordered by decreasing q9, so well depth should be increasing
            assert!(
                w[1].confinement.well_depth >= w[0].confinement.well_depth - 1e-6,
                "well depth fell as q9 dropped: {:.4} (q9={:.0}) -> {:.4} (q9={:.0})",
                w[0].confinement.well_depth,
                w[0].q9_au,
                w[1].confinement.well_depth,
                w[1].q9_au
            );
        }
        // The lowest-q9 well is meaningfully deeper than the canonical one.
        let first = rows.first().unwrap().confinement.well_depth;
        let last = rows.last().unwrap().confinement.well_depth;
        assert!(
            last > first + 0.05,
            "low-q9 well depth {last:.4} should exceed canonical {first:.4} by a clear margin"
        );
    }

    #[test]
    fn libration_island_tightens_as_perihelion_drops() {
        // As the P9 perihelion falls the apsidal libration island gets TIGHTER
        // (smaller half-width): the ETNO apses are confined to a narrower wedge
        // around the favoured apse. Monotone decrease in the half-width.
        let rows = run();
        for w in rows.windows(2) {
            assert!(
                w[1].confinement.libration_half_width
                    <= w[0].confinement.libration_half_width + 1e-6,
                "island widened as q9 dropped: {:.4} (q9={:.0}) -> {:.4} (q9={:.0})",
                w[0].confinement.libration_half_width,
                w[0].q9_au,
                w[1].confinement.libration_half_width,
                w[1].q9_au
            );
        }
        // Still a finite (sub-π) island at the lowest perihelion: confinement
        // persists, it does not dissolve into circulation.
        assert!(rows.last().unwrap().confinement.libration_half_width < std::f64::consts::PI);
    }

    #[test]
    fn favored_apse_is_consistent_across_sweep() {
        // The single-ring secular problem confines toward the SAME apse across
        // the whole low-q9 range — the mechanism is coherent, not flickering
        // between aligned/anti-aligned as q9 changes. (See the module-level
        // sign caveat: this bare test particle favours the aligned apse; the
        // OBSERVED anti-alignment needs P9 precession + the giants' field.)
        let rows = run();
        let first = rows[0].confinement.favored_apse;
        for r in &rows {
            assert_eq!(
                r.confinement.favored_apse, first,
                "favoured apse flipped at q9 = {:.0}",
                r.q9_au
            );
        }
    }
}
