//! Reproduction of Schwamb et al. (2023), "Tuning the Legacy Survey of Space
//! and Time (LSST) Observing Strategy for Solar System Science"
//! (arXiv:2303.02355).
//!
//! HEADLINE: LSST's cadence and footprint choices strongly affect the
//! discovery of distant solar-system bodies, including Planet Nine. Depth per
//! visit, the number of visits required to *link* a slow mover into a
//! tracklet/orbit, and the sky coverage all drive the completeness; a
//! well-tuned strategy maximises the discoverable fraction of slow distant
//! movers.
//!
//! WHAT THIS CRATE COMPUTES (no hard-coded answers): it draws a seeded
//! reference Planet Nine population (mass, distance, sky position,
//! reflected-light `V`) — reusing the `p9-2021-orbit` posterior emulator and
//! `p9_core::analysis::photometry` — and runs it through a real LSST detection
//! model built on the reused `p9-survey` Rubin/LSST footprint preset and the
//! shared survey depths. It then shows that the **discoverable fraction**
//!
//!   * is a probability in `(0, 1)`,
//!   * *decreases monotonically* when more visits are required for linking,
//!   * *decreases monotonically* as the footprint declination extent shrinks,
//!   * *increases* when the ten-year co-add depth replaces the single-visit
//!     depth.
//!
//! The published LSST parameters are kept as labelled reference constants in
//! [`strategy::published`]; residuals (where our reduced model cannot match a
//! full SSP/`MAF` simulation) are documented in `REPRODUCTION_NOTES.md`.
//!
//! Modules:
//!   * [`strategy`] — the LSST observing-strategy model and its knobs.
//!   * [`discoverable`] — the discoverable fraction of a P9 population.

pub mod discoverable;
pub mod strategy;

pub use discoverable::{
    dec_and_galactic_lat_deg, discoverable_fraction, discovery_probability, DiscoverableResult,
    P9Sample,
};
pub use strategy::{binomial_survival, LsstStrategy};

#[cfg(test)]
mod headline_tests {
    use super::*;
    use p9_core::data::reference_population::generate_reference_population;
    use p9_core::types::OrbitalElements;
    use rand::SeedableRng;

    /// Build a seeded reference Planet Nine population (sky position +
    /// reflected-light V) from the p9-2021-orbit posterior emulator.
    fn reference_population(n: usize, seed: u64) -> Vec<P9Sample> {
        let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
        generate_reference_population(n, &mut rng)
            .into_iter()
            .map(|p| P9Sample {
                elements: OrbitalElements {
                    a: p.a,
                    e: p.e,
                    i: p.i,
                    omega: p.omega,
                    omega_big: p.omega_big,
                    mean_anomaly: p.mean_anomaly,
                },
                v_magnitude: p.v_magnitude,
            })
            .collect()
    }

    /// The discoverable fraction is a real probability in (0, 1) for the
    /// baseline LSST strategy on a seeded reference population.
    #[test]
    fn baseline_discoverable_fraction_in_open_unit_interval() {
        let pop = reference_population(4000, 2023);
        let r = discoverable_fraction(&LsstStrategy::baseline(), &pop);
        assert!(
            r.fraction > 0.0 && r.fraction < 1.0,
            "baseline discoverable fraction = {:.4}",
            r.fraction
        );
        // Sanity: expected count is consistent with the fraction.
        assert!((r.expected_discovered - r.fraction * pop.len() as f64).abs() < 1e-9);
        // The footprint-only fraction is an upper bound on the discoverable
        // fraction (linking can only remove objects, not add them).
        assert!(r.fraction <= r.fraction_in_footprint + 1e-12);
        // The southern Rubin footprint with coverage 0.85 catches a
        // substantial but not unit fraction.
        assert!(r.fraction_in_footprint > 0.1 && r.fraction_in_footprint < 0.9);
    }

    /// Requiring more independent-night detections for linking monotonically
    /// reduces the discoverable fraction (the cadence/linking lever). This
    /// holds on the full reference population; the magnitude of the drop is
    /// modest there because — the honest finding — the nominal-orbit P9
    /// population is bright (V ≲ 23) relative to the single-visit depth, so
    /// most objects clear the per-visit threshold on nearly every night and
    /// linking is not the binding constraint (the footprint is). The lever
    /// bites hardest on the faint distant tail, exercised separately below.
    #[test]
    fn more_visits_for_linking_monotonically_reduces_discovery() {
        let pop = reference_population(4000, 7);
        let mut prev = f64::INFINITY;
        for n in [1u32, 2, 3, 5, 8, 12] {
            let strat = LsstStrategy::baseline().with_visits_for_linking(n);
            let f = discoverable_fraction(&strat, &pop).fraction;
            assert!(
                f <= prev + 1e-12,
                "discoverable fraction rose at n = {n}: {f:.4} > {prev:.4}"
            );
            prev = f;
        }
    }

    /// The linking lever has real teeth on the faint, distant movers — the
    /// regime Schwamb et al. (2023) care about. At a shallower single-visit
    /// depth (so a substantial fraction of the population sits near the
    /// per-visit limit), requiring more independent-night detections drops the
    /// discoverable fraction appreciably and monotonically.
    #[test]
    fn linking_lever_bites_hard_near_the_depth_limit() {
        let pop = reference_population(4000, 7);
        // A shallow per-visit cadence (e.g. time split across many filters):
        // depth ~ 21 puts much of the V≈19–23 population near the limit.
        let shallow = |n: u32| {
            LsstStrategy::baseline()
                .with_single_visit_depth(21.0)
                .with_visits_for_linking(n)
        };
        let mut prev = f64::INFINITY;
        let mut at1 = 0.0;
        let mut at12 = 0.0;
        for n in [1u32, 2, 3, 5, 8, 12] {
            let f = discoverable_fraction(&shallow(n), &pop).fraction;
            assert!(f <= prev + 1e-12, "rose at n = {n}: {f:.4} > {prev:.4}");
            if n == 1 {
                at1 = f;
            }
            if n == 12 {
                at12 = f;
            }
            prev = f;
        }
        assert!(
            at1 - at12 > 0.05,
            "linking threshold barely mattered near the limit: {at1:.4} -> {at12:.4}"
        );
    }

    /// Shrinking the footprint (lowering the northern declination limit)
    /// monotonically reduces the discoverable fraction (the footprint lever).
    #[test]
    fn shrinking_footprint_monotonically_reduces_discovery() {
        let pop = reference_population(4000, 99);
        // Sweep dec_max from the published +12° down toward the deep south.
        let mut prev = -1.0;
        let mut full = 0.0;
        let mut tiny = 0.0;
        for &dec_max in &[-60.0, -45.0, -30.0, -10.0, 12.0] {
            let strat = LsstStrategy::baseline().with_dec_max(dec_max);
            let f = discoverable_fraction(&strat, &pop).fraction;
            assert!(
                f >= prev - 1e-12,
                "fraction should rise as dec_max grows; dropped at {dec_max}: {f:.4} < {prev:.4}"
            );
            if dec_max == 12.0 {
                full = f;
            }
            if dec_max == -60.0 {
                tiny = f;
            }
            prev = f;
        }
        assert!(
            full - tiny > 0.02,
            "footprint extent barely mattered: tiny {tiny:.4} vs full {full:.4}"
        );
    }

    /// Deepening from the single-visit depth to the ten-year co-add can only
    /// raise the discoverable fraction (the depth lever): the stack links
    /// fainter, more distant solutions and never loses one. On the bright
    /// nominal-orbit population the gain is tiny (every object already clears
    /// the 24.5 single-visit depth — the honest finding that depth is not the
    /// binding constraint there); the gain is real and large when the
    /// per-visit cadence is shallow, exercised below.
    #[test]
    fn ten_year_stack_depth_never_lowers_discovery() {
        let pop = reference_population(4000, 314);
        let baseline = discoverable_fraction(&LsstStrategy::baseline(), &pop).fraction;
        let deep =
            discoverable_fraction(&LsstStrategy::baseline().with_stack_depth(), &pop).fraction;
        assert!(
            deep >= baseline - 1e-12,
            "stack depth lowered discovery: deep {deep:.4} vs single-visit {baseline:.4}"
        );
        assert!(deep > 0.0 && deep < 1.0);
    }

    /// When the per-visit cadence is shallow (so much of the population is
    /// near the single-visit limit), switching discovery to the ten-year
    /// co-add depth raises the discoverable fraction appreciably — the
    /// shift-and-stack depth lever Schwamb et al. highlight for distant
    /// movers.
    #[test]
    fn stack_depth_increases_discovery_near_the_limit() {
        let pop = reference_population(4000, 314);
        let shallow = LsstStrategy::baseline().with_single_visit_depth(21.0);
        let single = discoverable_fraction(&shallow, &pop).fraction;
        let deep = discoverable_fraction(&shallow.with_stack_depth(), &pop).fraction;
        assert!(
            deep - single > 0.05,
            "stack depth barely helped near the limit: {single:.4} -> {deep:.4}"
        );
    }

    /// The discoverable fraction is reproducible across seeds to within Monte
    /// Carlo scatter (the model is deterministic given a population; this pins
    /// the seed-to-seed stability of the seeded draw).
    #[test]
    fn discoverable_fraction_is_seed_stable() {
        let a = discoverable_fraction(&LsstStrategy::baseline(), &reference_population(4000, 1))
            .fraction;
        let b = discoverable_fraction(&LsstStrategy::baseline(), &reference_population(4000, 2))
            .fraction;
        assert!(
            (a - b).abs() < 0.05,
            "seed scatter too large: {a:.4} vs {b:.4}"
        );
    }
}
